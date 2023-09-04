import matplotlib.pyplot as plt
from pandas.io.sql import BaseEngine

from DualNumberModule import DualNumber as dn
from DualNumberModule import DualMatrix as dm
import ModuleOfRotationMatrices as rm
import numpy as np
from numpy import linalg as LA
from ModelAnimation import Begin
from ModelAnimation import Animation
#Введем начальные параметры по Д-Х для первого цикла
D = np.array([0.315, 0, 0.310, 0, 0.300, 0, 0.205])
aa = np.array([0.05, 0, 0, 0, 0, 0, 0])
alpha = np.array([( np.pi ) / 2, -(np.pi)/2, (np.pi)/2, -(np.pi)/2,
         (np.pi)/2, -(np.pi)/2,0])
#На этом моменте нужно реализовать решение задачи МНК
# Theta = np.array([-1.2444, 1.11285, 0.547214, -1.41126, 1.09214, -0.975915, -1.63531])
Theta = np.array([2.85644, -0.325443, 0.821165, 1.47851, -0.89284, 0.654015, 2.6179])

#Введем дополниительные вычислительные параметры для метода Ньютона , а именно матрицу ориентации Gd и
# построим прогнозную траекторию движения,по координатам x,y,z
Gd = np.array([[0, 0, 1],[-1, 0, 0],[0, -1, 0]])

#Для начала построим желаемую траекторию движения
coordinatesy = [0.2]
coordinatesz= [0.7]
radius = coordinatesy[0]
zhigh = coordinatesz[0]
xd = 0.3
angel = np.pi
#Разбиение
a=500
dangel = angel / a
Startangel = 0
while (Startangel < angel):
    Startangel = Startangel + dangel
    coordz = np.sin(Startangel) * radius + zhigh
    coordy = np.cos(Startangel) * radius
    coordinatesy.append(coordy)
    coordinatesz.append(coordz)

yd = coordinatesy[0]
zd = coordinatesz[0]
plt.plot(coordinatesy,coordinatesz)
plt.grid()

plt.xlabel('Y,m', size = 30)
plt.ylabel('Z,m', size = 30)
plt.title('Траектория движения в плоскости Y-Z', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)

ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()

n=len(coordinatesy)
a=500 #Разбиение на 200 точек
V = 0.1 # Модуль скорости
tt = [0]*(n)
Vxd = 0
Vy = [0]*(n-1)
Vz = [0]*(n-1)
for i in range(0,n-1):
    chisnormy = coordinatesy[i+1]-coordinatesy[i]
    chisnormz = coordinatesz[i + 1] - coordinatesz[i]
    chismass=[]
    chismass.append((chisnormy,chisnormz))
    T =LA.norm(chismass) / V
    tt[i+1]=tt[i]+T
    Vy[i] = chisnormy / T
    Vz[i] = chisnormz / T
omega = 0.1;

t =tt[:n-1]
deltaT = tt[n-1] / n
figure, axis = plt.subplots(2, 2)
axis[0, 0].plot(t, Vy)
axis[0, 0].set_title("X(t)")
axis[0, 0].set_ylabel('X,m')
axis[0, 0].set_xlabel('t,sec')

axis[0, 1].plot(t, Vz)
axis[0, 1].set_title("Y(t)")
axis[0, 1].set_ylabel('Y,m')
axis[0, 1].set_xlabel('t,sec')
plt.show()
ryy = yd
rx = xd
rzz = zd
Ryd = [ryy]
Rzd = [rzz]
Rxd = [rx]
for i in range(1,n-1):
    ryy = ryy + Vy[i] * deltaT
    rzz = rzz + Vz[i] * deltaT
    rxx = xd
    Ryd.append(ryy)
    Rzd.append(rzz)
    Rxd.append(rxx)
#Построение графиков
plt.plot(t, Rxd)
plt.grid()
plt.xlabel('t,sec', size = 30)
plt.ylabel('X,m', size = 30)
plt.title('X(t)', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)
ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()

plt.plot(t, Ryd)
plt.grid()
plt.xlabel('t,sec', size = 30)
plt.ylabel('Y,m', size = 30)
plt.title('Y(t)', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)
ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()

plt.plot(t, Rzd)
plt.grid()
plt.xlabel('t,sec', size = 30)
plt.ylabel('Z,m', size = 30)
plt.title('Z(t)', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)
ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()
q = []
q1 = []
q2 = []
q3 = []
q4 = []
q5 = []
q6 = []
q7 = []
Xnew = []
Ynew = []
Znew = []
q.append(Theta) #Массив обощенных координат

GG = []
Gabs = dm.eye(3)
    #Составление абсолютного верзора
for j in range(7):
    Rzm = rm.RMz(D[j], Theta[j])
    Rxm = rm.RMx(aa[j], alpha[j])
    Gper = Rzm * Rxm
    Gabs = Gabs*Gper
    GG.append(Gabs)
Matrix = dm.DTM(Gabs, False)
VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
x = VMatrix[2][1]
y = VMatrix[0][2]
z = VMatrix[1][0]
Xnew.append(x)
Ynew.append(y)
Znew.append(z)
#Реализация самого Алгоритма Ньютона
for i in range(1,n-1):
    Xprog = Rxd[i]
    Yprog = Ryd[i]
    Zprog = Rzd[i]
    # Вычислительный блок для решения прямой задачи кинематики
    #Сформируем скорости выходного звена

    GG = []
    Gabs = dm.eye(3)
    #Составление абсолютного верзора
    for j in range(7):
        Rzm = rm.RMz(D[j], Theta[j])
        Rxm = rm.RMx(aa[j], alpha[j])
        Gper = Rzm * Rxm
        Gabs = Gabs*Gper
        GG.append(Gabs)
    Matrix = dm.DTM(Gabs, False)
    VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
    x = VMatrix[2][1]
    y = VMatrix[0][2]
    z = VMatrix[1][0]
    #Окончание вычлительного блока для решения прямой задачи

    #Вычислительный блок для составления  матриц  якоби
    IdenRotVector = dm([[dn(0,0)],[dn(0,0)],[dn(1,0)]])
    Jac = []
    Jac.append(IdenRotVector)
    for j in range(len(GG)-1):
        StolbecJac = GG[j] * IdenRotVector
        Jac.append(StolbecJac)
    JacCon = dm.CrMaFrCol(Jac)
    Mat = dm.DTM(JacCon,False)
    ConMat = np.concatenate((np.array(Mat[0]),np.array(Mat[1])),axis=0)
    #Окончание вычислительного блока для составления матрицы якоби

    #Итерационный блок для решения обратной задачи методом Ньютона
    Gamma = []
    #Составление матрицы Gamma
    Gamma = (dm.DTM(Gabs, False))[0]
    #
    rAx = x
    rAy = y
    rAz = z
    Xprog = Rxd[i]
    Yprog = Ryd[i]
    Zprog = Rzd[i]
    dr = [Xprog - rAx, Yprog - rAy, Zprog - rAz]
    # Пересчитать 
    # dr = dr + ra* betta
    dG = Gd - Gamma
    dG = np.dot(dG, np.transpose(Gd))
    Betta = [dG[2, 1], dG[0, 2], dG[1, 0]]
    k=0
    while(((LA.norm(dr) > 0.001) or (LA.norm(Betta) > 0.001)) and k<10 ):
        k=k+1
        J0 = ConMat
        dU = np.concatenate((Betta, dr), axis=None)
        Theta = Theta - np.dot(LA.pinv(J0), -dU)

        GG = []
        Gabs = dm.eye(3)
        # Составление абсолютного верзора
        for j in range(7):
            Rzm = rm.RMz(D[j], Theta[j])
            Rxm = rm.RMx(aa[j], alpha[j])
            Gper = Rzm * Rxm
            Gabs = Gabs * Gper
            GG.append(Gabs)
        Matrix = dm.DTM(Gabs, False)
        VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
        x = VMatrix[2][1]
        y = VMatrix[0][2]
        z = VMatrix[1][0]
        # Окончание вычлительного блока для решения прямой задачи

        # Вычислительный блок для составления  матриц  якоби
        IdenRotVector = dm([[dn(0, 0)], [dn(0, 0)], [dn(1, 0)]])
        Jac = []
        Jac.append(IdenRotVector)
        for j in range(len(GG) - 1):
            StolbecJac = GG[j] * IdenRotVector
            Jac.append(StolbecJac)
        JacCon = dm.CrMaFrCol(Jac)
        Mat = dm.DTM(JacCon, False)
        ConMat = np.concatenate((np.array(Mat[0]), np.array(Mat[1])), axis=0)
        # Окончание вычислительного блока для составления матрицы якоби

        # Итерационный блок для решения обратной задачи методом Ньютона
        Gamma = []
        # Составление матрицы Gamma
        Gamma = (dm.DTM(Gabs, False))[0]
        #
        rAx = x
        rAy = y
        rAz = z
        Xprog = Rxd[i]
        Yprog = Ryd[i]
        Zprog = Rzd[i]
        dr = [Xprog - rAx, Yprog - rAy, Zprog - rAz]
        dG = Gd - Gamma
        dG = np.dot(dG, np.transpose(Gd))
        Betta = [dG[2, 1], dG[0, 2], dG[1, 0]]
    q.append(Theta)
    Xnew.append(x)
    Ynew.append(y)
    Znew.append(z)

for i in range(0,n-1):
    q1.append(q[i][0])
    q2.append(q[i][1])
    q3.append(q[i][2])
    q4.append(q[i][3])
    q5.append(q[i][4])
    q6.append(q[i][5])
    q7.append(q[i][6])
figure, axis = plt.subplots(3, 3)

axis[0, 0].plot(Ynew, Znew)
axis[0, 0].set_title("Траектория в Y-Z")
axis[0, 0].set_xlabel('Y,m')
axis[0, 0].set_ylabel('Z,m')

axis[0, 1].plot(t, q1)
axis[0, 1].set_title("q1(t)")


axis[1, 0].plot(t, q2)
axis[1, 0].set_title("q2(t)")


axis[1, 1].plot(t, q3)
axis[1, 1].set_title("q3(t)")

axis[2, 1].plot(t, q4)
axis[2, 1].set_title("q4(t)")

axis[1, 2].plot(t, q5)
axis[1, 2].set_title("q5(t)")

axis[2, 2].plot(t, q6)
axis[2, 2].set_title("q6(t)")

axis[0, 2].plot(t, q7)
axis[0, 2].set_title("q7(t)")
plt.show()

plt.plot(Ynew, Znew)
plt.grid()
plt.xlabel('Y,m', size = 30)
plt.ylabel('Z,m', size = 30)
plt.title('Отработанная траектория Y-Z', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)
ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()

plt.plot(t, q1)
plt.grid()
plt.xlabel('t,sec', size = 30)
plt.ylabel('q1,rad', size = 30)
plt.title('q1(t)', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)
ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()

plt.plot(t, q2)
plt.grid()
plt.xlabel('t,sec', size = 30)
plt.ylabel('q2,rad', size = 30)
plt.title('q2(t)', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)
ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()

plt.plot(t, q3)
plt.grid()
plt.xlabel('t,sec', size = 30)
plt.ylabel('q3,rad', size = 30)
plt.title('q3(t)', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)
ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()

plt.plot(t, q4)
plt.grid()
plt.xlabel('t,sec', size = 30)
plt.ylabel('q4,rad', size = 30)
plt.title('q4(t)', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)
ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()

plt.plot(t, q5)
plt.grid()
plt.xlabel('t,sec', size = 30)
plt.ylabel('q5,rad', size = 30)
plt.title('q5(t)', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)
ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()

plt.plot(t, q6)
plt.grid()
plt.xlabel('t,sec', size = 30)
plt.ylabel('q6,rad', size = 30)
plt.title('q6(t)', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)
ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()

plt.plot(t, q7)
plt.grid()
plt.xlabel('t,sec', size = 30)
plt.ylabel('q7,rad', size = 30)
plt.title('q7(t)', size = 30)
plt.xticks(size = 20)
plt.yticks(size = 20)
ax = plt.gca()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()

Begin(Theta,q1,q2,q3,q4,q5,q6,q7)


