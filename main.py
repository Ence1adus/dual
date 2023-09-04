import matplotlib.pyplot as plt

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
# Theta = np.array([ 0.744386, 2.14643, 1.58259, 1.88798, -0.826117, 0.588771, 0.328527])
Theta = np.array([0.423316, 0.66385, 2.30886, 1.49975, 0.392119, 0.729956, 0.781015])
# Theta = np.array([ 3, 1, 4, 5, 6, 7, 2])
#Введем дополниительные вычислительные параметры для метода Ньютона , а именно матрицу ориентации Gd и
# построим прогнозную траекторию движения,по координатам x,y,z
Gd = np.array([[0, 0, 1],[-1, 0, 0],[0, -1, 0]])
#Для начала построим желаемую траекторию движения

coordinatesy = [-0.2,0.0,0.3,-0.2]
coordinatesz= [0.7,0.8,0.7,0.7]
# coordinatesy = [-0.3,0.0,0.2,0.23]
# coordinatesz= [0.2,0.45,0.6,0.43]
xd = 0.3
yd = coordinatesy[0]
zd = coordinatesz[0]
# plt.plot(coordinatesy,coordinatesz)
# plt.grid()
#
# plt.xlabel('Y,m', size = 30)
# plt.ylabel('Z,m', size = 30)
# plt.title('Траектория движения в плоскости Y-Z', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
#
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()

n=len(coordinatesy)
a=500 #Разбиение на 200 точек
V = 0.1 # Модуль скорости
tt = [0]*n
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
a1 = int((tt[n-3]/tt[n-1])*a)
a2 = int(((tt[n-2]-tt[n-3])/tt[n-1])*a)
a3 = int(((tt[n-1]-tt[n-2])/tt[n-1])*a)
tcon = a1+a2+a3
t1 = np.linspace(0, tt[n-3] , a1)
t2 = np.linspace(tt[n-3] , tt[n-2] , a2)
t3 = np.linspace(tt[n-2] , tt[n-1] , a3)
t = np.concatenate((t1, t2,t3), axis=None)
dt = tt[n-1] /3*a
Vy1 =[]
Vy2 =[]
Vy3 =[]
Vyp1 =[]
Vyp2 =[]
Vyp3 =[]
Vz1 =[]
Vz2 =[]
Vz3 =[]
Vzp1 =[]
Vzp2 =[]
Vzp3 =[]
for i in range(0, a1):
    Vy1.append(Vy[0])
    Vyp1.append(0)
    Vz1.append(Vz[0])
    Vzp1.append(0)
for i in range(0, a2):
    Vy2.append(Vy[1])
    Vyp2.append((Vy[0]*tt[1]-Vy[1]*tt[1]))
    Vz2.append(Vz[1])
    Vzp2.append((Vz[0]*tt[1]-Vz[1]*tt[1]))
for i in range(0, a3):
    Vy3.append(Vy[2])
    Vyp3.append((Vy[0] * tt[1] - Vy[1] * tt[1]) + (Vy[1] * tt[2] - Vy[2] * tt[2]))
    Vz3.append(Vz[2])
    Vzp3.append((Vz[0]*tt[1]-Vz[1]*tt[1])+(Vz[1]*tt[2]-Vz[2]*tt[2]))

Vymas = np.concatenate((Vy1, Vy2,Vy3), axis=None)
Vypmas = np.concatenate((Vyp1, Vyp2,Vyp3), axis=None)
Vzmas = np.concatenate((Vz1, Vz2,Vz3), axis=None)
Vzpmas = np.concatenate((Vzp1, Vzp2,Vzp3), axis=None)

# for i in range(2*a,3*a):
#     print(Vymas[i])
#     print(Vypmas[i])
# Vyd = np.piecewise(t, [((t >= tt[0]) & (t < tt[1])), ((t >= tt[1]) & (t < tt[2])),((t >= tt[2]) & (t < tt[3]))], [Vy[0],Vy[1],Vy[2]])
# Vzd = np.piecewise(t, [((t >= tt[0]) & (t < tt[1])), ((t >= tt[1]) & (t < tt[2])),((t >= tt[2]) & (t < tt[3]))], [Vz[0],Vz[1],Vz[2]])
#Сформируем свободные слагаем для функции пишвайс
# Vydpom = np.piecewise(t,[((t >= tt[0]) & (t < tt[1])), ((t >= tt[1]) & (t < tt[2])),((t >= tt[2]) & (t < tt[3]))],
#                       [0,(Vy[0]*tt[1]-Vy[1]*tt[1]),(Vy[0]*tt[1]-Vy[1]*tt[1])+(Vy[1]*tt[2]-Vy[2]*tt[2])])
# Vzdpom = np.piecewise(t,[((t >= tt[0]) & (t < tt[1])), ((t >= tt[1]) & (t < tt[2])),((t >= tt[2]) & (t < tt[3]))],
#                       [0,(Vz[0]*tt[1]-Vz[1]*tt[1]),(Vz[0]*tt[1]-Vz[1]*tt[1])+(Vz[1]*tt[2]-Vz[2]*tt[2])])
#Создадим массив точек прогнозного движения
ry = yd
rx = xd
rz = zd
Ryd = [ry]
Rzd = [rz]
Rxd = [rx]
for i in range(1,tcon):
    ryy = ry + Vypmas[i] + Vymas[i] * t[i]
    rzz = rz + Vzpmas[i] + Vzmas[i] * t[i]
    rxx = xd
    Ryd.append(ryy)
    Rzd.append(rzz)
    Rxd.append(rxx)
# Построение графиков
# figure, axis = plt.subplots(2, 2)

# plt.xlabel('Y,m', size = 30)
# plt.ylabel('Z,m', size = 30)
# plt.title('Траектория движения в плоскости Y-Z', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)

# plt.plot(t, Rxd)
# plt.grid()
# plt.xlabel('t,sec', size = 30)
# plt.ylabel('X,m', size = 30)
# plt.title('X(t)', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()
#
# plt.plot(t, Ryd)
# plt.grid()
# plt.xlabel('t,sec', size = 30)
# plt.ylabel('Y,m', size = 30)
# plt.title('Y(t)', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()
#
# plt.plot(t, Rzd)
# plt.grid()
# plt.xlabel('t,sec', size = 30)
# plt.ylabel('Z,m', size = 30)
# plt.title('Z(t)', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()

# axis[0, 0].plot(t, Rxd)
# axis[0, 0].set_title("X(t)",size = 30)
# axis[0, 0].set_ylabel('X,m', size = 30)
# axis[0, 0].set_xlabel('t,sec', size = 30)
# for tick in axis.xaxis.get_major_ticks():
#     tick.label.set_fontsize(40)
#
# # axis[0, 0].set_xticks(size = 20)
# # axis[0, 0].set_yticks(size = 20)
#
# axis[0, 1].plot(t, Ryd)
# axis[0, 1].set_title("Y(t)")
# axis[0, 1].set_ylabel('Y,m')
# axis[0, 1].set_xlabel('t,sec')
#
# axis[1, 0].plot(t, Rzd)
# #axis[1, 0].set_title("Z(t)")
# axis[1, 0].set_ylabel('Z,m')
# axis[1, 0].set_xlabel('t,sec')
# plt.show()
#
# plt.scatter(t,Ryd)
# plt.grid()
# ax = plt.gca()
# # plot X - axis
# ax.axhline(y=0, color='k')
# # plot Y - axis
# ax.axvline(x=0, color='k')
# #display the graph
# plt.show()
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
for i in range(1,tcon):
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
    dG = Gd - Gamma
    dG = np.dot(dG, np.transpose(Gd))
    Betta = [dG[2, 1], dG[0, 2], dG[1, 0]]
    k=0
    while(((LA.norm(dr) > 0.001) or (LA.norm(Betta) > 0.001)) and k<10 ):
        k=k+1
        J0 = ConMat
        #Перевести точки dr к точек 0
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


# plt.scatter(Ynew, Znew)
# plt.grid()
# ax = plt.gca()
# # plot X - axis
# ax.axhline(y=0, color='r')
# # plot Y - axis
# ax.axvline(x=0, color='r')
# #display the graph
# plt.show()

for i in range(0,tcon):
    q1.append(q[i][0])
    q2.append(q[i][1])
    q3.append(q[i][2])
    q4.append(q[i][3])
    q5.append(q[i][4])
    q6.append(q[i][5])
    q7.append(q[i][6])

# plt.plot(t, q7)
# plt.grid()
# ax = plt.gca()
# # plot X - axis
# ax.axhline(y=0, color='r')
# # plot Y - axis
# ax.axvline(x=0, color='r')
# ax.set_xlabel('t,sec',fontsize=16)
# ax.set_ylabel('q7,rad',fontsize=16)
# #display the graph
# plt.show()
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
Begin(Theta,q1,q2,q3,q4,q5,q6,q7)
# Animation(q1,q2,q3,q4,q5,q6,q7)
#построение отдельных графиков
# plt.plot(Ynew, Znew)
# plt.grid()
# plt.xlabel('Y,m', size = 30)
# plt.ylabel('Z,m', size = 30)
# plt.title('Отработанная траектория Y-Z', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()
#
# plt.plot(t, q1)
# plt.grid()
# plt.xlabel('t,sec', size = 30)
# plt.ylabel('q1,rad', size = 30)
# plt.title('q1(t)', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()
#
# plt.plot(t, q2)
# plt.grid()
# plt.xlabel('t,sec', size = 30)
# plt.ylabel('q2,rad', size = 30)
# plt.title('q2(t)', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()
#
# plt.plot(t, q3)
# plt.grid()
# plt.xlabel('t,sec', size = 30)
# plt.ylabel('q3,rad', size = 30)
# plt.title('q3(t)', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()
#
# plt.plot(t, q4)
# plt.grid()
# plt.xlabel('t,sec', size = 30)
# plt.ylabel('q4,rad', size = 30)
# plt.title('q4(t)', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()
#
# plt.plot(t, q5)
# plt.grid()
# plt.xlabel('t,sec', size = 30)
# plt.ylabel('q5,rad', size = 30)
# plt.title('q5(t)', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()
#
# plt.plot(t, q6)
# plt.grid()
# plt.xlabel('t,sec', size = 30)
# plt.ylabel('q6,rad', size = 30)
# plt.title('q6(t)', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()
#
# plt.plot(t, q7)
# plt.grid()
# plt.xlabel('t,sec', size = 30)
# plt.ylabel('q7,rad', size = 30)
# plt.title('q7(t)', size = 30)
# plt.xticks(size = 20)
# plt.yticks(size = 20)
# ax = plt.gca()
# ax.axhline(y=0, color='k')
# ax.axvline(x=0, color='k')
# plt.show()





# D = [0.315, 0, 0.310, 0, 0.300, 0, 0.205]
# a = [0.05, 0, 0, 0, 0, 0, 0]
# q = []
# q.append(sympy.var("q1"))
# q.append(sympy.var("q2"))
# q.append(sympy.var("q3"))
# q.append(sympy.var("q4"))
# q.append(sympy.var("q5"))
# q.append(sympy.var("q6"))
# q.append(sympy.var("q7"))
# print(q)
#
# Theta =[ 0.744386, 2.14643, 1.58259, 1.88798, -0.826117, 0.588771, 0.328527]
# alpha = [(numpy.pi)/2 ,-(numpy.pi)/2,(numpy.pi)/2,-(numpy.pi)/2,
#          (numpy.pi)/2,-(numpy.pi)/2,0]
# d1 = 0.315
# theta1 = 10
# a1 = 0.05
# alpha1 = 4
# G = []
# GG = []
# Gabs = DualMatrix.eye(3)
#Составление абсолютного верзора
# for i in range(7):
#     Rzz = ModuleOfRotationMatrices.RMz(D[i],Theta[i])
#     Rxx = ModuleOfRotationMatrices.RMx(a[i],alpha[i])
#     Gper = Rzz*Rxx
#     Gabs = Gabs*Gper
#     G.append(Gper)
#     GG.append(Gabs)
# print(Gabs)
# Matrix = DualMatrix.DTM(Gabs,True)
# Part2 = Matrix[0]
# Part1 = Matrix[1]
# VMatrix = np.dot(Part1,np.transpose(Part2))
# print("|||||||||")
# print(VMatrix)
# print("|||||||||")
# print("x = " + str(VMatrix[2][1]))
# print("y = " + str(VMatrix[0][2]))
# print("z = " + str(VMatrix[1][0]))
# #Составление якобиана
# print("Составление дуального якобиана")
# IdenRotVector = DualMatrix([[DualNumber(0,0)],[DualNumber(0,0)],[DualNumber(1,0)]])
# Jac = []
# Jac.append(IdenRotVector)
# for i in range(len(GG)-1):
#     StolbecJac = GG[i] * IdenRotVector
#     Jac.append(StolbecJac)
# JacCon = DualMatrix.CrMaFrCol(Jac)
# print(JacCon)
# Mat = DualMatrix.DTM(JacCon,True)
# ConMat = np.concatenate((np.array(Mat[0]),np.array(Mat[1])),axis=0)

#Построим траекторию движения
# coordinatesx = [-0.3,0.5]
# coordinatesy= [0.2,0.6]
# coordinatesx = [-0.3,0.2,0.5,-0.3]
# coordinatesy= [0.2,0.6,0.2,0.2]
# xd = 0.3
# yd = coordinatesx[0]
# zd = coordinatesy[0]
# # plt.plot(coordinatesx,coordinatesy)
# # plt.grid()
# # ax = plt.gca()
# # # plot X - axis
# # ax.axhline(y=0, color='k')
# # # plot Y - axis
# # ax.axvline(x=0, color='k')
# # #display the graph
# # plt.show()
# n=len(coordinatesx)
# a=200
# V = 0.1
# tt = [0]*n
# Vxd = 0
# Vy = [0]*(n-1)
# Vz = [0]*(n-1)
# for i in range(0,n-1):
#     chisnormy = coordinatesx[i+1]-coordinatesx[i]
#     chisnormz = coordinatesy[i + 1] - coordinatesy[i]
#     chismass=[]
#     chismass.append((chisnormy,chisnormz))
#     T =LA.norm(chismass) / V
#     tt[i+1]=tt[i]+T
#     Vy[i] = chisnormy / T
#     Vz[i] = chisnormz / T
# t = np.linspace(0, tt[n-1], a)
# dt = tt[n-1] /a
# Vxd = np.linspace(0,0,a)
# Vyd = np.piecewise(t,[((t >= tt[0]) & (t < tt[1])), ((t >= tt[1]) & (t < tt[2])),((t >= tt[2]) & (t < tt[3]))], [Vy[0],Vy[1],Vy[2]])
# Vzd = np.piecewise(t,[((t >= tt[0]) & (t < tt[1])), ((t >= tt[1]) & (t < tt[2])),((t >= tt[2]) & (t < tt[3]))], [Vz[0],Vz[1],Vz[2]])
# #Сформируем свободные слагаем для функции пишвайс
# Yyd = np.piecewise(t,[((t >= tt[0]) & (t < tt[1])), ((t >= tt[1]) & (t < tt[2])),((t >= tt[2]) & (t < tt[3]))], [0,(Vy[0]*tt[1]-Vy[1]*tt[1]),(Vy[0]*tt[1]-Vy[1]*tt[1])+(Vy[1]*tt[2]-Vy[2]*tt[2])])
# Zyd = np.piecewise(t,[((t >= tt[0]) & (t < tt[1])), ((t >= tt[1]) & (t < tt[2])),((t >= tt[2]) & (t < tt[3]))], [0,(Vz[0]*tt[1]-Vz[1]*tt[1]),(Vz[0]*tt[1]-Vz[1]*tt[1])+(Vz[1]*tt[2]-Vz[2]*tt[2])])
# # plt.scatter(t,Zyd)
# # plt.grid()
# # ax = plt.gca()
# # # plot X - axis
# # ax.axhline(y=0, color='r')
# # # plot Y - axis
# # ax.axvline(x=0, color='r')
# # #display the graph
# # plt.show()
# # Vyd = np.piecewise(t,[((t >= tt[0]) & (t < tt[1])), ((t >= tt[1]) & (t < tt[2]))], [Vy[0],Vy[1]])
# # Vzd = np.piecewise(t,[((t >= tt[0]) & (t < tt[1])), ((t >= tt[1]) & (t < tt[2]))], [Vz[0],Vz[1]])
#
# # Vyd = np.piecewise(t,[((t >= tt[0]) & (t < tt[1]))], [Vy[0]])
# # Vzd = np.piecewise(t,[((t >= tt[0]) & (t < tt[1]))], [Vz[0]])
# #Сформируем массив координа
# ry = yd
# rx = xd
# rz = zd
# Ryd = [ry]
# Rzd = [rz]
# Rxd = [rx]
# for i in range(1,a):
#     ryy = ry + Yyd[i] + Vyd[i]*t[i]
#     rzz = rz + Zyd[i] + Vzd[i] * t[i]
#     rx = xd
#     Ryd.append(ryy)
#     Rzd.append(rzz)
#     Rxd.append(rx)
#
#

# Gd = [[0, 0, 1],
#       [-1, 0, 0],
#       [0, -1, 0]]
# dr = np.array([0,0,0])
# Betta =np.array([0,0,0])
# deltau =np.array([0,0,0,0,0,0])
# q0 =np.array([0,0,0,0,0,0,0])
# xx=[]
# yy=[]
# zz=[]
# x=0
# y=0
# z=0
# xad=0
# yad=0
# zad=0
# q1=[]
# q2=[]
# q3=[]
# q4=[]
# q5=[]
# q6=[]
# q7=[]
# pq =0
# q=[]
# Theta = np.array([0.744386, 2.14643, 1.58259, 1.88798, -0.826117, 0.588771, 0.328527])
# q0 = Theta
# q.append(q0)
# for i in range(0,len(Vyd)):
#     Gamma = []
#     D = [0.315, 0, 0.310, 0, 0.300, 0, 0.205]
#     a = [0.05, 0, 0, 0, 0, 0, 0]
#     alpha = [(np.pi) / 2, -(np.pi) / 2, (np.pi) / 2, -(np.pi) / 2,
#              (np.pi) / 2, -(np.pi) / 2, 0]
#     G = []
#     GG = []
#     Gabs = DualMatrix.eye(3)
#     # Составление абсолютного верзора
#     for k in range(7):
#         Rzz = ModuleOfRotationMatrices.RMz(D[k], Theta[k])
#         Rxx = ModuleOfRotationMatrices.RMx(a[k], alpha[k])
#         Gper = Rzz * Rxx
#         Gabs = Gabs * Gper
#         G.append(Gper)
#         GG.append(Gabs)
#      #Формирование матрицы Гамма
#     Gamma = (DualMatrix.DTM(Gabs,False))[0]
#
#     IdenRotVector = DualMatrix([[DualNumber(0,0)],[DualNumber(0,0)],[DualNumber(1,0)]])
#     Jac = []
#     Jac.append(IdenRotVector)
#     for k in range(len(GG)-1):
#         StolbecJac = GG[k] * IdenRotVector
#         Jac.append(StolbecJac)
#     JacCon = DualMatrix.CrMaFrCol(Jac)
#     # print(JacCon)
#     Mat = DualMatrix.DTM(JacCon,False)
#     #Окончательное формирование якобиана
#     ConMat = np.concatenate((np.array(Mat[0]),np.array(Mat[1])),axis=0)
#     J0 = ConMat
#     #u = [0, 0, 0, Vxd[i], Vyd[i], Vzd[i]]
#     #pq = np.array(np.dot(LA.pinv(ConMat),u))
#     #Theta = Theta +pq*dt
#     #q.append(Theta)
#     #Составление массива координат выходного звена
#     Matrix = DualMatrix.DTM(Gabs,False)
#     Part2 = Matrix[0]
#     Part1 = Matrix[1]
#     VMatrix = np.dot(Part1,np.transpose(Part2))
#     #Формирование вектора RAo
#     x = VMatrix[2][1]
#     y = VMatrix[0][2]
#     z = VMatrix[1][0]
#     # Формирование вектора RAd/t.->tm
#     xad = Rxd[i]
#     yad = Ryd[i]
#     zad = Rzd[i]
#     #Формирование вектора dr
#     dr = [xad-x, yad-y, zad-z]
#     dGamma = np.array(Gd) - Gamma
#     dGamma = np.dot(dGamma,np.transpose(dGamma))
#     Betta = [dGamma[2,1],dGamma[0,2],dGamma[1,0]]
#     j=0
#     while((LA.norm(Betta)>0.001 or LA.norm(dr)>0.001) and j<5):
#
#         deltau = np.concatenate((Betta,dr), axis=None)
#         vichet = np.dot(LA.pinv(J0),-deltau)
#         q0 = Theta - vichet
#         #q0 = Theta - LA.solve(J0,-deltau)
#         Theta = q0
#         G=[]
#         GG=[]
#         #Повторный вычислительный блок
#         # Составление абсолютного верзора
#         for k in range(7):
#             Rzz = ModuleOfRotationMatrices.RMz(D[k], Theta[k])
#             Rxx = ModuleOfRotationMatrices.RMx(a[k], alpha[k])
#             Gper = Rzz * Rxx
#             Gabs = Gabs * Gper
#             G.append(Gper)
#             GG.append(Gabs)
#         # Формирование матрицы Гамма
#         Gamma = (DualMatrix.DTM(Gabs, False))[0]
#
#         IdenRotVector = DualMatrix([[DualNumber(0, 0)], [DualNumber(0, 0)], [DualNumber(1, 0)]])
#         Jac = []
#         ConMat = []
#         Mat =[]
#         Jac.append(IdenRotVector)
#         for k in range(len(GG) - 1):
#             StolbecJac = GG[k] * IdenRotVector
#             Jac.append(StolbecJac)
#         JacCon = DualMatrix.CrMaFrCol(Jac)
#         # print(JacCon)
#         Mat = DualMatrix.DTM(JacCon, False)
#         # Окончательное формирование якобиана
#         ConMat = np.concatenate((np.array(Mat[0]), np.array(Mat[1])), axis=0)
#
#         # u = [0, 0, 0, Vxd[i], Vyd[i], Vzd[i]]
#         # pq = np.array(np.dot(LA.pinv(ConMat),u))
#         # Theta = Theta +pq*dt
#         # q.append(Theta)
#         # Составление массива координат выходного звена
#         Matrix = DualMatrix.DTM(Gabs, False)
#         Part2 = Matrix[0]
#         Part1 = Matrix[1]
#         VMatrix = np.dot(Part1, np.transpose(Part2))
#         #Конец
#         # Формирование вектора RAo
#         x = VMatrix[2][1]
#         y = VMatrix[0][2]
#         z = VMatrix[1][0]
#         # Формирование вектора RAd/t.->tm
#         # Формирование вектора dr
#         dr = [xad - x, yad - y, zad - z]
#         dGamma = np.array(Gd) - Gamma
#         dGamma = np.dot(dGamma, np.transpose(dGamma))
#         Betta = [dGamma[2, 1], dGamma[0, 2], dGamma[1, 0]]
#         j=j+1
#     q.append(q0)
#     xx.append(x)
#     yy.append(y)
#     zz.append(z)
#
#
# for i in range(0,len(Vyd)):
#     q1.append(q[i][0])
#     q2.append(q[i][1])
#     q3.append(q[i][2])
#     q4.append(q[i][3])
#     q5.append(q[i][4])
#     q6.append(q[i][5])
#     q7.append(q[i][6])
#     print("||||||||||||")
#     print((Rxd[i],Ryd[i],Rzd[i]))
#     print((xx[i],yy[i],zz[i]))
#     print("||||||||||||")
#
# # print(q2)
# # print(len(q2))
# #
# plt.scatter(yy,zz)
# plt.grid()
# ax = plt.gca()
# # plot X - axis
# ax.axhline(y=0, color='r')
# # plot Y - axis
# ax.axvline(x=0, color='r')
# #display the graph
# plt.show()



 # create function
# def f(y, t):
#     y1, y2 = y
#     return[y2, -y2 - y1]
#
# t = np.linspace( 0, 10, 41) # vector of time
# y0 = [0, 1] # start value
# w = odeint(f, y0, t) # solve eq.
# y1 = w[:,0]
# y2 = w[:,1]
# fig = plt.figure(facecolor='white')
# plt.plot(t, y1, '-o', t, y2, '-o', linewidth=2)
# plt.ylabel("z")
# plt.xlabel("t")
# plt.grid(True)
# plt.show() # di


# Rzz = ModuleOfRotationMatrices.RMz(D[1],Theta[1])
# Rxx = ModuleOfRotationMatrices.RMx(a[1],alpha[1])
# G12 = Rzz*Rxx
# Rzz = ModuleOfRotationMatrices.RMz(D[2],Theta[2])
# Rxx = ModuleOfRotationMatrices.RMx(a[2],alpha[2])
# G23 = Rzz*Rxx
# Rzz = ModuleOfRotationMatrices.RMz(D[3],Theta[3])
# Rxx = ModuleOfRotationMatrices.RMx(a[3],alpha[3])
# G34 = Rzz*Rxx
# Rzz = ModuleOfRotationMatrices.RMz(D[4],Theta[4])
# Rxx = ModuleOfRotationMatrices.RMx(a[4],alpha[4])
# G45 = Rzz*Rxx

#print(Rzz)
#print(Rxx)
#print(G01)
#DualMatrix.DTM(G01)
# dr = dr + np.array([0,0,1])
# Перевести точки dr к точек 0