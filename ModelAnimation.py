import numpy as np
import math as m
import matplotlib.pyplot as plt
import vpython as vp
from pyquaternion import Quaternion
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.animation as animation
from matplotlib.patches import Circle
from DualNumberModule import DualNumber as dn
from DualNumberModule import DualMatrix as dm
import ModuleOfRotationMatrices as rm
#Создание модели основной
def Calibrate(Theta,zv1_1,zv1_2,zv2_1,zv2_2,zv3_1,zv4_1,zv4_2,zv5_1,zv6_1,zv6_2,zv7_1,zv7_2):
    D = np.array([0.315, 0, 0.310, 0, 0.300, 0, 0.205])
    aa = np.array([0.05, 0, 0, 0, 0, 0, 0])
    alpha = np.array([(np.pi) / 2, -(np.pi) / 2, (np.pi) / 2, -(np.pi) / 2,
                      (np.pi) / 2, -(np.pi) / 2, 0])
    GG = []
    Gabs = dm.eye(3)
    for j in range(7):
        Rzm = rm.RMz(D[j], Theta[j])
        Rxm = rm.RMx(aa[j], alpha[j])
        Gper = Rzm * Rxm
        Gabs = Gabs * Gper
        GG.append(Gabs)
    # Вращение первого шарнира
    zv1_1.rotate(angle=Theta[0], axis=vp.vector(0, 0, 1))
    zv1_2.rotate(angle=Theta[0], axis=vp.vector(0, 0, 1))
    zv2_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
    zv2_2.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
    zv3_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
    zv4_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
    zv4_2.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
    zv5_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
    zv6_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
    zv6_2.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
    zv7_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
    zv7_2.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
    # Вращение второго шарнира
    vec = np.array([0, 0, 1])
    Matrix = dm.DTM(GG[0], False)
    VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
    xorigin = VMatrix[2][1]
    yorigin = VMatrix[0][2]
    zorigin = VMatrix[1][0]
    MatPov = np.dot(np.array(Matrix[0]), vec)
    zv2_1.rotate(angle=Theta[1], axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv2_2.rotate(angle=Theta[1], axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv3_1.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                       axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv4_1.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                       axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv4_2.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                       axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv5_1.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv6_1.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv6_2.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv7_1.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv7_2.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    # Вращение звена 3
    vec = np.array([0, 0, 1])
    Matrix = dm.DTM(GG[1], False)
    VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
    xorigin = VMatrix[2][1]
    yorigin = VMatrix[0][2]
    zorigin = VMatrix[1][0]
    MatPov = np.dot(np.array(Matrix[0]), vec)
    zv3_1.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                       axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv4_1.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                       axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv4_2.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                       axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv5_1.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv6_1.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv6_2.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv7_1.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv7_2.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    # Вразение 4го звена
    vec = np.array([0, 0, 1])
    Matrix = dm.DTM(GG[2], False)
    VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
    xorigin = VMatrix[2][1]
    yorigin = VMatrix[0][2]
    zorigin = VMatrix[1][0]
    MatPov = np.dot(np.array(Matrix[0]), vec)
    zv4_1.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                       axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv4_2.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                       axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv5_1.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv6_1.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv6_2.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv7_1.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv7_2.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    # Вразение 5го звена
    vec = np.array([0, 0, 1])
    Matrix = dm.DTM(GG[3], False)
    VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
    xorigin = VMatrix[2][1]
    yorigin = VMatrix[0][2]
    zorigin = VMatrix[1][0]
    MatPov = np.dot(np.array(Matrix[0]), vec)
    zv5_1.rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv6_1.rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv6_2.rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv7_1.rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv7_2.rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    # Вразение 6го звена
    vec = np.array([0, 0, 1])
    Matrix = dm.DTM(GG[4], False)
    VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
    xorigin = VMatrix[2][1]
    yorigin = VMatrix[0][2]
    zorigin = VMatrix[1][0]
    MatPov = np.dot(np.array(Matrix[0]), vec)
    zv6_1.rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv6_2.rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv7_1.rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv7_2.rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    # Вразение 7го звена
    vec = np.array([0, 0, 1])
    Matrix = dm.DTM(GG[5], False)
    VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
    xorigin = VMatrix[2][1]
    yorigin = VMatrix[0][2]
    zorigin = VMatrix[1][0]
    MatPov = np.dot(np.array(Matrix[0]), vec)
    zv7_1.rotate(angle=Theta[6], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
    zv7_2.rotate(angle=Theta[6], origin=vp.vector(xorigin, yorigin, zorigin),
                        axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))

def Animation(Theta1,Theta2,Theta3,Theta4,Theta5,Theta6,Theta7,zv1_1, zv1_2, zv2_1, zv2_2, zv3_1, zv4_1, zv4_2, zv5_1, zv6_1, zv6_2, zv7_1, zv7_2):
    koef = 1000
    zv1_1.visible = False
    zv1_2.visible = False
    zv2_1.visible = False
    zv2_2.visible = False
    zv3_1.visible = False
    zv4_1.visible = False
    zv4_2.visible = False
    zv5_1.visible = False
    zv6_1.visible = False
    zv6_2.visible = False
    zv7_1.visible = False
    zv7_2.visible = False

    # Приведение в начальное положение моделирования
    # Theta = np.array([0.744386, 2.14643, 1.58259, 1.88798, -0.826117, 0.588771, 0.328527])
    # Theta = np.array([ -0.744386, -2.14643, -1.58259, -1.88798, 0.826117, -0.588771, -0.328527])
    D = np.array([0.315, 0, 0.310, 0, 0.300, 0, 0.205])
    aa = np.array([0.05, 0, 0, 0, 0, 0, 0])
    alpha = np.array([(np.pi) / 2, -(np.pi) / 2, (np.pi) / 2, -(np.pi) / 2,
                      (np.pi) / 2, -(np.pi) / 2, 0])
    xd = 0.3
    yd = -0.3
    zd = 0.2
    # vp.curve(pos=[vp.vector(0, 0, 0), vp.vector(xd, yd, zd)], color=vp.color.cyan)
    Mas = []
    for i in range(len(Theta1)):
        vp.rate(10)
        koef = 1000
        a = []
        osnov1 = vp.cylinder(pos=vp.vector(0, 0, 0),
                             axis=vp.vector(0, 0, 1),
                             length=85 / koef,
                             radius=15 / koef,
                             color=vp.color.cyan)
        a.append(osnov1)
        osnov1.visible = False
        osnov2 = vp.cylinder(pos=vp.vector(0, 0, 85 / koef),
                             axis=vp.vector(0, 0, 1),
                             length=15 / koef,
                             radius=30 / koef,
                             color=vp.color.cyan)
        a.append(osnov2)
        osnov2.visible = False
        osnov3 = vp.cylinder(pos=vp.vector(0, 0, 100 / koef),
                             axis=vp.vector(0, 0, 1),
                             length=50 / koef,
                             radius=15 / koef,
                             color=vp.color.cyan)
        a.append(osnov3)
        osnov3.visible = False
        zv1_1 = vp.cylinder(pos=vp.vector(0, 0, 150 / koef),
                            axis=vp.vector(0, 0, 1),
                            length=(165 + 15) / koef,
                            radius=15 / koef,
                            color=vp.color.red)
        a.append(zv1_1)
        zv1_1.visible = False
        zv1_2 = vp.cylinder(pos=vp.vector(0, 0, 315 / koef),
                            axis=vp.vector(1, 0, 0),
                            length=(50) / koef,
                            radius=15 / koef,
                            color=vp.color.red)
        a.append(zv1_2)
        zv1_2.visible = False
        zv2_1 = vp.cylinder(pos=vp.vector(50 / koef, -17.5 / koef, 315 / koef),
                            axis=vp.vector(0, 1, 0),
                            length=35 / koef,  # Добавить +15
                            radius=15 / koef,
                            color=vp.color.yellow)
        a.append(zv2_1)
        zv2_1.visible = False
        zv2_2 = vp.cylinder(pos=vp.vector(50 / koef, 0, 315 / koef),
                            axis=vp.vector(0, 0, 1),
                            length=100 / koef,
                            radius=15 / koef,
                            color=vp.color.yellow)
        a.append(zv2_2)
        zv2_2.visible = False
        zv3_1 = vp.cylinder(pos=vp.vector(50 / koef, 0, 415 / koef),
                            axis=vp.vector(0, 0, 1),
                            length=(195 + 15) / koef,
                            radius=15 / koef,
                            color=vp.color.magenta)
        a.append(zv3_1)
        zv3_1.visible = False
        zv4_1 = vp.cylinder(pos=vp.vector(50 / koef, -17.5 / koef, 625 / koef),
                            axis=vp.vector(0, 1, 0),
                            length=35 / koef,
                            radius=15 / koef,
                            color=vp.color.magenta)
        a.append(zv4_1)
        zv4_1.visible = False
        zv4_2 = vp.cylinder(pos=vp.vector(50 / koef, 0, 625 / koef),
                            axis=vp.vector(0, 0, 1),
                            length=100 / koef,
                            radius=15 / koef,
                            color=vp.color.green)
        a.append(zv4_2)
        zv4_2.visible = False
        zv5_1 = vp.cylinder(pos=vp.vector(50 / koef, 0, 725 / koef),
                            axis=vp.vector(0, 0, 1),
                            length=200 / koef,
                            radius=15 / koef,
                            color=vp.color.orange)
        a.append(zv5_1)
        zv5_1.visible = False
        zv6_1 = vp.cylinder(pos=vp.vector(50 / koef, -17.5 / koef, 925 / koef),
                            axis=vp.vector(0, 1, 0),
                            length=35 / koef,
                            radius=15 / koef,
                            color=vp.color.orange)
        a.append(zv6_1)
        zv6_1.visible = False
        zv6_2 = vp.cylinder(pos=vp.vector(50 / koef, 0, 925 / koef),
                            axis=vp.vector(0, 0, 1),
                            length=100 / koef,
                            radius=15 / koef,
                            color=vp.color.green)
        a.append(zv6_2)
        zv6_2.visible = False
        zv7_1 = vp.cylinder(pos=vp.vector(50 / koef, 0, 1025 / koef),
                            axis=vp.vector(0, 0, 1),
                            length=90 / koef,
                            radius=15 / koef,
                            color=vp.color.blue)
        a.append(zv7_1)
        zv7_1.visible = False
        zv7_2 = vp.box(pos=vp.vector(50 / koef, 0, (1025 + 90) / koef),
                       length=(30 / m.sqrt(2)) / koef, height=(30 / m.sqrt(2)) / koef, width=30 / koef,
                       color=vp.color.white)
        a.append(zv7_2)
        zv7_2.visible = False
        Mas.append(a)
        Theta = np.array([Theta1[i],Theta2[i], Theta3[i], Theta4[i], Theta5[i], Theta6[i], Theta7[i]])
        D = np.array([0.315, 0, 0.310, 0, 0.300, 0, 0.205])
        aa = np.array([0.05, 0, 0, 0, 0, 0, 0])
        alpha = np.array([(np.pi) / 2, -(np.pi) / 2, (np.pi) / 2, -(np.pi) / 2,
                          (np.pi) / 2, -(np.pi) / 2, 0])
        GG = []
        Gabs = dm.eye(3)
        for j in range(7):
            Rzm = rm.RMz(D[j], Theta[j])
            Rxm = rm.RMx(aa[j], alpha[j])
            Gper = Rzm * Rxm
            Gabs = Gabs * Gper
            GG.append(Gabs)
        # Вращение первого шарнира
        (Mas[i])[3].rotate(angle=Theta[0], axis=vp.vector(0, 0, 1))
        (Mas[i])[4].rotate(angle=Theta[0], axis=vp.vector(0, 0, 1))
        (Mas[i])[5].rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
        (Mas[i])[6].rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
        (Mas[i])[7].rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
        (Mas[i])[8].rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
        (Mas[i])[9].rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
        (Mas[i])[10].rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
        (Mas[i])[11].rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
        (Mas[i])[12].rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
        (Mas[i])[13].rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
        (Mas[i])[14].rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
        # Вращение второго шарнира
        vec = np.array([0, 0, 1])
        Matrix = dm.DTM(GG[0], False)
        VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
        xorigin = VMatrix[2][1]
        yorigin = VMatrix[0][2]
        zorigin = VMatrix[1][0]
        MatPov = np.dot(np.array(Matrix[0]), vec)
        (Mas[i])[5].rotate(angle=Theta[1], axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[6].rotate(angle=Theta[1], axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[7].rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[8].rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[9].rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[10].rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[11].rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[12].rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[13].rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[14].rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        # Вращение звена 3
        vec = np.array([0, 0, 1])
        Matrix = dm.DTM(GG[1], False)
        VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
        xorigin = VMatrix[2][1]
        yorigin = VMatrix[0][2]
        zorigin = VMatrix[1][0]
        MatPov = np.dot(np.array(Matrix[0]), vec)
        (Mas[i])[7].rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[8].rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[9].rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[10].rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[11].rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[12].rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[13].rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[14].rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        # Вразение 4го звена
        vec = np.array([0, 0, 1])
        Matrix = dm.DTM(GG[2], False)
        VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
        xorigin = VMatrix[2][1]
        yorigin = VMatrix[0][2]
        zorigin = VMatrix[1][0]
        MatPov = np.dot(np.array(Matrix[0]), vec)
        (Mas[i])[8].rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[9].rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[10].rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[11].rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[12].rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[13].rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[14].rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        # Вразение 5го звена
        vec = np.array([0, 0, 1])
        Matrix = dm.DTM(GG[3], False)
        VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
        xorigin = VMatrix[2][1]
        yorigin = VMatrix[0][2]
        zorigin = VMatrix[1][0]
        MatPov = np.dot(np.array(Matrix[0]), vec)
        (Mas[i])[10].rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[11].rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[12].rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[13].rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[14].rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        # Вразение 6го звена
        vec = np.array([0, 0, 1])
        Matrix = dm.DTM(GG[4], False)
        VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
        xorigin = VMatrix[2][1]
        yorigin = VMatrix[0][2]
        zorigin = VMatrix[1][0]
        MatPov = np.dot(np.array(Matrix[0]), vec)
        (Mas[i])[11].rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[12].rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[13].rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[14].rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        # Вразение 7го звена
        vec = np.array([0, 0, 1])
        Matrix = dm.DTM(GG[5], False)
        VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
        xorigin = VMatrix[2][1]
        yorigin = VMatrix[0][2]
        zorigin = VMatrix[1][0]
        MatPov = np.dot(np.array(Matrix[0]), vec)
        (Mas[i])[13].rotate(angle=Theta[6], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        (Mas[i])[14].rotate(angle=Theta[6], origin=vp.vector(xorigin, yorigin, zorigin),
                     axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
        Matrix = dm.DTM(Gabs, False)
        VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
        x = VMatrix[2][1]
        y = VMatrix[0][2]
        z = VMatrix[1][0]
        vp.sphere(pos=vp.vector(x, y, z), radius= 1/ koef, color = vp.color.cyan)
        (Mas[i])[0].visible = True
        (Mas[i])[1].visible = True
        (Mas[i])[2].visible = True
        (Mas[i])[3].visible = True
        (Mas[i])[4].visible = True
        (Mas[i])[5].visible = True
        (Mas[i])[6].visible = True
        (Mas[i])[7].visible = True
        (Mas[i])[8].visible = True
        (Mas[i])[9].visible = True
        (Mas[i])[10].visible = True
        (Mas[i])[11].visible = True
        (Mas[i])[12].visible = True
        (Mas[i])[13].visible = True
        (Mas[i])[14].visible = True
        if i>0:
            (Mas[i - 1])[0].visible = False
            (Mas[i - 1])[1].visible = False
            (Mas[i - 1])[2].visible = False
            (Mas[i - 1])[3].visible = False
            (Mas[i - 1])[4].visible = False
            (Mas[i - 1])[5].visible = False
            (Mas[i - 1])[6].visible = False
            (Mas[i - 1])[7].visible = False
            (Mas[i - 1])[8].visible = False
            (Mas[i - 1])[9].visible = False
            (Mas[i - 1])[10].visible = False
            (Mas[i - 1])[11].visible = False
            (Mas[i - 1])[12].visible = False
            (Mas[i - 1])[13].visible = False
            (Mas[i - 1])[14].visible = False

        # zv1_1.visible = False
        # zv1_2.visible = False
        # zv2_1.visible = False
        # zv2_2.visible = False
        # zv3_1.visible = False
        # zv4_1.visible = False
        # zv4_2.visible = False
        # zv5_1.visible = False
        # zv6_1.visible = False
        # zv6_2.visible = False
        # zv7_1.visible = False
        # zv7_2.visible = False

def Begin(Theta,Theta1,Theta2,Theta3,Theta4,Theta5,Theta6,Theta7):
    vp.canvas(background=vp.color.white)
    koef = 1000
    # vp.scene.camera.pos = vp.vector(50/koef,50/koef,50/koef)
    # vp.scene.camera.axis = -vp.vector(50/koef,50/koef,50/koef)
    # vp.scene.center = vp.vector(50/koef , 250/koef, 250/koef)
    osnov1 = vp.cylinder(pos=vp.vector(0, 0, 0),
                         axis=vp.vector(0, 0, 1),
                         length=85 / koef,
                         radius=15 / koef,
                         color=vp.color.cyan)
    osnov2 = vp.cylinder(pos=vp.vector(0, 0, 85 / koef),
                         axis=vp.vector(0, 0, 1),
                         length=15 / koef,
                         radius=30 / koef,
                         color=vp.color.cyan)
    osnov3 = vp.cylinder(pos=vp.vector(0, 0, 100 / koef),
                         axis=vp.vector(0, 0, 1),
                         length=50 / koef,
                         radius=15 / koef,
                         color=vp.color.cyan)

    zv1_1 = vp.cylinder(pos=vp.vector(0, 0, 150 / koef),
                        axis=vp.vector(0, 0, 1),
                        length=(165 + 15) / koef,
                        radius=15 / koef,
                        color=vp.color.red)
    zv1_2 = vp.cylinder(pos=vp.vector(0, 0, 315 / koef),
                        axis=vp.vector(1, 0, 0),
                        length=(50) / koef,
                        radius=15 / koef,
                        color=vp.color.red)

    zv2_1 = vp.cylinder(pos=vp.vector(50 / koef, -17.5 / koef, 315 / koef),
                        axis=vp.vector(0, 1, 0),
                        length=35 / koef,  # Добавить +15
                        radius=15 / koef,
                        color=vp.color.yellow)
    zv2_2 = vp.cylinder(pos=vp.vector(50 / koef, 0, 315 / koef),
                        axis=vp.vector(0, 0, 1),
                        length=100 / koef,
                        radius=15 / koef,
                        color=vp.color.yellow)

    zv3_1 = vp.cylinder(pos=vp.vector(50 / koef, 0, 415 / koef),
                        axis=vp.vector(0, 0, 1),
                        length=(195 + 15) / koef,
                        radius=15 / koef,
                        color=vp.color.magenta)
    #
    zv4_1 = vp.cylinder(pos=vp.vector(50 / koef, -17.5 / koef, 625 / koef),
                        axis=vp.vector(0, 1, 0),
                        length=35 / koef,
                        radius=15 / koef,
                        color=vp.color.magenta)
    zv4_2 = vp.cylinder(pos=vp.vector(50 / koef, 0, 625 / koef),
                        axis=vp.vector(0, 0, 1),
                        length=100 / koef,
                        radius=15 / koef,
                        color=vp.color.green)

    zv5_1 = vp.cylinder(pos=vp.vector(50 / koef, 0, 725 / koef),
                        axis=vp.vector(0, 0, 1),
                        length=200 / koef,
                        radius=15 / koef,
                        color=vp.color.orange)

    zv6_1 = vp.cylinder(pos=vp.vector(50 / koef, -17.5 / koef, 925 / koef),
                        axis=vp.vector(0, 1, 0),
                        length=35 / koef,
                        radius=15 / koef,
                        color=vp.color.orange)
    zv6_2 = vp.cylinder(pos=vp.vector(50 / koef, 0, 925 / koef),
                        axis=vp.vector(0, 0, 1),
                        length=100 / koef,
                        radius=15 / koef,
                        color=vp.color.green)

    zv7_1 = vp.cylinder(pos=vp.vector(50 / koef, 0, 1025 / koef),
                        axis=vp.vector(0, 0, 1),
                        length=90 / koef,
                        radius=15 / koef,
                        color=vp.color.blue)
    zv7_2 = vp.box(pos=vp.vector(50 / koef, 0, (1025 + 90) / koef),
                   length=(30 / m.sqrt(2)) / koef, height=(30 / m.sqrt(2)) / koef, width=30 / koef,
                   color=vp.color.white)
    # Задаим координатные глобальные оси
    position = vp.vector(0, 0, 0)
    Xvector = vp.curve(pos=[position, vp.vector(1500 / koef, 0, 0)], color=vp.color.green)
    Yvector = vp.curve(pos=[position, vp.vector(0, 1500 / koef, 0)], color=vp.color.red)
    Zvector = vp.curve(pos=[position, vp.vector(0, 0, 1500 / koef)], color=vp.color.blue)
    vp.sleep(4)
    # Calibrate(Theta, zv1_1, zv1_2, zv2_1, zv2_2, zv3_1, zv4_1, zv4_2, zv5_1, zv6_1, zv6_2, zv7_1, zv7_2)
    vp.sleep(2)
    Animation(Theta1,Theta2,Theta3,Theta4,Theta5,Theta6,Theta7,zv1_1, zv1_2, zv2_1, zv2_2, zv3_1, zv4_1, zv4_2, zv5_1, zv6_1, zv6_2, zv7_1, zv7_2)

# koef = 1000
# osnov1 = vp.cylinder(pos = vp.vector(0, 0, 0),
#           axis=vp.vector(0,0,1),
#           length = 85/koef,
#           radius = 15/koef,
#           color = vp.color.cyan)
# osnov2 = vp.cylinder(pos = vp.vector(0, 0, 85/koef),
#           axis=vp.vector(0,0,1),
#           length = 15/koef,
#           radius = 30/koef,
#           color = vp.color.cyan)
# osnov3 = vp.cylinder(pos = vp.vector(0, 0, 100/koef),
#           axis=vp.vector(0,0,1),
#           length = 50/koef,
#           radius = 15/koef,
#           color = vp.color.cyan)
#
# zv1_1 = vp.cylinder(pos = vp.vector(0, 0, 150/koef),
#           axis=vp.vector(0,0,1),
#           length = (165+15)/koef,
#           radius = 15/koef,
#           color = vp.color.red)
# zv1_2 = vp.cylinder(pos = vp.vector(0, 0,315/koef),
#           axis=vp.vector(1,0,0),
#           length = (50)/koef,
#           radius = 15/koef,
#           color = vp.color.red)
#
# zv2_1 = vp.cylinder(pos = vp.vector(50/koef, -17.5/koef,315/koef),
#           axis=vp.vector(0,1,0),
#           length = 35/koef,# Добавить +15
#           radius = 15/koef,
#           color = vp.color.yellow)
# zv2_2 = vp.cylinder(pos = vp.vector(50/koef, 0, 315/koef),
#           axis=vp.vector(0,0,1),
#           length = 100/koef,
#           radius = 15/koef,
#           color = vp.color.yellow)
#
# zv3_1 = vp.cylinder(pos = vp.vector(50/koef, 0, 415/koef),
#           axis=vp.vector(0,0,1),
#           length = (195+15)/koef,
#           radius = 15/koef,
#           color = vp.color.magenta)
# #
# zv4_1 = vp.cylinder(pos = vp.vector(50/koef, -17.5/koef, 625/koef),
#           axis=vp.vector(0,1,0),
#           length =35/koef ,
#           radius = 15/koef,
#           color = vp.color.magenta)
# zv4_2 = vp.cylinder(pos = vp.vector(50/koef, 0, 625/koef),
#           axis=vp.vector(0,0,1),
#           length = 100/koef,
#           radius = 15/koef,
#           color = vp.color.green)
#
# zv5_1 = vp.cylinder(pos = vp.vector(50/koef, 0, 725/koef),
#           axis=vp.vector(0,0,1),
#           length = 200/koef,
#           radius = 15/koef,
#           color = vp.color.orange)
#
# zv6_1 = vp.cylinder(pos = vp.vector(50/koef, -17.5/koef, 925/koef),
#           axis=vp.vector(0,1,0),
#           length =35/koef ,
#           radius = 15/koef,
#           color = vp.color.orange)
# zv6_2 = vp.cylinder(pos = vp.vector(50/koef, 0, 925/koef),
#           axis=vp.vector(0,0,1),
#           length = 100/koef,
#           radius = 15/koef,
#           color = vp.color.green)
#
# zv7_1 = vp.cylinder(pos = vp.vector(50/koef, 0, 1025/koef),
#           axis=vp.vector(0,0,1),
#           length = 90/koef,
#           radius = 15/koef,
#           color = vp.color.blue)
# zv7_2 = vp.box(pos = vp.vector(50/koef, 0, (1025+90)/koef),
#           length=(30/m.sqrt(2))/koef, height=(30/m.sqrt(2))/koef, width=30/koef,
#           color = vp.color.white)
# #Зададим относительный поворот
# dt = 0.1
# omega = 0.0001
# #Задаим координатные глобальные оси
# position = vp.vector(0,0,0)
# Xvector = vp.curve (pos=[position,vp.vector(1500/koef,0,0)],color = vp.color.green)
# Yvector = vp.curve (pos=[position,vp.vector(0,1500/koef,0)],color = vp.color.red)
# Zvector = vp.curve (pos=[position,vp.vector(0,0,1500/koef)],color = vp.color.blue)
#
# #Приведение в начальное положение моделирования
# Theta = np.array([ 0.744386, 2.14643, 1.58259, 1.88798, -0.826117, 0.588771, 0.328527])
# # Theta = np.array([ -0.744386, -2.14643, -1.58259, -1.88798, 0.826117, -0.588771, -0.328527])
# D = np.array([0.315, 0, 0.310, 0, 0.300, 0, 0.205])
# aa = np.array([0.05, 0, 0, 0, 0, 0, 0])
# alpha = np.array([( np.pi ) / 2, -(np.pi)/2, (np.pi)/2, -(np.pi)/2,
#          (np.pi)/2, -(np.pi)/2,0])
# xd = 0.3
# yd = -0.3
# zd = 0.2
# vp.curve (pos=[vp.vector(0,0,0),vp.vector(xd, yd,zd)],color = vp.color.cyan)
# GG = []
# Gabs = dm.eye(3)
# for j in range(7):
#     Rzm = rm.RMz(D[j], Theta[j])
#     Rxm = rm.RMx(aa[j], alpha[j])
#     Gper = Rzm * Rxm
#     Gabs = Gabs*Gper
#     GG.append(Gabs)
# #Вращение первого шарнира
# zv1_1.rotate(angle=Theta[0], axis=vp.vector(0,0,1))
# zv1_2.rotate(angle=Theta[0], axis=vp.vector(0,0,1))
# zv2_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
# zv2_2.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
# zv3_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
# zv4_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
# zv4_2.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
# zv5_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
# zv6_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
# zv6_2.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
# zv7_1.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
# zv7_2.rotate(angle=Theta[0], origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
# #Вращение второго шарнира
# vec = np.array([0,0,1])
# Matrix = dm.DTM(GG[0], False)
# VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
# xorigin = VMatrix[2][1]
# yorigin = VMatrix[0][2]
# zorigin = VMatrix[1][0]
# MatPov = np.dot(np.array(Matrix[0]),vec)
# zv2_1.rotate(angle=Theta[1], axis=vp.vector(MatPov[0], MatPov[1],MatPov[2]))
# zv2_2.rotate(angle=Theta[1], axis=vp.vector(MatPov[0], MatPov[1],MatPov[2]))
# zv3_1.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1],MatPov[2]))
# zv4_1.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1],MatPov[2]))
# zv4_2.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1],MatPov[2]))
# zv5_1.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1],MatPov[2]))
# zv6_1.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1],MatPov[2]))
# zv6_2.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1],MatPov[2]))
# zv7_1.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1],MatPov[2]))
# zv7_2.rotate(angle=Theta[1], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1],MatPov[2]))
# #Вращение звена 3
# vec = np.array([0,0,1])
# Matrix = dm.DTM(GG[1], False)
# VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
# xorigin = VMatrix[2][1]
# yorigin = VMatrix[0][2]
# zorigin = VMatrix[1][0]
# MatPov = np.dot(np.array(Matrix[0]),vec)
# zv3_1.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv4_1.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv4_2.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv5_1.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv6_1.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv6_2.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv7_1.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv7_2.rotate(angle=Theta[2], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# #Вразение 4го звена
# vec = np.array([0,0,1])
# Matrix = dm.DTM(GG[2], False)
# VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
# xorigin = VMatrix[2][1]
# yorigin = VMatrix[0][2]
# zorigin = VMatrix[1][0]
# MatPov = np.dot(np.array(Matrix[0]),vec)
# zv4_1.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv4_2.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv5_1.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv6_1.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv6_2.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv7_1.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv7_2.rotate(angle=Theta[3], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# #Вразение 5го звена
# vec = np.array([0,0,1])
# Matrix = dm.DTM(GG[3], False)
# VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
# xorigin = VMatrix[2][1]
# yorigin = VMatrix[0][2]
# zorigin = VMatrix[1][0]
# MatPov = np.dot(np.array(Matrix[0]),vec)
# zv5_1.rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv6_1.rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv6_2.rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv7_1.rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv7_2.rotate(angle=Theta[4], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# #Вразение 6го звена
# vec = np.array([0,0,1])
# Matrix = dm.DTM(GG[4], False)
# VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
# xorigin = VMatrix[2][1]
# yorigin = VMatrix[0][2]
# zorigin = VMatrix[1][0]
# MatPov = np.dot(np.array(Matrix[0]),vec)
# zv6_1.rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv6_2.rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv7_1.rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv7_2.rotate(angle=Theta[5], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# #Вразение 7го звена
# vec = np.array([0,0,1])
# Matrix = dm.DTM(GG[5], False)
# VMatrix = np.dot(np.array(Matrix[1]), np.transpose(Matrix[0]))
# xorigin = VMatrix[2][1]
# yorigin = VMatrix[0][2]
# zorigin = VMatrix[1][0]
# MatPov = np.dot(np.array(Matrix[0]),vec)
# zv7_1.rotate(angle=Theta[6], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# zv7_2.rotate(angle=Theta[6], origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))

# while(True):
#     #Вращение вокург 1го звена
#     zv1_1.rotate(angle=omega * dt, axis=vp.vector(0,0,1))
#     zv1_2.rotate(angle=omega * dt, axis=vp.vector(0,0,1))
#     zv2_1.rotate(angle=omega * dt, origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
#     zv2_2.rotate(angle=omega * dt, origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
#     zv3_1.rotate(angle=omega * dt, origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
#     zv4_1.rotate(angle=omega * dt, origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
#     zv4_2.rotate(angle=omega * dt, origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
#     zv5_1.rotate(angle=omega * dt, origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
#     zv6_1.rotate(angle=omega * dt, origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
#     zv6_2.rotate(angle=omega * dt, origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
#     zv7_1.rotate(angle=omega * dt, origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))
#     zv7_2.rotate(angle=omega * dt, origin=vp.vector(0, 0, 1), axis=vp.vector(0, 0, 1))

# while(True):
#     # Вращение вокург 2го звена
#     zv2_1.rotate(angle=omega * dt, axis=vp.vector(xorigin, yorigin, 0))
#     zv2_2.rotate(angle=omega * dt, axis=vp.vector(xorigin, yorigin, 0))
#     zv3_1.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(xorigin, yorigin, 0))
#     zv4_1.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(xorigin, yorigin, 0))
#     zv4_2.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(xorigin, yorigin, 0))
#     zv5_1.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(xorigin, yorigin, 0))
#     zv6_1.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(xorigin, yorigin, 0))
#     zv6_2.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(xorigin, yorigin, 0))
#     zv7_1.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(xorigin, yorigin, 0))
#     zv7_2.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin), axis=vp.vector(xorigin, yorigin, 0))

# while(True):
#     zv3_1.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin),
#                  axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
#     zv4_1.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin),
#                  axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
#     zv4_2.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin),
#                  axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
#     zv5_1.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin),
#                  axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
#     zv6_1.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin),
#                  axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
#     zv6_2.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin),
#                  axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
#     zv7_1.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin),
#                  axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
#     zv7_2.rotate(angle=omega * dt, origin=vp.vector(xorigin, yorigin, zorigin),
#                  axis=vp.vector(MatPov[0], MatPov[1], MatPov[2]))
# while(True):
#     # Вращение вокург 4го звена
#     zv4_1.rotate(angle=omega * dt, origin=vp.vector(Vorin3[0], Vorin3[1],Vorin3[2]), axis=vp.vector(Veorinnow1[0], Veorinnow1[1], Veorinnow1[2]))
#     zv4_2.rotate(angle=omega * dt, origin=vp.vector(Vorin3[0], Vorin3[1],Vorin3[2]), axis=vp.vector(Veorinnow1[0], Veorinnow1[1], Veorinnow1[2]))
#     zv5_1.rotate(angle=omega * dt, origin=vp.vector(Vorin3[0], Vorin3[1],Vorin3[2]), axis=vp.vector(Veorinnow1[0], Veorinnow1[1], Veorinnow1[2]))
#     zv6_1.rotate(angle=omega * dt, origin=vp.vector(Vorin3[0], Vorin3[1],Vorin3[2]), axis=vp.vector(Veorinnow1[0], Veorinnow1[1], Veorinnow1[2]))
#     zv6_2.rotate(angle=omega * dt, origin=vp.vector(Vorin3[0], Vorin3[1],Vorin3[2]), axis=vp.vector(Veorinnow1[0], Veorinnow1[1], Veorinnow1[2]))
#     zv7_1.rotate(angle=omega * dt, origin=vp.vector(Vorin3[0], Vorin3[1],Vorin3[2]), axis=vp.vector(Veorinnow1[0], Veorinnow1[1], Veorinnow1[2]))
#     zv7_2.rotate(angle=omega * dt, origin=vp.vector(Vorin3[0], Vorin3[1],Vorin3[2]), axis=vp.vector(Veorinnow1[0], Veorinnow1[1], Veorinnow1[2]))
# while(True):
#     # Вращение вокург 5го звена
#     zv5_1.rotate(angle=omega * dt, origin=vp.vector(50/koef, 0, 725/koef), axis=vp.vector(0, 0, 1))
#     zv6_1.rotate(angle=omega * dt, origin=vp.vector(50/koef, 0, 725/koef), axis=vp.vector(0, 0, 1))
#     zv6_2.rotate(angle=omega * dt, origin=vp.vector(50/koef, 0, 725/koef), axis=vp.vector(0, 0, 1))
#     zv7_1.rotate(angle=omega * dt, origin=vp.vector(50/koef, 0, 725/koef), axis=vp.vector(0, 0, 1))
#     zv7_2.rotate(angle=omega * dt, origin=vp.vector(50/koef, 0, 725/koef), axis=vp.vector(0, 0, 1))
# while(True):
#     # Вращение вокург 6го звена
#     zv6_1.rotate(angle=omega * dt, origin=vp.vector(50/koef, 0, 925/koef), axis=vp.vector(1, 0, 0))
#     zv6_2.rotate(angle=omega * dt, origin=vp.vector(50/koef, 0, 925/koef), axis=vp.vector(1, 0, 0))
#     zv7_1.rotate(angle=omega * dt, origin=vp.vector(50/koef, 0, 925/koef), axis=vp.vector(1, 0, 0))
#     zv7_2.rotate(angle=omega * dt, origin=vp.vector(50/koef, 0, 925/koef), axis=vp.vector(1, 0, 0))
# while(True):
#     # Вращение вокург 7го звена
#     zv7_1.rotate(angle=omega * dt, origin=vp.vector(50/koef, 0, 825/koef), axis=vp.vector(0, 0, 1))
#     zv7_2.rotate(angle=omega * dt, origin=vp.vector(50/koef, 0, 825/koef), axis=vp.vector(0, 0, 1))




