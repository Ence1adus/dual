import math
import numpy as np
import cmath
from sys import stdin
from copy import deepcopy
import DualNumberModule
import sympy
from DualNumberModule import DualNumber
from DualNumberModule import DualMatrix

def RMx(d,theta):
    Rx = DualMatrix([[DualNumber(1, 0), DualNumber(0, 0), DualNumber(0, 0)],
                [DualNumber(0, 0), DualNumber.Cos(DualNumber(theta, d)),
                 DualNumber(-1, 0) * DualNumber.Sin(DualNumber(theta, d))],
                [DualNumber(0, 0), DualNumber.Sin(DualNumber(theta, d)),
                 DualNumber.Cos(DualNumber(theta, d))]])
    return Rx

def RMz(d,theta):
    Rz = DualMatrix([[DualNumber.Cos(DualNumber(theta, d)),
                        DualNumber(-1, 0) * DualNumber.Sin(DualNumber(theta, d)), DualNumber(0, 0)],
                       [DualNumber.Sin(DualNumber(theta, d)), DualNumber.Cos(DualNumber(theta, d)),
                        DualNumber(0, 0)],
                       [DualNumber(0, 0), DualNumber(0, 0), DualNumber(1, 0)]])
    return Rz

def RMy(d,theta):
    Ry = DualMatrix([[DualNumber.Cos(DualNumber(theta, d)),
                      DualNumber(0, 0),DualNumber.Sin(DualNumber(theta, d))],
                       [DualNumber(0, 0), DualNumber(1, 0),DualNumber(0, 0)],
                       [DualNumber(-1, 0) * DualNumber.Sin(DualNumber(theta, d)),
                        DualNumber(0, 0),DualNumber.Cos(DualNumber(theta, d))]])
    return Ry