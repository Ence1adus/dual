import math
import numpy as np
import cmath
from sys import stdin
import sympy
from copy import deepcopy

class DualNumber():
    """Реализуем класс  дуальных чисел"""
    def __init__(self,real,dual):
        """Зададим само дуальное число"""
        self.real = real
        self.dual = dual

    def __add__(self, other):
        return DualNumber(self.real + other.real , self.dual + other.dual)

    def __sub__(self, other):
        return DualNumber(self.real - other.real, self.dual - other.dual)

    def __mul__(self, other):
        return DualNumber(self.real * other.real, self.dual * other.real + self.real * other.dual)

    def __truediv__(self, other):
        return DualNumber(self.real / other.real,(self.dual * other.real - self.real * other.dual) / (other.real * other.real))

    def __str__(self):
        if self.dual == 0:
            result = "%.2f+0.00\u03B5" % (self.real)
        elif self.real == 0:
            if self.dual >= 0:
                result = "0.00+%.2f\u03B5" % (self.dual)
            else:
                result = "0.00-%.2f\u03B5" % (abs(self.dual))
        elif self.dual > 0:
            result = "%.2f+%.2f\u03B5" % (self.real, self.dual)
        else:
            result = "%.2f-%.2f\u03B5" % (self.real, abs(self.dual))
        return result

    def Sin(self):
        """Возвращает значение дуального Синуса"""
        return DualNumber(math.sin(self.real), self.dual * math.cos(self.real))
    def Cos(self):
        """Возвращает значение дуального Косинуса"""
        return DualNumber(math.cos(self.real), - self.dual * math.sin(self.real))

    def InfoAboutNumber(self):
        """Возвращает значение Реальной и Дуальной части числа"""
        print("Реальная часть числа = " + str(self.real))
        print("Дуальная часть числа = " + str(self.dual))

class DualMatrix(DualNumber):
    """Реализуем класс  дуальных Матриц"""

    def __init__(self, list_of_lists):
        self.matrix = deepcopy(list_of_lists)


    def __str__(self):
        mystring = '\n'.join('\t'.join(map(str, row))
                         for row in self.matrix)
        mystring += '\n//////////////////////////////////'
        return mystring

    def size(self):
        sizepair = (len(self.matrix), len(self.matrix[0]))
        return sizepair


    def __getitem__(self, idx):
        return self.matrix[idx]

    def __add__(self, other):
        """Производит суммирование матриц с алгеброй деальных чисел"""
        # other = DualMatrix(other)
        sizeFirstMatrix = self.size()
        LengthString = sizeFirstMatrix[0]
        LengthColumn = sizeFirstMatrix[1]
        sizeSecondMatrix = other.size()
        LengthColumnSM = sizeSecondMatrix[1]
        result_matrix = [[0 for i in range(LengthColumn)] for i in range(LengthString)]
        for i in range(LengthString):
            for j in range(LengthColumn):
                element = DualNumber(0, 0)
                result_matrix[i][j] = element
                summa = DualNumber(0, 0)
                summa.real = self.matrix[i][j].real + other[i][j].real
                summa.dual = self.matrix[i][j].dual + other[i][j].dual
                result_matrix[i][j].real = result_matrix[i][j].real + summa.real
                result_matrix[i][j].dual = result_matrix[i][j].dual + summa.dual
        return DualMatrix(result_matrix)

    def __sub__(self, other):
        """Производит вычитание матриц с алгеброй деальных чисел"""
        # other = DualMatrix(other)
        sizeFirstMatrix = self.size()
        LengthString = sizeFirstMatrix[0]
        LengthColumn = sizeFirstMatrix[1]
        sizeSecondMatrix = other.size()
        LengthColumnSM = sizeSecondMatrix[1]
        result_matrix = [[0 for i in range(LengthColumn)] for i in range(LengthString)]
        for i in range(LengthString):
            for j in range(LengthColumn):
                element = DualNumber(0, 0)
                result_matrix[i][j] = element
                summa = DualNumber(0, 0)
                summa.real = self.matrix[i][j].real - other[i][j].real
                summa.dual = self.matrix[i][j].dual - other[i][j].dual
                result_matrix[i][j].real = result_matrix[i][j].real + summa.real
                result_matrix[i][j].dual = result_matrix[i][j].dual + summa.dual
        return DualMatrix(result_matrix)

    def __mul__(self, other):
        """Производит умножение матриц с алгеброй деальных чисел"""
        #other = DualMatrix(other)
        sizeFirstMatrix = self.size()
        LengthString = sizeFirstMatrix[0]
        LengthColumn = sizeFirstMatrix[1]
        sizeSecondMatrix = other.size()
        LengthColumnSM = sizeSecondMatrix[1]

        result_matrix = [[0 for i in range(LengthColumnSM)] for i in range(LengthString)]
        for i in range(LengthString):
            for j in range(LengthColumnSM):
                element = DualNumber(0, 0)
                result_matrix[i][j] = element
                for k in range(LengthColumn):
                    element = DualNumber(0, 0)
                    element.real = self.matrix[i][k].real * other[k][j].real
                    element.dual = self.matrix[i][k].dual * other[k][j].real + self.matrix[i][k].real * other[k][j].dual
                    result_matrix[i][j].real = result_matrix[i][j].real + element.real
                    result_matrix[i][j].dual = result_matrix[i][j].dual + element.dual
        return DualMatrix(result_matrix)

    def DTM(self,per):
        """Возвращает даульную матрицу в разложенном виде"""
        sizeFirstMatrix = self.size()
        LengthString = sizeFirstMatrix[0]
        LengthColumn = sizeFirstMatrix[1]

        result_matrix1 = [[0 for i in range(LengthColumn)] for i in range(LengthString)]
        result_matrix2 = [[0 for i in range(LengthColumn)] for i in range(LengthString)]
        for i in range(LengthString):
            for j in range(LengthColumn):
                result_matrix1[i][j] = self[i][j].real
                result_matrix2[i][j]= self[i][j].dual
        if per :
            print("Действительная часть")
            for i in range(0, len(result_matrix1)):
                for i2 in range(0, len(result_matrix1[i])):
                    print("{0:.6f}".format(result_matrix1[i][i2]), end=' ')
                print()
            print("Дуальная часть")
            for i in range(0, len(result_matrix2)):
                for i2 in range(0, len(result_matrix2[i])):
                    print("{0:.6f}".format(result_matrix2[i][i2]), end=' ')
                print()
        return result_matrix1,result_matrix2
    def eye(a):
        result_matrix = [[0 for i in range(a)] for i in range(a)]
        for i in range(a):
            for j in range(a):
                element = DualNumber(0, 0)
                DualNum = DualNumber(1, 0)
                result_matrix[i][j] = element
                if i==j:
                    result_matrix[i][j] = DualNum
        return DualMatrix(result_matrix)
    def Transpose(self):
        """Трансопирует матрицу"""
        sizeFirstMatrix = self.size()
        LengthString = sizeFirstMatrix[0]
        LengthColumn = sizeFirstMatrix[1]
        result_matrix = [[0 for i in range(LengthString)] for i in range(LengthColumn)]
        for i in range(LengthColumn):
            for j in range(LengthString):
                element = DualNumber(0, 0)
                result_matrix[i][j] = element
                result_matrix[i][j].real = self[j][i].real
                result_matrix[i][j].dual = self[j][i].dual

        return DualMatrix(result_matrix)
    def CM(self,a,b,c,d):
        """Введите сначала параметры выреза по строке а потом по стоблцу"""
        LengthString = b-a+1
        LengthColumn = d-c+1
        result_matrix = [[0 for i in range(LengthString)] for i in range(LengthColumn)]
        for i in range(LengthColumn):
            for j in range(LengthString):
                element = DualNumber(0, 0)
                result_matrix[i][j] = element
                result_matrix[i][j].real = self[a+i][j+c].real
                result_matrix[i][j].dual =self[a+i][j+c].dual
        return DualMatrix(result_matrix)

    def CrMaFrCol (Jac):
        """Возвращает матрицу на основании совокупности столбцов"""
        LenghtMassive = len(Jac)
        LenghtElement = (Jac[0].size())[0]
        result_matrix = [[0 for i in range(LenghtMassive)] for i in range(LenghtElement)]
        for i in range(LenghtElement):
            for j in range(LenghtMassive):
                element = DualNumber(0, 0)
                result_matrix[i][j] = element
                result_matrix[i][j].real = ((Jac[j])[i][0]).real
                result_matrix[i][j].dual = ((Jac[j])[i][0]).dual
        return DualMatrix(result_matrix)





