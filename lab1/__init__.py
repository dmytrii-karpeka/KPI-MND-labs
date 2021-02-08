from random import randint
import numpy as np

# Лабораторна робота №1 "ЗАГАЛЬНІ ПРИНЦИПИ ОРГАНІЗАЦІЇ ЕКСПЕРИМЕНТІВ З ДОВІЛЬНИМИ ЗНАЧЕННЯМИ ФАКТОРІВ" з предмету МОПЕ
# Варіянт №210 Карпека Дмитрій - min(Y)
#-------------------------------1--------------------------------------------------------------------------------------
max_point = 20
matrix_plan = [[randint(0, max_point) for j in range(8)] for i in range(3)]
print("Генерація значення факторів у точках експерименту:\n", np.matrix(matrix_plan))
#-------------------------------2--------------------------------------------------------------------------------------
a0 = 0.2
a1 = 1.3
a2 = 1.7
a3 = 6


def reaction(a0, a1, a2, a3, x1, x2, x3):
    y = a0 + a1*x1 + a2*x2 + a3*x3
    return y


y_list = []
for i in range(8):
    vector = []
    for j in range(3):
        vector.append(matrix_plan[j][i])
    y_list.append(round((reaction(a0, a1, a2, a3, vector[0], vector[1], vector[2])), 3))
print("Обчислення функції відгуку в кожній точці:\n", y_list)
#-------------------------------3--------------------------------------------------------------------------------------
# Пошук центру експерименту
x0_list = []
for i in range(3):
    x0 = (min(matrix_plan[i]) + max(matrix_plan[i]))/2
    x0_list.append(x0)
print("Центри експерименту:\n", x0_list)
# Пошук інтервала зміни фактора
dx_list = []
for i in range(3):
    dx = x0_list[i] - min(matrix_plan[i])
    dx_list.append(dx)
print("Інтервали зміни фактора:\n", dx_list)
# Пошук нормованого значення Xн
xn_matrix = []
for i in range(3):
    xn = []
    for j in range(8):
        xn.append(round( ( (matrix_plan[i][j] - x0_list[i]) / dx_list[i]), 3))
    xn_matrix.append(xn)
print("Матриця нормованих значень Xн:\n", np.matrix(xn_matrix))
# Пошук Yет
yet = reaction(a0, a1, a2, a3, x0_list[0], x0_list[1], x0_list[2])
print("Значення функції відгуку при нульових рівнях факторів:\n", yet)
#-------------------------------4--------------------------------------------------------------------------------------
answer = [matrix_plan[i][y_list.index(min(y_list))] for i in range(3)]
print("Точка плану, що задовільняє критерій вибору оптимальности:\n", answer)
