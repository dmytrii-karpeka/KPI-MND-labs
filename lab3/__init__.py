from random import randint
import prettytable
import math
import numpy
# Лабораторна робота №2 "ПРОВЕДЕННЯ ДВОФАКТОРНОГО ЕКСПЕРИМЕНТУ З
# ВИКОРИСТАННЯМ ЛІНІЙНОГО РІВНЯННЯ РЕГРЕСІЇ" з предмету МОПЕ
# Варіянт №210 Карпека Дмитрій

#--------------------------------------------------Початкові умови-----------------------------------------------------
variant = 210

y_max = (30 - variant) * 10
y_min = (20 - variant) * 10

#--------------------------------------------------1-------------------------------------------------------------------
# y = b0 + b1X1 + b2X2
#--------------------------------------------------2-------------------------------------------------------------------

x0 = 1
x1_min = -25
x1_max = -5
x2_min = -70
x2_max = -10
m = 6

pt = prettytable.PrettyTable()
pt.field_names = ["X1", "X2"] + ["Y" + str(x) for x in range(1, m+1)]

#--------------------------------------------------3-------------------------------------------------------------------

def matrix_plan(m, ymin, ymax, n=3):
    return [[randint(ymin, ymax) for _ in range(m)] for _ in range(n)]


matrix_y = matrix_plan(m, y_min, y_max)
matrix_y = [[x1_min, x2_min] + matrix_y[0],
            [x1_max, x2_min] + matrix_y[1],
            [x1_min, x2_max] + matrix_y[2]]
pt.add_rows(matrix_y)

print("3) Заповнена матриця планування експерименту для m={0}:".format(m))
print(pt)
print("Значення факторів не нормовані")
#--------------------------------------------------4-------------------------------------------------------------------
#4.1
def average(list):
    return sum(list) / len(list)


y1 = round(average(matrix_y[0][2:len(matrix_y[0])]), 2)
y2 = round(average(matrix_y[1][2:len(matrix_y[1])]), 2)
y3 = round(average(matrix_y[2][2:len(matrix_y[2])]), 2)
print("4) Перевірка дисперсії за критерієм Романовського:")
print("Середнє значення ф-ції відгуку у рядку:\ny1 = {0}\ny2 = {1}\ny3 = {2}".format(y1, y2, y3))
#4.2
def disp(row_y, y_average):
    sigma_squared = 0
    for y in row_y:
        sigma_squared += (y - y_average)**2
    sigma_squared = sigma_squared/len(row_y)
    return sigma_squared


d1 = round(disp(matrix_y[0][2:len(matrix_y[0])], y1), 2)
d2 = round(disp(matrix_y[1][2:len(matrix_y[1])], y2), 2)
d3 = round(disp(matrix_y[2][2:len(matrix_y[2])], y3), 2)
print("Знайдемо дисперсію по рядках:\nD1 = {0}\nD2 = {1}\nD3 = {2}".format(d1, d2, d3))
#4.3
funddev = round(math.sqrt((2*(2*m - 2))/(m*(m-4))), 3)
print("Основне відхилення: ", funddev)
#4.4
def fuvn(dn, dm):
    if dn >= dm:
        fuv = dn/dm
    else:
        fuv = dm/dn
    return fuv


fuv1 = round(fuvn(d1, d2), 3)
fuv2 = round(fuvn(d2, d3), 3)
fuv3 = round(fuvn(d1, d3), 3)
print("Обчислення Fuv:\nFuv1 = {0}\nFuv2 = {1}\nFuv3 = {2}".format(fuv1, fuv2, fuv3))
#4.5
def tetta_uvn(fuv, m=m):
    return ((m - 2)/m * fuv)


tetta_uv1 = round(tetta_uvn(fuv1), 4)
tetta_uv2 = round(tetta_uvn(fuv2), 4)
tetta_uv3 = round(tetta_uvn(fuv3), 4)
print("Обчислення 0uv:\n0uv1 = {0}\n0uv2 = {1}\n0uv3 = {2}".format(tetta_uv1, tetta_uv2, tetta_uv3))
#4.6
def ruvn(tetta, fundev=funddev):
    return abs(tetta - 1)/fundev

ruv1 = ruvn(tetta_uv1)
ruv2 = ruvn(tetta_uv2)
ruv3 = ruvn(tetta_uv3)
print("Обчислення Ruv:\nRuv1 = {0}\nRuv2 = {1}\nRuv3 = {2}".format(ruv1, ruv2, ruv3))
#4.7
def criteria(r1, r2, r3, probability, m=m):
    table = [[0, 2, 6, 8, 10, 12, 15, 20],
             [0.99, 1.73, 2.16, 2.43, 2.62, 2.75, 2.9, 3.08],
             [0.98, 1.72, 2.13, 2.37, 2.54, 2.66, 2.8, 2.96],
             [0.95, 1.71, 2.1, 2.27, 2.41, 2.52, 2.64, 2.78],
             [0.9, 1.69, 2, 2.17, 2.29, 2.39, 2.49, 2.62]]
    rkr = 0
    for i in table[1:len(table)]:
        if i[0] == probability:
            rkr = i[table[0].index(m)]
    if r1 < rkr and r2 < rkr and r3 < rkr:
        print("Дисперсія однорідна")
    else:
        print("Дисперсія неоднорідна")

criteria(ruv1, ruv2, ruv3, 0.9)
#5
x1min_n, x2min_n = -1, -1
x1max_n, x2max_n = 1, 1
matrix_norm = [[x1min_n, x2min_n],
               [x1max_n, x2min_n],
               [x1min_n, x2max_n]]

mx1 = sum([i[0] for i in matrix_norm])/3
mx2 = sum([i[1] for i in matrix_norm])/3
my = sum([y1, y2, y3])/3
a1 = sum(i[0]**2 for i in matrix_norm) / 3
a2 = sum([i[0]*i[1] for i in matrix_norm])/3
a3 = sum(i[1]**2 for i in matrix_norm) / 3

a11 = sum([matrix_norm[0][0]*y1, matrix_norm[1][0]*y2, matrix_norm[2][0]*y3])/3
a22 = sum([matrix_norm[0][1]*y1, matrix_norm[1][1]*y2, matrix_norm[2][1]*y3])/3

b0 = round(numpy.linalg.det([[my, mx1, mx2], [a11, a1, a2], [a22, a2, a3]]) / numpy.linalg.det([[1, mx1, mx2], [mx1, a1, a2], [mx2, a2, a3]]), 2)
b1 = round(numpy.linalg.det([[1, my, mx2], [mx1, a11, a2], [mx2, a22, a3]]) / numpy.linalg.det([[1, mx1, mx2], [mx1, a1, a2], [mx2, a2, a3]]), 2)
b2 = round(numpy.linalg.det([[1, mx1, my], [mx1, a1, a11], [mx2, a2, a22]]) / numpy.linalg.det([[1, mx1, mx2], [mx1, a1, a2], [mx2, a2, a3]]), 2)
print("5) Розрахунок нормованих коефіцієнтів рівняння регресії:\nmx1 = {0}\nmx2 = {1}\nmy = {2}\na1 = {3}\na2 = {4}\na3 = {5}\nb0 = {6}\nb1 = {7}\nb2 = {8}".format(round(mx1, 2), round(mx2, 2), round(my, 2), round(a1, 2), round(a2, 2), round(a3, 2), round(b0, 2), round(b1, 2), round(b2, 2)))
print("Нормоване рівняння регресії:\ny = {0} + {1}*x1 + {2}*x2".format(b0, b1, b2))

y1_check = round(b0 + b1*matrix_norm[0][0] + b2*matrix_norm[0][1], 2)
y2_check = round(b0 + b1*matrix_norm[1][0] + b2*matrix_norm[1][1], 2)
y3_check = round(b0 + b1*matrix_norm[2][0] + b2*matrix_norm[2][1], 2)
print("Перевірка y1:\n{0} + {1}*{2} + {3}*{4} = {5} = {6}".format(b0, b1, matrix_norm[0][0], b2, matrix_norm[0][1], y1_check, y1))
print("Перевірка y2:\n{0} + {1}*{2} + {3}*{4} = {5} = {6}".format(b0, b1, matrix_norm[1][0], b2, matrix_norm[1][1], y2_check, y2))
print("Перевірка y3:\n{0} + {1}*{2} + {3}*{4} = {5} = {6}".format(b0, b1, matrix_norm[2][0], b2, matrix_norm[2][1], y3_check, y3))

#6
#TODO: end 6) and commit lab, pin a link in GT, notify teacher
print("Натуралізація коефіцієнтів:")
delta_x1 = abs(x1_max - x1_min)/2
delta_x2 = abs(x2_max - x2_min)/2
print("delta x1 = {0}\ndelta x2 = {1}".format(delta_x1, delta_x2))
x10 = (x1_max + x1_min)/2
x20 = (x2_max + x2_min)/2
print("x10 = {0}\nx20 = {1}".format(x10, x20))
a0 = round(b0 - b1*x10/delta_x1 - b2*x20/delta_x2, 2)
print("a0 = b0 - b1*x10/delta_x - b2*x20/delta_x2 =", a0)
a1 = round(b1/delta_x1, 3)
a2 = round(b2/delta_x2, 3)
print("a1 = {0}\na2 = {1}".format(a1, a2))
print("Натуралізоване рівняння регресії:\ny = {0} + {1}*x1 + {2}*x2".format(a0, a1, a2))

ych1 = round(a0 + a1*x1_min + a2*x2_min, 2)
ych2 = round(a0 + a1*x1_max + a2*x2_min, 2)
ych3 = round(a0 + a1*x1_min + a2*x2_max, 2)
print("Перевірка по рядках:\ny1 = {0}, y2 = {1}, y3 = {2}".format(ych1, ych2, ych3))
if abs(y1_check - ych1)/ych1 < 0.05 and abs(y2_check - ych2)/ych2 < 0.05 and abs(y3_check - ych3)/ych3 < 0.05:
    print("Коефіцієнти натуралізованого рівняння регресії вірні")