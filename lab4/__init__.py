from random import randint
import prettytable
import math
from scipy.stats import f, t
from functools import partial

# Лабораторна робота №4 "ПРОВЕДЕННЯ ТРЬОХФАКТОРНОГО ЕКСПЕРИМЕНТУ З ВИКОРИСТАННЯМ
# ЛІНІЙНОГО РІВНЯННЯ РЕГРЕСІЇ З УРАХУВАННЯМ ЕФЕКТУ ВЗАЄМОДІЇ" з предмету МОПЕ
# Варіянт №210 Карпека Дмитрій

#--------------------------------------------------Початкові умови-----------------------------------------------------
variant = 210

min_x = [-20, -30, -30]
max_x = [15, 45, -15]

x0 = [1]
norm_x = [[-1, -1, -1],
        [-1, -1, 1],
        [-1, 1, -1],
        [-1, 1, 1],
        [1, -1, -1],
        [1, -1, 1],
        [1, 1, -1],
        [1, 1, 1]]
natur_x = [ [min_x[0], min_x[1], min_x[2]],
            [min_x[0], min_x[1], max_x[2]],
            [min_x[0], max_x[1], min_x[2]],
            [min_x[0], max_x[1], max_x[2]],
            [max_x[0], min_x[1], min_x[2]],
            [max_x[0], min_x[1], max_x[2]],
            [max_x[0], max_x[1], min_x[2]],
            [max_x[0], max_x[1], max_x[2]] ]

def experiment(m=3, n=8):

    regression_str = 'y = {} + {} * x1 + {} * x2 + {} * x3 + {} * x1x2 + {} * x1x3 + {} * x2x3 + {} * x1x2x3'

    def matrix_plan(m, ymin, ymax, n):
        return [[randint(ymin, ymax) for _ in range(m)] for _ in range(n)]

    def multiplication(a, b):
        return a * b

    def average(list):
        return sum(list) / len(list)

    def round_to_2(number):
        return round(number, 2)

    def dispersion(list_y, aver_list_y):
        return [round_to_2(average(list(map(lambda y: (y - aver_list_y[i]) ** 2, list_y[i])))) for i in range(len(list_y))]

    def cochrane_criteria(S_y):
        global m
        print("\nКритерій Кохрена\n")
        Gp = max(S_y) / sum(S_y)
        q = 0.05
        q_ = q / f2
        chr = f.ppf(q=1 - q_, dfn=f1, dfd=(f2 - 1) * f1)
        Gt = chr / (chr + f2 - 1)
        print("Тест Кохрена: Gr = " + str(round(Gp, 3)))
        if Gp < Gt:
            print("Дисперсії однорідні з імовірністю 0.95")
            pass
        else:
            print("\nДисперсії неоднорідні.\nПовтор експерименту для m + 1\n")
            m = m + 1
            experiment(m)

    def student_criteria(S_y, d):
        print("\nКритерій Ст'юдента\n")
        bettaList = [sum(S_y) * x0[0] / n,
                     average(list(map(multiplication, S_y, x1i))),
                     average(list(map(multiplication, S_y, x2i))),
                     average(list(map(multiplication, S_y, x3i))),
                     average(list(map(multiplication, S_y, norm_x12))),
                     average(list(map(multiplication, S_y, norm_x13))),
                     average(list(map(multiplication, S_y, norm_x23))),
                     average(list(map(multiplication, S_y, norm_x123)))]
        bettaList = [round_to_2(i) for i in bettaList]

        list_t = [bettaList[i] * S for i in range(n)]

        for i in range(n):
            if list_t[i] < t.ppf(q=0.975, df=f3):
                list_b[i] = 0
                d -= 1
                print('Коефіцієнт b' + str(i) + ' незначимий, тому виключається із рівняння регресії')
        print("\nСкореговане рівняння регресії:")
        print(regression_str.format(*map(round_to_2, list_b)))

    def fisher_criteria(d):
        global m
        print("\nКритерій Фішера\n")
        f4 = n - d
        S_ad = (m * sum(
            [(list_b[0] + list_b[1] * x1i[i] + list_b[2] * x2i[i] + list_b[3] * x3i[i] + list_b[4] * norm_x12[i] +
              list_b[5] * norm_x13[i] + list_b[6] * norm_x23[i] + list_b[7] * norm_x123[i]
              - average_list_y[i]) ** 2 for i in range(n)]) / f4)
        Fp = S_ad / Sb

        if Fp > f.ppf(q=0.95, dfn=f4, dfd=f3):  # перевірка критерію Фішера з використанням scipy
            print('Математична модель неадекватна експериментальним даним на рівні значимості 0.05.\nПовтор експерименту для m+1')
            m = m + 1
            experiment(m)
        else:
            print('Математична модель адекватна експериментальним даним на рівні значущості 0.05')

    def printed_matrixes():
        pt1 = prettytable.PrettyTable()
        pt2 = prettytable.PrettyTable()
        pt1.field_names = ["X0", "X1", "X2", "X3"] + ["X12", "X13", "X23", "X123"] + ["Y" + str(x) for x in range(1, m + 1)] + ["Aver Y"] + ["S_y"]
        pt2.field_names = ["X0", "X1", "X2", "X3"] + ["X12", "X13", "X23", "X123"] + ["Y" + str(x) for x in range(1, m + 1)] + ["Aver Y"] + ["S_y"]

        print("Матриця повного факторного експерименту з натуралізованими значеннями:\n")
        pt1.add_rows([x0 + natur_x[i] + natur_x12[i] + natur_x13[i] + natur_x23[i] + natur_x123[i] + matrix_y[i] + [average_list_y[i]] + [S_y[i]] for i in range(n)])
        print(pt1)

        print("\nМатриця повного факторного експерименту з нормалізованими значеннями:\n")
        pt2.add_rows([x0 + norm_x[i] + [norm_x12[i]] + [norm_x13[i]] + [norm_x23[i]] + [norm_x123[i]] + matrix_y[i] + [average_list_y[i]] + [S_y[i]] for i in range(n)])
        print(pt2)

    m = m
    n = 8

    x_average_max = sum(max_x) / 3
    x_average_min = sum(min_x) / 3

    y_max = round(200 + x_average_max)
    y_min = round(200 + x_average_min)

    matrix_y = matrix_plan(m, y_min, y_max, n)
    average_list_y = [round(average(matrix_y[i]), 2) for i in range(len(matrix_y))]

    S_y = dispersion(matrix_y, average_list_y)

    f1 = m - 1
    f2 = n
    f3 = f1 * f2
    d = 4

    Sb = sum(S_y) / n
    S = math.sqrt(Sb / (n * m))

    norm_x12 = [norm_x[i][0] * norm_x[i][1] for i in range(len(norm_x))]
    norm_x13 = [norm_x[i][0] * norm_x[i][2] for i in range(len(norm_x))]
    norm_x23 = [norm_x[i][1] * norm_x[i][2] for i in range(len(norm_x))]
    norm_x123 = [norm_x[i][0] * norm_x[i][1] * norm_x[i][2] for i in range(len(norm_x))]

    natur_x12 = [[natur_x[i][0] * natur_x[i][1]] for i in range(len(natur_x))]
    natur_x13 = [[natur_x[i][0] * natur_x[i][2]] for i in range(len(natur_x))]
    natur_x23 = [[natur_x[i][1] * natur_x[i][2]] for i in range(len(natur_x))]
    natur_x123 = [[natur_x[i][0] * natur_x[i][1] * natur_x[i][2]] for i in range(len(natur_x))]

    x1i = [norm_x[i][0] for i in range(n)]
    x2i = [norm_x[i][1] for i in range(n)]
    x3i = [norm_x[i][2] for i in range(n)]

    list_b = [0] * n # b0, b1, b2, b3, b12, b13, b23, b123
    list_b[0] = average(average_list_y)
    list_b[1] = average([average_list_y[i] * x1i[i] for i in range(n)])
    list_b[2] = average([average_list_y[i] * x2i[i] for i in range(n)])
    list_b[3] = average([average_list_y[i] * x3i[i] for i in range(n)])
    list_b[4] = average([average_list_y[i] * x1i[i] * x2i[i] for i in range(n)])
    list_b[5] = average([average_list_y[i] * x1i[i] * x3i[i] for i in range(n)])
    list_b[6] = average([average_list_y[i] * x2i[i] * x3i[i] for i in range(n)])
    list_b[7] = average([average_list_y[i] * x1i[i]* x2i[i] * x3i[i] for i in range(n)])

    printed_matrixes()
    print("\nРівняння\n" + regression_str.format(*map(round_to_2, list_b)))
    cochrane_criteria(S_y)
    student_criteria(S_y, d)
    fisher_criteria(d)
#--------------------
m = 3
experiment(m)
