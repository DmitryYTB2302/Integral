import math
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

# нижняя и верхняя граница
a = 0
b = 2.5

# исходная функция интеграла
def f(x):
    return (1 + 2 * x) / ((2 + x) * np.sqrt(2 + x))

# метод Ньютона - Лейбница
def newton_Liebnitz(a, b):
    return (14 + 4 * b) / np.sqrt(2 + b) - (14 + 4 * a) / np.sqrt(2 + a)

# метод Ньютона - Котеса для интегрирования
def newton_Kotes(a, b, k):
    # задаем коэффициенты Котеса
    H = [0.0659722, 0.2604166, 0.1736111, 0.1736111, 0.2604166, 0.0659722]
    # количество узлов
    n = 5
    # вычисление шага интегрирования
    h = float(b - a) / float(k * n)
    sum = 0
    i = a
    # запускаем цикл по всем отрезкам
    while (i < b - (n - 2) * h):
        # запускаем цикл по каждму весу
        for j in range(0, n + 1):
            # считаем интеграл по формуле Котеса
            sum += ((H[j]) * f(i + (j * h)))
        # прохидмся по всем отрезам
        i += h * n
    return sum * h * n

# метод Гаусса
def gauss(a, b, k):
    # весы и узлы для n = 5
    X = [0.0469101, 0.2307653, 0.5, 0.7692347, 0.9530899]
    A = [0.1184634, 0.2393143, 0.2844444, 0.2393143, 0.1184634]
    # количество узлов
    n = 5
    # вычисление шага интегрирования
    h = float(b - a) / float(k)
    sum = 0
    i = a
    # глобальный цикл для вычислений
    while (i < b - (h / 2)):
        # цикл по весам
        for j in range(0, n):
            # считаем интеграл по формуле Гаусса
            sum += ((A[j]) * f(i + (h * X[j])))
        # проходимся по всем отрезкам
        i += h
    return sum * h

def newton_Kotes_test():
    print("Метод Ньютона - Котесса")
    for i in range(1, 11):
        print(f"Число шагов {i}: ", newton_Kotes(a, b, i))
    return ''

def gauss_test():
    print("Метод Гаусса")
    for i in range(1, 11):
        print(f"Число шагов {i}: ", gauss(a, b, i))
    return ''

def chart1():
    fig, ax = plt.subplots()
    ax.set_title('График на подынтегральном промежутке')
    ax.set_ylabel('Ось y')
    ax.set_xlabel('Ось x')
    x = np.arange(0, 2.5, 0.01)
    y = f(x)
    ax.plot(x, y, 'r-')
    ax.grid()
    plt.ylim([0.35, 0.650])
    plt.xlim([0, 2.5])
    plt.show()

def chart2():
    exact_value = newton_Liebnitz(a, b)
    steps = list(range(1, 21))
    errors_Kotes = []
    errors_Gauss = []
    for i in steps:
        approx_value_Kotes = newton_Kotes(a, b, i)
        approx_value_Gauss = gauss(a, b, i)
        error_Kotes = abs(exact_value - approx_value_Kotes)
        error_Gauss = abs(exact_value - approx_value_Gauss)
        errors_Kotes.append(error_Kotes)
        errors_Gauss.append(error_Gauss)
    plt.plot(steps, errors_Kotes, '-r', label='Ошибка метода Ньютона - Котеса')
    plt.plot(steps, errors_Gauss, '-b', label='Ошибка метода Гаусса')
    plt.xlabel('Количество шагов интегрирования')
    plt.ylabel('Фактическая ошибка')
    plt.title('Зависимость ошибки от числа шагов')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.show()

def chart3():
    exact_value = newton_Liebnitz(a, b)
    steps = list(range(1, 21))
    errors_Kotes = []
    errors_Gauss = []
    h_values = []
    for i in steps:
        h = (b - a) / i
        h_values.append(h)
        approx_value_Kotes = newton_Kotes(a, b, i)
        approx_value_Gauss = gauss(a, b, i)
        error_Kotes = abs(exact_value - approx_value_Kotes)
        error_Gauss = abs(exact_value - approx_value_Gauss)
        errors_Kotes.append(error_Kotes)
        errors_Gauss.append(error_Gauss)
    plt.plot(h_values, errors_Kotes, '-r', label='Ошибка метода Ньютона - Котеса')
    plt.plot(h_values, errors_Gauss, '--b', label='Ошибка метода Гаусса')
    plt.xlabel('Размер шага')
    plt.ylabel('Фактическая ошибка')
    plt.title('Зависимость ошибки от размера шага')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.show()

def chart4():
    exact_value = newton_Liebnitz(a, b)
    steps = list(range(1, 21))
    values_Kotes = []
    values_Gauss = []
    for i in steps:
        approx_value_Kotes = newton_Kotes(a, b, i)
        approx_value_Gauss = gauss(a, b, i)
        values_Kotes.append(approx_value_Kotes)
        values_Gauss.append(approx_value_Gauss)
    plt.plot(steps, values_Kotes, '-g', label='Значение метода Ньютона - Котеса')
    plt.plot(steps, values_Gauss, '--b', label='Значение метода Гаусса')
    plt.plot(steps, [exact_value] * len(steps), 'r--', label='Точное значение')
    plt.xlabel('Количество шагов интегрирования')
    plt.ylabel('Значение интеграла')
    #plt.title('Зависимость значения интеграла от числа шагов')
    plt.legend()
    plt.grid()
    plt.show()

#newton_Kotes_test()
#gauss_test()

#chart1()
#chart2()
chart3()
#chart4()