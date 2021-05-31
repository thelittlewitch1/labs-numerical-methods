import math

if __name__ == '__main__':
    eps = 10 ** -9
    i = 0
    x1 = 0
    x2 = 0
    y1 = 1
    y2 = 0
    f1 = x1 * x1 - x2 + 1
    f2 = x1 - math.cos(math.pi * x2 / 2)

    print('I, X1, X2, F1, F2\n')
    print(i, ', ', y1, ', ', y2, ', ', f1, ', ', f2, ';\n')

    while (math.sqrt((y1 - x1) ** 2 + (y2 - x2) ** 2)) >= eps:
        x1 = y1
        x2 = y2

        i = i + 1

        h = (x2 - 1 + x1 ** 2 - 2 * x1 * math.cos(math.pi * x2 / 2)) / (- math.pi * x1 * math.sin(math.pi * x2 / 2) - 1)
        g = - x1 + math.cos(math.pi * x2 / 2) - math.pi * math.sin(math.pi * x2 / 2)*h / 2

        y1 = x1 + g
        y2 = x2 + h

        f1 = y1 ** 2 - y2 + 1
        f2 = y1 - math.cos(math.pi * y2 / 2)
        print(i, ', ', y1, ', ', y2, ', ', f1, ', ', f2, ';\n')
    print('x1 = ', y1, ', x2 = ', y2, ',\nf1 = ', f1, ', f2 = ', f2, ';\n')