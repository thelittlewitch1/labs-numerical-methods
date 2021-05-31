import numpy as np
import math as m


class SourceData:
    def __init__(self):
        self.n = 2
        omega = 25
        self.a = 2.5 + omega / 40
        self.u0 = np.array([
            [0], [-0.412]
        ])
        self.T = 1
        self.epsAdd = [10 ** -3, 10 ** -5]
        self.tayMax = 1
        self.tayMin = 0.1

    def calc_f(self, u, t):
        return np.array([
            [-u[1] * u[2] + m.sin(t) / t],
            [-u[2] ** 2 + self.a * t / (1 + t ** 2)]
        ])

    # TODO: дописать
    def Newton(self):
        eps = min(self.epsAdd)
        i = 0
        x1 = 0
        x2 = 0
        y1 = 1
        y2 = 0
        f = np.array([
            []
        ])
        while (m.sqrt((y1 - x1) ** 2 + (y2 - x2) ** 2)) >= eps:
            x1 = y1
            x2 = y2

            i = i + 1

            # h =
            # g =

            #y1 = x1 + g
            #y2 = x2 + h

            # f =  np.array([
            #        []
            # ])
        return f


class InitialCondition:
    def __init__(self, source):
        self.tk = 0
        self.yk = source.u0
        self.yk0 = source.u0
        self.yk1 = source.u0
        self.tayk = source.taymin
        self.tayk0 = source.taymin
        self.tayk1 = source.taymin
        self.eps = source.epsadd

    def calc_tay(self, source, f):
        res = np.zeros(source.n, 1)
        for i in range(0, source.n):
            res[i] = source.epsAdd[i] / (abs(f[i]) + source.epsAdd[i] / source.tayMax)
        return res

    def calc_tay(self, source):
        tayk = np.array([np.zeros(source.n)])
        for i in range(source.n):
            if abs(self.eps(i)) > source.epsAdd(i):
                tayk[i] = self.tayk / 2
            elif (abs(self.eps[i]) > source.epsAdd / 4) & (abs(self.eps[i]) <= source.epsAdd[i]):
                tayk[i] = self.tayk
            else:
                tayk[i] = self.tayk * 2

    def calc_eps(self):
        return -(self.tayk / (self.tayk + self.tayk0)) * (
                self.yk1 - self.yk - self.tayk / self.tayk0 * (self.yk - self.yk0))


def forward_euler_method(s, i):
    print('yk = ', i.yk, 'tk = ', i.tk, '\n---\n')
    while i.tk < s.T:
        f = s.calc_f(i.yk, i.tk)
        i.tayk = min(i.calc_tay(s, f))
        i.yk = i.yk + i.tayk * f
        i.tk = i.tk + i.tayk
        print('tayk = ', i.tayk, 'yk = ', i.yk, 'tk = ', i.tk, '\n---\n')


def implicit_euler_method(s, i):
    while i.tk < s.T:
        i.tk1 = i.tk + i.tayk
        i.yk1 = s.Newton()
        i.eps = i.calc_eps
        flag = False
        for i in range(s.n):
            if abs(i.eps(i)) > s.epsAdd(i):
                i.tayk = i.tayk / 2
                i.tk1 = i.tk
                i.yk1 = i.yk
                flag = True
        if flag:
            continue
        i.tayk1 = i.calc_tay(s)
        if i.tayk1 > s.tayMax:
            i.tayk1 = s.tayMax
        print(i.yk1, i.tk1)
        i.yk0 = i.yk
        i.yk = i.yk1
        i.tayk0 = i.tayk
        i.tayk = i.tayk1
        i.tk = i.tk1


def shikhman_method (s, i):
    return 0;