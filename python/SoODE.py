import numpy as np
import math as m


class SourceData:
    n = 2
    omega = 25
    a = 2.5 + omega / 40
    u0 = np.array([[0], [-0.412]])
    T = 0.2
    epsAdd = [10 ** -3, 10 ** -5]
    tayMax = 1
    tayMin = 0.1

    def __init__(self):
        self.omega = 25

    def calc_f(self, u, t):
        return np.array([
            -u[0] * u[1] + m.sin(t) / t,
            -u[1] ** 2 + self.a * t / (1 + t ** 2)
        ])

    def calc_df(self, u, t):
        return np.array([
            [-u[1], -u[0]],
            [0, -2 * u[1]]
        ])

    def Newton(self, init, t, flag):
        if flag:
            alfa0 = 1
            alfa1 = 0
            beta = init.tayk
        else:
            alfa0 = ((init.tayk + init.tayk0) ** 2) / (init.tayk0 * (2 * init.tayk + init.tayk0))
            alfa1 = -init.tayk ** 2 / (init.tayk0 * (2 * init.tayk + init.tayk0))
            beta = (init.tayk * (init.tayk + init.tayk0)) / (2 * init.tayk + init.tayk0)
        eps = min(self.epsAdd)
        i = 0
        x = np.array([[0], [0]])
        nx = np.array([[0], [1]])
        y = init.yk
        z = init.yk0
        while (m.sqrt((nx[0] - x[0]) ** 2 + (nx[1] - x[1]) ** 2)) >= eps:
            x = nx
            initF = self.calc_f(self, x, t)
            f = x - alfa1 * z - alfa0 * y - beta * initF
            dF = self.calc_df(self, x, t)

            i = i + 1

            h = (f[1,0] - dF[1, 0] / dF[0, 0] * f[0,0]) / (dF[1, 1] - dF[1, 0] * dF[0, 1] / dF[0, 0])
            g = f[0,0] / dF[0, 0] - dF[0, 1] / dF[0, 0] * h

            nx[0] = x[0] + g
            nx[1] = x[1] + h
        F = np.array(nx)
        return F


class InitialCondition:
    def __init__(self, source):
        self.tk = 0.01
        self.yk = source.u0
        self.yk0 = source.u0
        self.yk1 = source.u0
        self.tayk = source.tayMin
        self.tayk0 = source.tayMin
        self.tayk1 = source.tayMin
        self.eps = source.epsAdd

    def calc_tay_forward(self, source, f):
        res = np.zeros([source.n, 1])
        for i in range(0, source.n):
            # print(source.epsAdd[i], "\n",  f[i], "\n", source.tayMax, "\n", source.epsAdd[i] / (abs(f[i]) + source.epsAdd[i] / source.tayMax))
            res[i,0] = source.epsAdd[i] / (abs(f[i]) + source.epsAdd[i] / source.tayMax)
        return res

    def calc_tay_implicit(self, source):
        tayk = np.zeros(source.n)
        for i in range(source.n):
            if abs(self.eps[i]) > source.epsAdd[i]:
                tayk[i] = self.tayk / 2
            elif (abs(self.eps[i]) > source.epsAdd[i]/4) & (abs(self.eps[i]) <= source.epsAdd[i]):
                tayk[i] = self.tayk
            else:
                tayk[i] = self.tayk * 2
            if tayk[i] > source.tayMax:
                tayk[i] = source.tayMax
        return tayk.min()

    def calc_eps(self):
        res = -(self.tayk / (self.tayk + self.tayk0)) * (
                self.yk1 - self.yk - self.tayk / self.tayk0 * (self.yk - self.yk0))
        return res


def forward_euler_method(s, i):
    print('yk = ', i.yk, '\ntk = ', i.tk, '\n---\n')
    while i.tk < s.T:
        f = s.calc_f(self=s, u=i.yk, t=i.tk)
        i.tayk = min(i.calc_tay_forward(s, f))
        i.yk = i.yk + i.tayk * f
        i.tk = i.tk + i.tayk
        print('tayk = ', i.tayk, '\nyk = ', i.yk, '\ntk = ', i.tk, '\n---\n')


def implicit_euler_method(s, init):
    while init.tk < s.T:
        init.tk1 = init.tk + init.tayk
        init.yk1 = s.Newton(s, init, init.tk, True)
        init.eps = init.calc_eps()
        flag = False
        for j in range(s.n):
            if abs(init.eps[j]) > s.epsAdd[j]:
                init.tayk = init.tayk / 2
                init.tk1 = init.tk
                init.yk1 = init.yk
                flag = True
        if flag:
            continue
        init.tayk1 = init.calc_tay_implicit(s)
        print('tayk = ', init.tayk, '\ny(k+1) = ', init.yk1, '\nt(k+1) = ', init.tk1, '\n---\n')
        init.yk0 = init.yk
        init.yk = init.yk1
        init.tayk0 = init.tayk
        init.tayk = init.tayk1
        init.tk = init.tk1


def shikhman_method(s, i):
    i.yk = s.Newton(s,i, i.tk, True)
    i.yk1 = s.Newton(s, i, i.tk, True)
    while i.tk < s.T:
        i.tk = i.tk + i.tayk
        i.yk0 = i.yk
        i.yk = i.yk1
        i.yk1 = s.Newton(s, i, i.tk, False)
        print('tayk = ', init.tayk, '\ny(k+1) = ', init.yk1, '\nt(k+1) = ', init.tk, '\n---\n')


if __name__ == '__main__':
    sorse = SourceData
    init = InitialCondition(sorse)
    print("Solving of SoODE:\n\tf1(u1,u2,t) = -u1u2 + sin(t)/t\n\tf2(x1,x2,t)=-u2^2 + at/(1+t^2)\n")
    a = int(input("Please, choose the method: 1 - forward euler, 2 - implicit euler, 3 - shikman.\n>"))
    if a == 1:
        print("You've choose forward euler method.\n")
        forward_euler_method(sorse, init)
    elif a == 2:
        print("You've choose implicit euler method.\n")
        implicit_euler_method(sorse, init)
    elif a == 3:
        print("You've choose shikhman method.\n")
        shikhman_method(sorse, init)
    else:
        print("Ding-Dong, you are wrong")
