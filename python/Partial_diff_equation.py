import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def implicit():
    c = 0.5
    T = 6
    Tau = 0.01
    Lx = 5
    Ly = 3
    hx = 0.5
    hy = 0.5
    nx = int(Lx / hx)
    ny = int(Ly / hy)
    U0 = 30
    t = 0
    Ucur = np.ones(shape=(nx, ny)) * U0
    Tx0 = 30
    Ty0 = 30
    TxL = 100
    TyL = 30
    Unext = Ucur
    bounds = np.zeros((nx + 2, ny + 2))
    bounds[nx + 1, :] = 30
    bounds[0, :] = 30
    bounds[:, 0] = 30
    bounds[:, ny + 1] = 100
    bounds[1:nx + 1, 1:ny + 1] = Ucur
    Ax = np.ones(nx + 2) * (c / (hx ** 2))
    Bx = np.ones(nx + 2) * (2 * c / (hx ** 2) + 1 / Tau)
    Cx = np.ones(nx + 2) * (c / (hx ** 2))
    Ay = np.ones(nx + 2) * (c / (hy ** 2))
    By = np.ones(nx + 2) * (2 * c / (hy ** 2) + 1 / Tau)
    Cy = np.ones(nx + 2) * (c / (hy ** 2))
    alfx = np.zeros(nx + 2)
    alfx[0] = 0
    bettax = np.zeros(nx + 2)
    bettax[0] = Tx0
    alfy = np.zeros(ny + 2)
    alfy[0] = 0
    bettay = np.zeros(ny + 2)
    bettay[0] = 30
    NewT = np.zeros(nx + 2)
    NewTy = np.zeros(ny + 2)
    Fx = np.zeros(nx + 2)
    Fy = np.zeros(ny + 2)
    Newbounds = bounds
    while t < T:
        for i in range(1, nx):
            for j in range(0, ny + 1):
                Fy[j] = -bounds[j, i-1] / Tau
            for j in range(1, ny + 1):
                alfy[j] = Ay[j] / (By[j] - Cy[j] * alfy[j - 1])
                bettay[j] = (Cy[j] * bettay[j - 1] - Fy[j]) / (By[j] - Cy[j] * alfy[j - 1])
            NewTy[ny + 1] = TyL
            for j in range(0, ny):
                NewTy[ny + 1 - j] = alfy[ny + 1 - j] * NewTy[ny + 1 - j] + bettay[ny + 1 - j]
            for k in range(1, ny):
                Newbounds[k, i-1] = NewT[k]

        for i in range(1, ny):
            for j in range(0, nx + 1):
                Fx[j] = -Newbounds[i, j] / Tau - f(t)
            for j in range(1, nx + 1):
                alfx[j] = Ax[j] / (Bx[j] - Cx[j] * alfx[j - 1])
                bettax[j] = (Cx[j] * bettax[j - 1] - Fx[j]) / (Bx[j] - Cx[j] * alfx[j - 1])
            NewT[nx + 1] = TxL
            for j in range(0, nx):
                NewT[nx + 2 - j] = alfx[nx + 2 - j] * NewT[nx + 3 - j] + bettax[nx + 2 - j]
            bounds[i, 1:nx + 1] = NewT[1, nx + 1]
        plot(Ucur, t, Tau, T)
        t = t + Tau
    plot(Ucur, t, Tau, T)


def explicit():
    c = 0.5
    T = 6
    Tau = 0.01
    Lx = 5
    Ly = 3
    h = 0.5
    nx = int(Lx / h)
    ny = int(Ly / h)
    T0 = 30
    t = 0
    if c*Tau/(h**2)<1/2:
        Ucur = np.ones(shape=(nx, ny)) * T0
        Unext = Ucur
        bounds = np.zeros((nx + 2, ny + 2))
        bounds[nx + 1, :] = 30
        bounds[0, :] = 30
        bounds[:, 0] = 30
        bounds[:, ny + 1] = 100
        bounds[1:nx + 1, 1:ny + 1] = Ucur
        while t < T:
            for i in range(nx):
                for j in range(ny):
                    I = i + 1
                    J = j + 1
                    Unext[i, j] = Tau * c * \
                                  (1 / (h ** 2) * (bounds[I, J + 1] - 2 * bounds[I, J] + bounds[I, J - 1]) + 1 / (h ** 2) *
                                   (bounds[I + 1, J] - 2 * bounds[I, J] + bounds[I - 1, J])) + Tau * f(t) + Unext[i, j]
            bounds[1:nx + 1, 1:ny + 1] = Unext
            plot(Ucur, t, Tau, T)
            t = t + Tau
        plot(Ucur, t, Tau, T)
    else:
        print("Something went wrong.")


def plot(Ucur, t, Tau, T):
    for i in range(T):
        if ((t - i <= Tau) and (t - i > 0)) or (t-T >= 0):
            ax = sns.heatmap(Ucur, vmin=25, vmax=105)
            ax.set_title("heatmap at t = %.2f" % t)
            ax.set(xlabel='X axe', ylabel='Y axe')
            # plt.savefig("plots/heatmap_at_" + str(t) + ".png")
            plt.show()
            plt.clf()


def f(t):
    return t ** 2 + 2 * t + 1


if __name__ == '__main__':
    print("Solving of partial differential equation:")
    # a = int(input("Please, choose the method: 1 - explicit scheme, 2 - implicit scheme.\n>"))
    a = 1
    if a == 1:
        print("You've choose explicit scheme.\n")
        explicit()
    elif a == 2:
        print("You've choose implicit scheme.\n")
        implicit()
    else:
        print("Ding-Dong, you are wrong")