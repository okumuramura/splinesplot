import matplotlib.pyplot as plt
import numpy as np

import random


def power_list(x, n):
    t = 1
    for i in range(n):
        yield t
        t = t * x


def line_spline(x, y):
    output_formuls = []
    for i in range(len(x) - 1):
        f = lambda _x, y = y, x = x, i = i : y[i] + (_x - x[i]) * (y[i+1] - y[i]) / (x[i + 1] - x[i])
        output_formuls.append(f)
    return output_formuls

def cube_spline(ax, ay):
    # si     = ai + bix + cix^2 + dix^3
    # si'    = bi + 2cix + 3dix^2
    # si''   = 2ci + 6dix
    power = 3
    n = len(ax)
    v = np.zeros(4*(n - 1)) # vector of solutions
    mat: np.ndarray = np.zeros((4*(n - 1), 4*(n - 1))) # matrix of values
    spline_num = len(ax) - 1 # num of splines
    # f
    # si(xi) = yi and s(i+1)(x(i+1)) = y(i+1)
    for i in range(spline_num):
        # point 1
        v[i*2] = y[i]
        spos = 4 * i
        mat[i*2][spos:spos+4] = [t for t in power_list(x[i], 4)] # generate powers of x[i]
        # point 2
        v[i*2+1] = y[i+1]
        mat[i*2+1][spos:spos+4] = [t for t in power_list(x[i+1], 4)] # generate powers of x[i+1]

    # f'
    adjust = spline_num * 2
    for i in range(1, spline_num):
        v[adjust+i-1] = 0
        spos = 4 * (i - 1)
        xx = x[i] * x[i]
        mat[adjust+i-1][spos:spos+8] = [0, 1, 2 * x[i], 3 * xx, 0, -1, -2 * x[i], -3 * xx]
    # f''
    adjust = adjust + spline_num - 1
    for i in range(1, spline_num):
        v[adjust+i-1] = 0
        spos = 4 * (i - 1)
        mat[adjust+i-1][spos:spos+8] = [0, 0, 2, 6*x[i], 0, 0, -2, -6 * x[i]]

    v[-1] = 0
    mat[-1][0:4] = [0, 0, 2, 6 * x[0]] # f''(x0) = 0
    v[-2] = 0
    mat[-2][-4:] = [0, 0, 2, 6 * x[-1]] # f''(xn) = 0

    params = np.linalg.solve(mat, v) # solve linear system

    formuls = []
    for i in range(spline_num):
        abcd = params[i*4:i*4+4] # get formul
        f = lambda x, p = abcd: p[0] + p[1] * x + p[2] * (x ** 2) + p[3] * (x ** 3) # generating function
        formuls.append(f)

    return formuls


def compute_points(x, fs, num = 5):
    x1 = []
    y1 = []
        
    for xi in range(len(x) - 1):
        for xj in np.linspace(x[xi], x[xi+1], num = num):
            x1.append(xj)
            y1.append(fs[xi](xj))
    return x1, y1


if __name__ == "__main__":
    x = [-1, 1, 2, 4] 
    y = [ 4, 9, 1, 6]

    # random data
    x = []
    y = []
    t = -10
    for i in range(20):
        t = t + random.randint(1, 5)
        x.append(t)
        y.append(random.random())

    splines1 = line_spline(x, y)
    splines3 = cube_spline(x, y)
    x1, y1 = compute_points(x, splines1, 5) # (x - x'es, splines - functions list, 5 - num of dots for computing)
    x3, y3 = compute_points(x, splines3, 30)

    plt.plot(x1, y1, "-", label = "Линейный сплайн")
    plt.plot(x3, y3, "-", label = "Кубический сплайн")
    plt.plot(x, y, ".", label = "Исходные данные", mew = 4)
    plt.legend()
    plt.tight_layout()
    plt.savefig("splines.png", dpi = 300)
    plt.show()

