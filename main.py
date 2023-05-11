import numpy as np
import sympy as sm
import scipy as sc
import matplotlib.pyplot as plt

# Мантисса
k = 100

# модуль Юнга (упругости) и коэффициент Пуассона для стали
E = sm.N(2 * 10 ** 11, k)
v = sm.N(0.3, k)

# Размеры пластины (в мм)
h = sm.N(0.01, k)
a = sm.N(1, k)
b = sm.N(0.5, k)
a1 = a / 2
a2 = a / 2
b1 = b
b2 = b

# Цилиндрическая жесткость пластинки
D = (E * h ** 3) / (12 * (1 - v ** 2))

# Распределенная равномерная нагрузка q = const
q = 1

x = sm.var('x')
y = sm.var('y')
beta = sm.var('beta')
alfa = sm.var('alfa')

M = int(input("Введите размерность: "))
print('\n')

wxn = sm.zeros(M, 1)
txn = sm.zeros(M, 1)
Mxn = sm.zeros(M, 1)
Vxn = sm.zeros(M, 1)

wxna = sm.zeros(M, 1)
txna = sm.zeros(M, 1)
Mxna = sm.zeros(M, 1)
Vxna = sm.zeros(M, 1)

wym = sm.zeros(M, 1)
tym = sm.zeros(M, 1)
Mym = sm.zeros(M, 1)
Vym = sm.zeros(M, 1)

wyma = sm.zeros(M, 1)
tyma = sm.zeros(M, 1)
Myma = sm.zeros(M, 1)
Vyma = sm.zeros(M, 1)


for i in range(M):
    wxn[i] = sm.var('wxn' + str(i))
    txn[i] = sm.var('txn' + str(i))
    Mxn[i] = sm.var('Mxn' + str(i))
    Vxn[i] = sm.var('Vxn' + str(i))
    wym[i] = sm.var('wym' + str(i))
    tym[i] = sm.var('tym' + str(i))
    Mym[i] = sm.var('Mym' + str(i))
    Vym[i] = sm.var('Vym' + str(i))

    wxna[i] = sm.var('wxna' + str(i))
    txna[i] = sm.var('txna' + str(i))
    Mxna[i] = sm.var('Mxna' + str(i))
    Vxna[i] = sm.var('Vxna' + str(i))
    wyma[i] = sm.var('wyma' + str(i))
    tyma[i] = sm.var('tyma' + str(i))
    Myma[i] = sm.var('Myma' + str(i))
    Vyma[i] = sm.var('Vyma' + str(i))


X1 = sm.cosh(beta * x) + 1 / 2 * (sm.sinh(beta * x) * x * beta * v) - 1 / 2 * (sm.sinh(beta * x) * x * beta)
X2 = -(beta * sm.cosh(beta * x) * x * v - beta * sm.cosh(beta * x) * x - sm.sinh(beta * x) * v - sm.sinh(beta * x)) / (2 * beta)
X3 = -(sm.sinh(beta * x) * x) / (2 * D * beta)
X4 = -(beta * sm.cosh(beta * x) * x - sm.sinh(beta * x)) / (2 * D * beta ** 3)

Y1 = sm.cosh(alfa * y) + 1 / 2 * sm.sinh(alfa * y) * y * alfa * v - 1 / 2 * sm.sinh(alfa * y) * y * alfa
Y2 = -(alfa * sm.cosh(alfa * y) * y * v - alfa * sm.cosh(alfa * y) * y - sm.sinh(alfa * y) * v - sm.sinh(alfa * y)) / (2 * alfa)
Y3 = -(sm.sinh(alfa * y) * y) / (2 * D * alfa)
Y4 = -(alfa * sm.cosh(alfa * y) * y - sm.sinh(alfa * y)) / (2 * D * alfa ** 3)


wx = sm.zeros(M, 1)
wxa = sm.zeros(M, 1)

wy = sm.zeros(M, 1)
wya = sm.zeros(M, 1)


for i in range(M):

    wx[i] = sm.var('wx' + str(i))
    wxa[i] = sm.var('wxa' + str(i))

    wy[i] = sm.var('wy' + str(i))
    wya[i] = sm.var('wya' + str(i))


sumwx = 0
sumwxa = 0

sumwy = 0
sumwya = 0

sumwpar = 0
sumwpara = 0


for n in range(0, M):

    wx[n] = (X1 * wxn[n] + X2 * txn[n] + X3 * Mxn[n] + X4 * Vxn[n]) * sm.sin(beta * y)
    wx[n] = wx[n].subs(beta, (((n+1) * sm.pi.evalf(k)) / b1))
    sumwx += wx[n]

    wxa[n] = (X1 * wxna[n] + X2 * txna[n] + X3 * Mxna[n] + X4 * Vxna[n]) * sm.sin(beta * y)
    wxa[n] = wxa[n].subs(beta, (((n + 1) * sm.pi.evalf(k)) / b2))
    sumwxa += wxa[n]


for m in range(0, M):

    wy[m] = (Y1 * wym[m] + Y2 * tym[m] + Y3 * Mym[m] + Y4 * Vym[m]) * sm.sin(alfa * x)
    wy[m] = wy[m].subs(alfa, (((m+1) * sm.pi.evalf(k)) / a1))
    sumwy += wy[m]

    wya[m] = (Y1 * wyma[m] + Y2 * tyma[m] + Y3 * Myma[m] + Y4 * Vyma[m]) * sm.sin(alfa * x)
    wya[m] = wya[m].subs(alfa, (((m+1) * sm.pi.evalf(k)) / a2))
    sumwya += wya[m]


for n in range(1, M + 1):
    for m in range(1, M + 1):

        wpar = (4 * ((-1) ** (m + n) - (-1) ** m - (-1) ** n + 1) * sm.sin(m * x * sm.pi.evalf(k) / a1) * sm.sin(
            n * y * sm.pi.evalf(k) / b1)) / (sm.pi.evalf(k) ** 6 * m * n * D * ((m ** 2) / (a1 ** 2) + (n ** 2) / (b1 ** 2)) ** 2)
        sumwpar += wpar

        wpara = (4 * ((-1) ** (m + n) - (-1) ** m - (-1) ** n + 1) * sm.sin(m * x * sm.pi.evalf(k) / a2) * sm.sin(
            n * y * sm.pi.evalf(k) / b2)) / (sm.pi.evalf(k) ** 6 * m * n * D * ((m ** 2) / (a2 ** 2) + (n ** 2) / (b2 ** 2)) ** 2)
        sumwpara += wpara


w = sumwx + sumwy + sumwpar

wa = sumwxa + sumwya + sumwpara


tx = sm.diff(w, x, 1)
ty = sm.diff(w, y, 1)
txa = sm.diff(wa, x, 1)
tya = sm.diff(wa, y, 1)

Mx = -D * (sm.diff(w, x, 2) + v * sm.diff(w, y, 2))
My = -D * (sm.diff(w, y, 2) + v * sm.diff(w, x, 2))
Mxa = -D * (sm.diff(wa, x, 2) + v * sm.diff(wa, y, 2))
Mya = -D * (sm.diff(wa, y, 2) + v * sm.diff(wa, x, 2))

Vx = -D * (sm.diff(w, x, 3) + (2 - v) * sm.diff(w, y, 2, x, 1))
Vy = -D * (sm.diff(w, y, 3) + (2 - v) * sm.diff(w, x, 2, y, 1))
Vxa = -D * (sm.diff(wa, x, 3) + (2 - v) * sm.diff(wa, y, 2, x, 1))
Vya = -D * (sm.diff(wa, y, 3) + (2 - v) * sm.diff(wa, x, 2, y, 1))


ww0 = sm.zeros(M, 1)
ww1 = sm.zeros(M, 1)
ww2 = sm.zeros(M, 1)
ww3 = sm.zeros(M, 1)

txtx0 = sm.zeros(M, 1)
txtx1 = sm.zeros(M, 1)
tyty0 = sm.zeros(M, 1)
tyty1 = sm.zeros(M, 1)

ww0a = sm.zeros(M, 1)
ww2a = sm.zeros(M, 1)
ww3a = sm.zeros(M, 1)

txtx0a = sm.zeros(M, 1)
tyty0a = sm.zeros(M, 1)
tyty1a = sm.zeros(M, 1)

p1 = sm.zeros(M, 1)
p2 = sm.zeros(M, 1)
p3 = sm.zeros(M, 1)
p4 = sm.zeros(M, 1)


for i in range(M):

    #ww0[i] = sm.var('ww0' + str(i))
    ## ww1[i] = sm.var('ww1' + str(i))
    #ww2[i] = sm.var('ww2' + str(i))
    #ww3[i] = sm.var('ww3' + str(i))

    #ww0a[i] = sm.var('ww0a' + str(i))
    #ww2a[i] = sm.var('ww2a' + str(i))
    #ww3a[i] = sm.var('ww3a' + str(i))

    #txtx0[i] = sm.var('txtx0' + str(i))
    ## txtx1[i] = sm.var('txtx1' + str(i))
    #tyty0[i] = sm.var('tyty0' + str(i))
    #tyty1[i] = sm.var('tyty1' + str(i))

    #txtx0a[i] = sm.var('txtx0a' + str(i))
    #tyty0a[i] = sm.var('tyty0a' + str(i))
    #tyty1a[i] = sm.var('tyty1a' + str(i))

    a1i = (a1 / (M + 1)) * (i + 1)
    a2i = (a2 / (M + 1)) * (i + 1)
    b1i = (b1 / (M + 1)) * (i + 1)
    b2i = (b2 / (M + 1)) * (i + 1)

    ww0[i] = w.subs([(x, 0), (y, b1i)])
    ww1[i] = w.subs([(x, a1), (y, b1i)])
    ww2[i] = w.subs([(x, a1i), (y, 0)])
    ww3[i] = w.subs([(x, a1i), (y, b1)])

    txtx0[i] = tx.subs([(x, 0), (y, b1i)])
    txtx1[i] = tx.subs([(x, a1), (y, b1i)])
    tyty0[i] = ty.subs([(x, a1i), (y, 0)])
    tyty1[i] = ty.subs([(x, a1i), (y, b1)])

    ww0a[i] = wa.subs([(x, a2), (y, b2i)])
    ww2a[i] = wa.subs([(x, a2i), (y, 0)])
    ww3a[i] = wa.subs([(x, a2i), (y, b2)])

    txtx0a[i] = txa.subs([(x, a2), (y, b2i)])
    tyty0a[i] = tya.subs([(x, a2i), (y, 0)])
    tyty1a[i] = tya.subs([(x, a2i), (y, b2)])

    wsubs = w.subs([(x, a1), (y, b1i)])
    txsubs = tx.subs([(x, a1), (y, b1i)])
    Mxsubs = Mx.subs([(x, a1), (y, b1i)])
    Vxsubs = Vx.subs([(x, a1), (y, b1i)])
    wasubs = wa.subs([(x, 0), (y, b2i)])
    txasubs = txa.subs([(x, 0), (y, b2i)])
    Mxasubs = Mxa.subs([(x, 0), (y, b2i)])
    Vxasubs = Vxa.subs([(x, 0), (y, b2i)])

    p1[i] = wsubs - wasubs
    p2[i] = txsubs - txasubs
    p3[i] = Mxsubs - Mxasubs
    p4[i] = Vxsubs - Vxasubs

system = sm.zeros(16 * M, 1)
symbols = sm.zeros(16 * M, 1)

"""
for j in range(M):

    system[j] = ww0[j]
    symbols[j] = wxn[j]
    system[j + M] = ww1[j]
    symbols[j + M] = txn[j]
    system[j + 2 * M] = ww2[j]
    symbols[j + 2 * M] = Mxn[j]
    system[j + 3 * M] = ww3[j]
    symbols[j + 3 * M] = Vxn[j]
    system[j + 4 * M] = txtx0[j]
    symbols[j + 4 * M] = wym[j]
    system[j + 5 * M] = txtx1[j]
    symbols[j + 5 * M] = tym[j]
    system[j + 6 * M] = tyty0[j]
    symbols[j + 6 * M] = Mym[j]
    system[j + 7 * M] = tyty1[j]
    symbols[j + 7 * M] = Vym[j]
"""

for j in range(M):

    system[j] = ww0[j]
    system[j + M] = ww2[j]
    system[j + 2 * M] = ww3[j]
    system[j + 3 * M] = txtx0[j]
    system[j + 4 * M] = tyty0[j]
    system[j + 5 * M] = tyty1[j]
    system[j + 6 * M] = ww0a[j]
    system[j + 7 * M] = ww2a[j]
    system[j + 8 * M] = ww3a[j]
    system[j + 9 * M] = txtx0a[j]
    system[j + 10 * M] = tyty0a[j]
    system[j + 11 * M] = tyty1a[j]
    system[j + 12 * M] = p1[j]
    system[j + 13 * M] = p2[j]
    system[j + 14 * M] = p3[j]
    system[j + 15 * M] = p4[j]

    symbols[j] = wxn[j]
    symbols[j + M] = txn[j]
    symbols[j + 2 * M] = Mxn[j]
    symbols[j + 3 * M] = Vxn[j]
    symbols[j + 4 * M] = wym[j]
    symbols[j + 5 * M] = tym[j]
    symbols[j + 6 * M] = Mym[j]
    symbols[j + 7 * M] = Vym[j]
    symbols[j + 8 * M] = wxna[j]
    symbols[j + 9 * M] = txna[j]
    symbols[j + 10 * M] = Mxna[j]
    symbols[j + 11 * M] = Vxna[j]
    symbols[j + 12 * M] = wyma[j]
    symbols[j + 13 * M] = tyma[j]
    symbols[j + 14 * M] = Myma[j]
    symbols[j + 15 * M] = Vyma[j]


mtx = sm.linear_eq_to_matrix(system, *symbols)

#mtx1 = sm.matrix2numpy(mtx[0], dtype=np.longdouble)
#mtx2 = sm.matrix2numpy(mtx[1], dtype=np.longdouble)

#slv = sc.linalg.solve(mtx1, mtx2)
slv = sm.linsolve((mtx[0], mtx[1]), *symbols)

print('Решение системы уравнений: \n\n')
print(np.asarray(slv.args).T, '\n\n')


for i in range(M):

    wxn = sm.var('wxn' + str(i))
    txn = sm.var('txn' + str(i))
    Mxn = sm.var('Mxn' + str(i))
    Vxn = sm.var('Vxn' + str(i))
    wym = sm.var('wym' + str(i))
    tym = sm.var('tym' + str(i))
    Mym = sm.var('Mym' + str(i))
    Vym = sm.var('Vym' + str(i))

    wxna = sm.var('wxna' + str(i))
    txna = sm.var('txna' + str(i))
    Mxna = sm.var('Mxna' + str(i))
    Vxna = sm.var('Vxna' + str(i))
    wyma = sm.var('wyma' + str(i))
    tyma = sm.var('tyma' + str(i))
    Myma = sm.var('Myma' + str(i))
    Vyma = sm.var('Vyma' + str(i))


    w1 = slv.args[0][i]
    t1 = slv.args[0][i + M]
    M1 = slv.args[0][i + 2 * M]
    V1 = slv.args[0][i + 3 * M]
    w2 = slv.args[0][i + 4 * M]
    t2 = slv.args[0][i + 5 * M]
    M2 = slv.args[0][i + 6 * M]
    V2 = slv.args[0][i + 7 * M]

    w1a = slv.args[0][i + 8 * M]
    t1a = slv.args[0][i + 9 * M]
    M1a = slv.args[0][i + 10 * M]
    V1a = slv.args[0][i + 11 * M]
    w2a = slv.args[0][i + 12 * M]
    t2a = slv.args[0][i + 13 * M]
    M2a = slv.args[0][i + 14 * M]
    V2a = slv.args[0][i + 15 * M]

    w = w.subs([
        (wxn, w1),
        (txn, t1),
        (Mxn, M1),
        (Vxn, V1),
        (wym, w2),
        (tym, t2),
        (Mym, M2),
        (Vym, V2)
    ])
    wa = wa.subs([
        (wxna, w1a),
        (txna, t1a),
        (Mxna, M1a),
        (Vxna, V1a),
        (wyma, w2a),
        (tyma, t2a),
        (Myma, M2a),
        (Vyma, V2a)
    ])
    Mx = Mx.subs([
        (wxn, w1),
        (txn, t1),
        (Mxn, M1),
        (Vxn, V1),
        (wym, w2),
        (tym, t2),
        (Mym, M2),
        (Vym, V2)
    ])
    Mxa = Mxa.subs([
        (wxna, w1a),
        (txna, t1a),
        (Mxna, M1a),
        (Vxna, V1a),
        (wyma, w2a),
        (tyma, t2a),
        (Myma, M2a),
        (Vyma, V2a)
    ])
    Vx = Vx.subs([
        (wxn, w1),
        (txn, t1),
        (Mxn, M1),
        (Vxn, V1),
        (wym, w2),
        (tym, t2),
        (Mym, M2),
        (Vym, V2)
    ])
    Vxa = Vxa.subs([
        (wxna, w1a),
        (txna, t1a),
        (Mxna, M1a),
        (Vxna, V1a),
        (wyma, w2a),
        (tyma, t2a),
        (Myma, M2a),
        (Vyma, V2a)
    ])


ww = (w * D) / (q * a1 ** 4)

wwa = (wa * D) / (q * a2 ** 4)

MMx = Mx / (q * a1 ** 2)
MMxa = Mxa / (q * a2 ** 2)

VVx = Vx / (q * a1)
VVxa = Vxa / (q * a2)


www = sm.zeros(M + 2, 1)
wwwa = sm.zeros(M + 2, 1)
MMMx = sm.zeros(M + 2, 1)
MMMxa = sm.zeros(M + 2, 1)
VVVx = sm.zeros(M + 2, 1)
VVVxa = sm.zeros(M + 2, 1)
xi = sm.zeros(M + 2, 1)
xai = sm.zeros(M + 2, 1)
yi = sm.zeros(M + 2, 1)
yai = sm.zeros(M + 2, 1)


print(ww.subs([(x, a1), (y, b1 / 2)]))
print(wwa.subs([(x, 0), (y, b2 / 2)]))


for i in range(M + 2):

    xi[i] = (a1 / (M + 1)) * i
    xai[i] = (a2 / (M + 1)) * i
    yi[i] = b1 / 2
    yai[i] = b2 / 2
    www[i] = ww.subs([(x, xi[i]), (y, yi[i])])
    wwwa[i] = wwa.subs([(x, xai[i]), (y, yai[i])])

"""
for i in range(M + 2):
    xi[i] = 0
    yi[i] = (b1 / (M + 1)) * i
    xai[i] = a2
    yai[i] = (b2 / (M + 1)) * i
    MMMx[i] = MMx.subs([(x, xi[i]), (y, yi[i])])
    MMMxa[i] = MMxa.subs([(x, xai[i]), (y, yai[i])])
"""
"""
for i in range(M + 2):
    xi[i] = 0
    yi[i] = (b1 / (M + 1)) * i
    xai[i] = a2
    yai[i] = (b2 / (M + 1)) * i
    VVVx[i] = VVx.subs([(x, xi[i]), (y, yi[i])])
    VVVxa[i] = VVxa.subs([(x, xai[i]), (y, yai[i])])
"""

plt.subplot(121)
plt.plot(xi, www)
plt.subplot(122)
plt.plot(xai, wwwa)

"""
plt.subplot(121)
plt.plot(yi, VVVx)
plt.subplot(122)
plt.plot(yai, VVVxa)
"""
plt.show()
