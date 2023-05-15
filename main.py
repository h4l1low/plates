import numpy as np
import sympy as sm
import matplotlib.pyplot as plt

# Мантисса
k = 100

# модуль Юнга (упругости) и коэффициент Пуассона для стали
E = sm.N(2 * 10 ** 11, k)
v = sm.N(0.3, k)

# Размеры пластины (в мм)
h = sm.N(0.01, k)
a = sm.N(3, k)
b = sm.N(2, k)
a1 = a / 3
a2 = a / 3
a3 = a / 3
a4 = a / 3
a5 = a / 3
b1 = b / 2
b2 = b / 2
b3 = b / 2
b4 = b / 2
b5 = b / 2

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

wxna = sm.zeros(M, 1)   # Эти записи и дальнейшие создают матрицу нулей размера M,1
txna = sm.zeros(M, 1)
Mxna = sm.zeros(M, 1)
Vxna = sm.zeros(M, 1)

wxnb = sm.zeros(M, 1)
txnb = sm.zeros(M, 1)
Mxnb = sm.zeros(M, 1)
Vxnb = sm.zeros(M, 1)

wxnc = sm.zeros(M, 1)
txnc = sm.zeros(M, 1)
Mxnc = sm.zeros(M, 1)
Vxnc = sm.zeros(M, 1)

wxnd = sm.zeros(M, 1)
txnd = sm.zeros(M, 1)
Mxnd = sm.zeros(M, 1)
Vxnd = sm.zeros(M, 1)

wxne = sm.zeros(M, 1)
txne = sm.zeros(M, 1)
Mxne = sm.zeros(M, 1)
Vxne = sm.zeros(M, 1)

wyma = sm.zeros(M, 1)
tyma = sm.zeros(M, 1)
Myma = sm.zeros(M, 1)
Vyma = sm.zeros(M, 1)

wymb = sm.zeros(M, 1)
tymb = sm.zeros(M, 1)
Mymb = sm.zeros(M, 1)
Vymb = sm.zeros(M, 1)

wymc = sm.zeros(M, 1)
tymc = sm.zeros(M, 1)
Mymc = sm.zeros(M, 1)
Vymc = sm.zeros(M, 1)

wymd = sm.zeros(M, 1)
tymd = sm.zeros(M, 1)
Mymd = sm.zeros(M, 1)
Vymd = sm.zeros(M, 1)

wyme = sm.zeros(M, 1)
tyme = sm.zeros(M, 1)
Myme = sm.zeros(M, 1)
Vyme = sm.zeros(M, 1)

for i in range(M):  # Эти циклы и дальнейшие создают в этих пустых матрицах наборы переменных, по M на каждое значение

    wxna[i] = sm.var('wxna' + str(i))   # Вид: wxna0, wxna1...
    txna[i] = sm.var('txna' + str(i))
    Mxna[i] = sm.var('Mxna' + str(i))
    Vxna[i] = sm.var('Vxna' + str(i))
    wyma[i] = sm.var('wyma' + str(i))
    tyma[i] = sm.var('tyma' + str(i))
    Myma[i] = sm.var('Myma' + str(i))
    Vyma[i] = sm.var('Vyma' + str(i))

    wxnb[i] = sm.var('wxnb' + str(i))
    txnb[i] = sm.var('txnb' + str(i))
    Mxnb[i] = sm.var('Mxnb' + str(i))
    Vxnb[i] = sm.var('Vxnb' + str(i))
    wymb[i] = sm.var('wymb' + str(i))
    tymb[i] = sm.var('tymb' + str(i))
    Mymb[i] = sm.var('Mymb' + str(i))
    Vymb[i] = sm.var('Vymb' + str(i))

    wxnc[i] = sm.var('wxnc' + str(i))
    txnc[i] = sm.var('txnc' + str(i))
    Mxnc[i] = sm.var('Mxnc' + str(i))
    Vxnc[i] = sm.var('Vxnc' + str(i))
    wymc[i] = sm.var('wymc' + str(i))
    tymc[i] = sm.var('tymc' + str(i))
    Mymc[i] = sm.var('Mymc' + str(i))
    Vymc[i] = sm.var('Vymc' + str(i))

    wxnd[i] = sm.var('wxnd' + str(i))
    txnd[i] = sm.var('txnd' + str(i))
    Mxnd[i] = sm.var('Mxnd' + str(i))
    Vxnd[i] = sm.var('Vxnd' + str(i))
    wymd[i] = sm.var('wymd' + str(i))
    tymd[i] = sm.var('tymd' + str(i))
    Mymd[i] = sm.var('Mymd' + str(i))
    Vymd[i] = sm.var('Vymd' + str(i))

    wxne[i] = sm.var('wxne' + str(i))
    txne[i] = sm.var('txne' + str(i))
    Mxne[i] = sm.var('Mxne' + str(i))
    Vxne[i] = sm.var('Vxne' + str(i))
    wyme[i] = sm.var('wyme' + str(i))
    tyme[i] = sm.var('tyme' + str(i))
    Myme[i] = sm.var('Myme' + str(i))
    Vyme[i] = sm.var('Vyme' + str(i))

print('Считаем X,Y...')

X1 = sm.cosh(beta * x) + 1 / 2 * (sm.sinh(beta * x) * x * beta * v) - 1 / 2 * (sm.sinh(beta * x) * x * beta)
X2 = -(beta * sm.cosh(beta * x) * x * v - beta * sm.cosh(beta * x) * x - sm.sinh(beta * x) * v - sm.sinh(beta * x)) / (
            2 * beta)
X3 = -(sm.sinh(beta * x) * x) / (2 * D * beta)
X4 = -(beta * sm.cosh(beta * x) * x - sm.sinh(beta * x)) / (2 * D * beta ** 3)

Y1 = sm.cosh(alfa * y) + 1 / 2 * sm.sinh(alfa * y) * y * alfa * v - 1 / 2 * sm.sinh(alfa * y) * y * alfa
Y2 = -(alfa * sm.cosh(alfa * y) * y * v - alfa * sm.cosh(alfa * y) * y - sm.sinh(alfa * y) * v - sm.sinh(alfa * y)) / (
            2 * alfa)
Y3 = -(sm.sinh(alfa * y) * y) / (2 * D * alfa)
Y4 = -(alfa * sm.cosh(alfa * y) * y - sm.sinh(alfa * y)) / (2 * D * alfa ** 3)

print('Готово\n')

wxa = sm.zeros(M, 1)
wxb = sm.zeros(M, 1)
wxc = sm.zeros(M, 1)
wxd = sm.zeros(M, 1)
wxe = sm.zeros(M, 1)

wya = sm.zeros(M, 1)
wyb = sm.zeros(M, 1)
wyc = sm.zeros(M, 1)
wyd = sm.zeros(M, 1)
wye = sm.zeros(M, 1)

for i in range(M):

    wxa[i] = sm.var('wxa' + str(i))
    wxb[i] = sm.var('wxb' + str(i))
    wxc[i] = sm.var('wxc' + str(i))
    wxd[i] = sm.var('wxd' + str(i))
    wxe[i] = sm.var('wxe' + str(i))

    wya[i] = sm.var('wya' + str(i))
    wyb[i] = sm.var('wyb' + str(i))
    wyc[i] = sm.var('wyc' + str(i))
    wyd[i] = sm.var('wyd' + str(i))
    wye[i] = sm.var('wye' + str(i))

sumwxa = 0
sumwxb = 0
sumwxc = 0
sumwxd = 0
sumwxe = 0

sumwya = 0
sumwyb = 0
sumwyc = 0
sumwyd = 0
sumwye = 0

sumwpara = 0
sumwparb = 0
sumwparc = 0
sumwpard = 0
sumwpare = 0

print('Считаем w...')

for n in range(0, M):   # Суммы по x
    wxa[n] = (X1 * wxna[n] + X2 * txna[n] + X3 * Mxna[n] + X4 * Vxna[n]) * sm.sin(beta * y)
    wxa[n] = wxa[n].subs(beta, (((n + 1) * sm.pi.evalf(k)) / b1))
    sumwxa += wxa[n]

    wxb[n] = (X1 * wxnb[n] + X2 * txnb[n] + X3 * Mxnb[n] + X4 * Vxnb[n]) * sm.sin(beta * y)
    wxb[n] = wxb[n].subs(beta, (((n + 1) * sm.pi.evalf(k)) / b2))
    sumwxb += wxb[n]

    wxc[n] = (X1 * wxnc[n] + X2 * txnc[n] + X3 * Mxnc[n] + X4 * Vxnc[n]) * sm.sin(beta * y)
    wxc[n] = wxc[n].subs(beta, (((n + 1) * sm.pi.evalf(k)) / b3))
    sumwxc += wxc[n]

    wxd[n] = (X1 * wxnd[n] + X2 * txnd[n] + X3 * Mxnd[n] + X4 * Vxnd[n]) * sm.sin(beta * y)
    wxd[n] = wxd[n].subs(beta, (((n + 1) * sm.pi.evalf(k)) / b4))
    sumwxd += wxd[n]

    wxe[n] = (X1 * wxne[n] + X2 * txne[n] + X3 * Mxne[n] + X4 * Vxne[n]) * sm.sin(beta * y)
    wxe[n] = wxe[n].subs(beta, (((n + 1) * sm.pi.evalf(k)) / b5))
    sumwxe += wxe[n]

for m in range(0, M):   # Суммы по y
    wya[m] = (Y1 * wyma[m] + Y2 * tyma[m] + Y3 * Myma[m] + Y4 * Vyma[m]) * sm.sin(alfa * x)
    wya[m] = wya[m].subs(alfa, (((m + 1) * sm.pi.evalf(k)) / a1))
    sumwya += wya[m]

    wyb[m] = (Y1 * wymb[m] + Y2 * tymb[m] + Y3 * Mymb[m] + Y4 * Vymb[m]) * sm.sin(alfa * x)
    wyb[m] = wyb[m].subs(alfa, (((m + 1) * sm.pi.evalf(k)) / a2))
    sumwyb += wyb[m]

    wyc[m] = (Y1 * wymc[m] + Y2 * tymc[m] + Y3 * Mymc[m] + Y4 * Vymc[m]) * sm.sin(alfa * x)
    wyc[m] = wyc[m].subs(alfa, (((m + 1) * sm.pi.evalf(k)) / a3))
    sumwyc += wyc[m]

    wyd[m] = (Y1 * wymd[m] + Y2 * tymd[m] + Y3 * Mymd[m] + Y4 * Vymd[m]) * sm.sin(alfa * x)
    wyd[m] = wyd[m].subs(alfa, (((m + 1) * sm.pi.evalf(k)) / a4))
    sumwyd += wyd[m]

    wye[m] = (Y1 * wyme[m] + Y2 * tyme[m] + Y3 * Myme[m] + Y4 * Vyme[m]) * sm.sin(alfa * x)
    wye[m] = wye[m].subs(alfa, (((m + 1) * sm.pi.evalf(k)) / a5))
    sumwye += wye[m]

for n in range(1, M + 1):   # Частные значения
    for m in range(1, M + 1):
        wpara = (4 * ((-1) ** (m + n) - (-1) ** m - (-1) ** n + 1) * sm.sin(m * x * sm.pi.evalf(k) / a1) * sm.sin(
            n * y * sm.pi.evalf(k) / b1)) / (
                            sm.pi.evalf(k) ** 6 * m * n * D * ((m ** 2) / (a1 ** 2) + (n ** 2) / (b1 ** 2)) ** 2)
        sumwpara += wpara

        wparb = (4 * ((-1) ** (m + n) - (-1) ** m - (-1) ** n + 1) * sm.sin(m * x * sm.pi.evalf(k) / a2) * sm.sin(
            n * y * sm.pi.evalf(k) / b2)) / (
                            sm.pi.evalf(k) ** 6 * m * n * D * ((m ** 2) / (a2 ** 2) + (n ** 2) / (b2 ** 2)) ** 2)
        sumwparb += wparb

        wparc = (4 * ((-1) ** (m + n) - (-1) ** m - (-1) ** n + 1) * sm.sin(m * x * sm.pi.evalf(k) / a3) * sm.sin(
            n * y * sm.pi.evalf(k) / b3)) / (
                        sm.pi.evalf(k) ** 6 * m * n * D * ((m ** 2) / (a3 ** 2) + (n ** 2) / (b3 ** 2)) ** 2)
        sumwparc += wparc

        wpard = (4 * ((-1) ** (m + n) - (-1) ** m - (-1) ** n + 1) * sm.sin(m * x * sm.pi.evalf(k) / a4) * sm.sin(
            n * y * sm.pi.evalf(k) / b4)) / (
                        sm.pi.evalf(k) ** 6 * m * n * D * ((m ** 2) / (a4 ** 2) + (n ** 2) / (b4 ** 2)) ** 2)
        sumwpard += wpard

        wpare = (4 * ((-1) ** (m + n) - (-1) ** m - (-1) ** n + 1) * sm.sin(m * x * sm.pi.evalf(k) / a5) * sm.sin(
            n * y * sm.pi.evalf(k) / b5)) / (
                        sm.pi.evalf(k) ** 6 * m * n * D * ((m ** 2) / (a5 ** 2) + (n ** 2) / (b5 ** 2)) ** 2)
        sumwpare += wpare

wa = sumwxa + sumwya + sumwpara

wb = sumwxb + sumwyb + sumwparb

wc = sumwxc + sumwyc + sumwparc

wd = sumwxd + sumwyd + sumwpard

we = sumwxe + sumwye + sumwpare

print('Готово\n')

print('Считаем дифференциалы...')

txa = sm.diff(wa, x, 1)
tya = sm.diff(wa, y, 1)
txb = sm.diff(wb, x, 1)
tyb = sm.diff(wb, y, 1)
txc = sm.diff(wc, x, 1)
tyc = sm.diff(wc, y, 1)
txd = sm.diff(wd, x, 1)
tyd = sm.diff(wd, y, 1)
txe = sm.diff(we, x, 1)
tye = sm.diff(we, y, 1)

Mxa = -D * (sm.diff(wa, x, 2) + v * sm.diff(wa, y, 2))
Mya = -D * (sm.diff(wa, y, 2) + v * sm.diff(wa, x, 2))
Mxb = -D * (sm.diff(wb, x, 2) + v * sm.diff(wb, y, 2))
Myb = -D * (sm.diff(wb, y, 2) + v * sm.diff(wb, x, 2))
Mxc = -D * (sm.diff(wc, x, 2) + v * sm.diff(wc, y, 2))
Myc = -D * (sm.diff(wc, y, 2) + v * sm.diff(wc, x, 2))
Mxd = -D * (sm.diff(wd, x, 2) + v * sm.diff(wd, y, 2))
Myd = -D * (sm.diff(wd, y, 2) + v * sm.diff(wd, x, 2))
Mxe = -D * (sm.diff(we, x, 2) + v * sm.diff(we, y, 2))
Mye = -D * (sm.diff(we, y, 2) + v * sm.diff(we, x, 2))

Vxa = -D * (sm.diff(wa, x, 3) + (2 - v) * sm.diff(wa, y, 2, x, 1))
Vya = -D * (sm.diff(wa, y, 3) + (2 - v) * sm.diff(wa, x, 2, y, 1))
Vxb = -D * (sm.diff(wb, x, 3) + (2 - v) * sm.diff(wb, y, 2, x, 1))
Vyb = -D * (sm.diff(wb, y, 3) + (2 - v) * sm.diff(wb, x, 2, y, 1))
Vxc = -D * (sm.diff(wc, x, 3) + (2 - v) * sm.diff(wc, y, 2, x, 1))
Vyc = -D * (sm.diff(wc, y, 3) + (2 - v) * sm.diff(wc, x, 2, y, 1))
Vxd = -D * (sm.diff(wd, x, 3) + (2 - v) * sm.diff(wd, y, 2, x, 1))
Vyd = -D * (sm.diff(wd, y, 3) + (2 - v) * sm.diff(wd, x, 2, y, 1))
Vxe = -D * (sm.diff(we, x, 3) + (2 - v) * sm.diff(we, y, 2, x, 1))
Vye = -D * (sm.diff(we, y, 3) + (2 - v) * sm.diff(we, x, 2, y, 1))

print('Готово\n')

wa1 = sm.zeros(M, 1)
wa2 = sm.zeros(M, 1)
txa1 = sm.zeros(M, 1)
tya1 = sm.zeros(M, 1)

wb1 = sm.zeros(M, 1)
tyb1 = sm.zeros(M, 1)
Myb1 = sm.zeros(M, 1)
Vyb1 = sm.zeros(M, 1)

wc1 = sm.zeros(M, 1)
wc2 = sm.zeros(M, 1)
txc1 = sm.zeros(M, 1)
tyc1 = sm.zeros(M, 1)
Mxc1 = sm.zeros(M, 1)
Vxc1 = sm.zeros(M, 1)

wd1 = sm.zeros(M, 1)
wd2 = sm.zeros(M, 1)
tyd1 = sm.zeros(M, 1)
txd1 = sm.zeros(M, 1)

we1 = sm.zeros(M, 1)
we2 = sm.zeros(M, 1)
txe1 = sm.zeros(M, 1)
tye1 = sm.zeros(M, 1)
Mxe1 = sm.zeros(M, 1)
Vxe1 = sm.zeros(M, 1)

wab = sm.zeros(M, 1)
txab = sm.zeros(M, 1)
Mxab = sm.zeros(M, 1)
Vxab = sm.zeros(M, 1)

wac = sm.zeros(M, 1)
tyac = sm.zeros(M, 1)
Myac = sm.zeros(M, 1)
Vyac = sm.zeros(M, 1)

wbd = sm.zeros(M, 1)
txbd = sm.zeros(M, 1)
Mxbd = sm.zeros(M, 1)
Vxbd = sm.zeros(M, 1)

wde = sm.zeros(M, 1)
tyde = sm.zeros(M, 1)
Myde = sm.zeros(M, 1)
Vyde = sm.zeros(M, 1)

print('Составляем систему из граничных условий...')

for i in range(M):
    a1i = (a1 / (M + 1)) * (i + 1)
    a2i = (a2 / (M + 1)) * (i + 1)
    a3i = (a3 / (M + 1)) * (i + 1)
    a4i = (a4 / (M + 1)) * (i + 1)
    a5i = (a5 / (M + 1)) * (i + 1)
    b1i = (b1 / (M + 1)) * (i + 1)
    b2i = (b2 / (M + 1)) * (i + 1)
    b3i = (b3 / (M + 1)) * (i + 1)
    b4i = (b4 / (M + 1)) * (i + 1)
    b5i = (b5 / (M + 1)) * (i + 1)

    # 1
    wa1[i] = wa.subs([(x, 0), (y, b1i)])
    wa2[i] = wa.subs([(x, a1i), (y, 0)])
    txa1[i] = txa.subs([(x, 0), (y, b1i)])
    tya1[i] = tya.subs([(x, a1i), (y, 0)])

    # 2
    wb1[i] = wb.subs([(x, a2i), (y, 0)])
    tyb1[i] = tyb.subs([(x, a2i), (y, 0)])
    Myb1[i] = Myb.subs([(x, a2i), (y, b2)])
    Vyb1[i] = Vyb.subs([(x, a2i), (y, b2)])

    # 3
    wc1[i] = wc.subs([(x, 0), (y, b3i)])
    wc2[i] = wc.subs([(x, a3i), (y, b3)])
    txc1[i] = txc.subs([(x, 0), (y, b3i)])
    tyc1[i] = tyc.subs([(x, a3i), (y, b3)])
    Mxc1[i] = Mxc.subs([(x, a3), (y, b3i)])
    Vxc1[i] = Vxc.subs([(x, a3), (y, b3i)])

    # 4
    wd1[i] = wd.subs([(x, a4i), (y, 0)])
    wd2[i] = wd.subs([(x, a4), (y, b4i)])
    tyd1[i] = tyd.subs([(x, a4i), (y, 0)])
    txd1[i] = txd.subs([(x, a4), (y, b4i)])

    # 5
    we1[i] = we.subs([(x, a5), (y, b5i)])
    we2[i] = we.subs([(x, a5i), (y, b5)])
    txe1[i] = txe.subs([(x, a5), (y, b5i)])
    tye1[i] = tye.subs([(x, a5i), (y, b5)])
    Mxe1[i] = Mxe.subs([(x, 0), (y, b5i)])
    Vxe1[i] = Vxe.subs([(x, 0), (y, b5i)])

    # 1-2
    wab[i] = wa.subs([(x, a1), (y, b1i)]) - wb.subs([(x, 0), (y, b2i)])
    txab[i] = txa.subs([(x, a1), (y, b1i)]) - txb.subs([(x, 0), (y, b2i)])
    Mxab[i] = Mxa.subs([(x, a1), (y, b1i)]) - Mxb.subs([(x, 0), (y, b2i)])
    Vxab[i] = Vxa.subs([(x, a1), (y, b1i)]) - Vxb.subs([(x, 0), (y, b2i)])

    # 1-3
    wac[i] = wa.subs([(x, a1i), (y, b1)]) - wc.subs([(x, a2i), (y, 0)])
    tyac[i] = tya.subs([(x, a1i), (y, b1)]) - tyc.subs([(x, a2i), (y, 0)])
    Myac[i] = Mya.subs([(x, a1i), (y, b1)]) - Myc.subs([(x, a2i), (y, 0)])
    Vyac[i] = Vya.subs([(x, a1i), (y, b1)]) - Vyc.subs([(x, a2i), (y, 0)])

    # 2-4
    wbd[i] = wb.subs([(x, a2), (y, b2i)]) - wd.subs([(x, 0), (y, b4i)])
    txbd[i] = txb.subs([(x, a2), (y, b2i)]) - txd.subs([(x, 0), (y, b4i)])
    Mxbd[i] = Mxb.subs([(x, a2), (y, b2i)]) - Mxd.subs([(x, 0), (y, b4i)])
    Vxbd[i] = Vxb.subs([(x, a2), (y, b2i)]) - Vxd.subs([(x, 0), (y, b4i)])

    # 4-5
    wde[i] = wd.subs([(x, a4i), (y, b4)]) - we.subs([(x, a5i), (y, 0)])
    tyde[i] = tyd.subs([(x, a4i), (y, b4)]) - tye.subs([(x, a5i), (y, 0)])
    Myde[i] = Myd.subs([(x, a4i), (y, b4)]) - Mye.subs([(x, a5i), (y, 0)])
    Vyde[i] = Vyd.subs([(x, a4i), (y, b4)]) - Vye.subs([(x, a5i), (y, 0)])

print('Готово\n')

system = sm.zeros(40 * M, 1)
symbols = sm.zeros(40 * M, 1)

for i in range(M):
    system[i] = wa1[i]
    system[i + M] = wa2[i]
    system[i + 2 * M] = txa1[i]
    system[i + 3 * M] = tya1[i]
    system[i + 4 * M] = wb1[i]
    system[i + 5 * M] = tyb1[i]
    system[i + 6 * M] = Myb1[i]
    system[i + 7 * M] = Vyb1[i]
    system[i + 8 * M] = wc1[i]
    system[i + 9 * M] = wc2[i]
    system[i + 10 * M] = txc1[i]
    system[i + 11 * M] = tyc1[i]
    system[i + 12 * M] = Mxc1[i]
    system[i + 13 * M] = Vxc1[i]
    system[i + 14 * M] = wd1[i]
    system[i + 15 * M] = wd2[i]
    system[i + 16 * M] = tyd1[i]
    system[i + 17 * M] = txd1[i]
    system[i + 18 * M] = we1[i]
    system[i + 19 * M] = we2[i]
    system[i + 20 * M] = txe1[i]
    system[i + 21 * M] = tye1[i]
    system[i + 22 * M] = Mxe1[i]
    system[i + 23 * M] = Vxe1[i]
    system[i + 24 * M] = wab[i]
    system[i + 25 * M] = txab[i]
    system[i + 26 * M] = Mxab[i]
    system[i + 27 * M] = Vxab[i]
    system[i + 28 * M] = wac[i]
    system[i + 29 * M] = tyac[i]
    system[i + 30 * M] = Myac[i]
    system[i + 31 * M] = Vyac[i]
    system[i + 32 * M] = wbd[i]
    system[i + 33 * M] = txbd[i]
    system[i + 34 * M] = Mxbd[i]
    system[i + 35 * M] = Vxbd[i]
    system[i + 36 * M] = wde[i]
    system[i + 37 * M] = tyde[i]
    system[i + 38 * M] = Myde[i]
    system[i + 39 * M] = Vyde[i]

    symbols[i] = wxna[i]
    symbols[i + M] = txna[i]
    symbols[i + 2 * M] = Mxna[i]
    symbols[i + 3 * M] = Vxna[i]
    symbols[i + 4 * M] = wyma[i]
    symbols[i + 5 * M] = tyma[i]
    symbols[i + 6 * M] = Myma[i]
    symbols[i + 7 * M] = Vyma[i]
    symbols[i + 8 * M] = wxnb[i]
    symbols[i + 9 * M] = txnb[i]
    symbols[i + 10 * M] = Mxnb[i]
    symbols[i + 11 * M] = Vxnb[i]
    symbols[i + 12 * M] = wymb[i]
    symbols[i + 13 * M] = tymb[i]
    symbols[i + 14 * M] = Mymb[i]
    symbols[i + 15 * M] = Vymb[i]
    symbols[i + 16 * M] = wxnc[i]
    symbols[i + 17 * M] = txnc[i]
    symbols[i + 18 * M] = Mxnc[i]
    symbols[i + 19 * M] = Vxnc[i]
    symbols[i + 20 * M] = wymc[i]
    symbols[i + 21 * M] = tymc[i]
    symbols[i + 22 * M] = Mymc[i]
    symbols[i + 23 * M] = Vymc[i]
    symbols[i + 24 * M] = wxnd[i]
    symbols[i + 25 * M] = txnd[i]
    symbols[i + 26 * M] = Mxnd[i]
    symbols[i + 27 * M] = Vxnd[i]
    symbols[i + 28 * M] = wymd[i]
    symbols[i + 29 * M] = tymd[i]
    symbols[i + 30 * M] = Mymd[i]
    symbols[i + 31 * M] = Vymd[i]
    symbols[i + 32 * M] = wxne[i]
    symbols[i + 33 * M] = txne[i]
    symbols[i + 34 * M] = Mxne[i]
    symbols[i + 35 * M] = Vxne[i]
    symbols[i + 36 * M] = wyme[i]
    symbols[i + 37 * M] = tyme[i]
    symbols[i + 38 * M] = Myme[i]
    symbols[i + 39 * M] = Vyme[i]

mtx = sm.linear_eq_to_matrix(system, *symbols)

# mtx1 = sm.matrix2numpy(mtx[0], dtype=np.longdouble)
# mtx2 = sm.matrix2numpy(mtx[1], dtype=np.longdouble)

# slv = sc.linalg.solve(mtx1, mtx2)

print('Решаем систему...')

slv = sm.linsolve((mtx[0], mtx[1]), *symbols)

print('Готово\n')

print('Решение системы уравнений: \n\n')
print(np.asarray(slv.args).T, '\n\n')

print('Подставляем решение системы...')

for i in range(M):
    wxna = sm.var('wxna' + str(i))
    txna = sm.var('txna' + str(i))
    Mxna = sm.var('Mxna' + str(i))
    Vxna = sm.var('Vxna' + str(i))
    wyma = sm.var('wyma' + str(i))
    tyma = sm.var('tyma' + str(i))
    Myma = sm.var('Myma' + str(i))
    Vyma = sm.var('Vyma' + str(i))

    wxnb = sm.var('wxnb' + str(i))
    txnb = sm.var('txnb' + str(i))
    Mxnb = sm.var('Mxnb' + str(i))
    Vxnb = sm.var('Vxnb' + str(i))
    wymb = sm.var('wymb' + str(i))
    tymb = sm.var('tymb' + str(i))
    Mymb = sm.var('Mymb' + str(i))
    Vymb = sm.var('Vymb' + str(i))

    wxnc = sm.var('wxnc' + str(i))
    txnc = sm.var('txnc' + str(i))
    Mxnc = sm.var('Mxnc' + str(i))
    Vxnc = sm.var('Vxnc' + str(i))
    wymc = sm.var('wymc' + str(i))
    tymc = sm.var('tymc' + str(i))
    Mymc = sm.var('Mymc' + str(i))
    Vymc = sm.var('Vymc' + str(i))

    wxnd = sm.var('wxnd' + str(i))
    txnd = sm.var('txnd' + str(i))
    Mxnd = sm.var('Mxnd' + str(i))
    Vxnd = sm.var('Vxnd' + str(i))
    wymd = sm.var('wymd' + str(i))
    tymd = sm.var('tymd' + str(i))
    Mymd = sm.var('Mymd' + str(i))
    Vymd = sm.var('Vymd' + str(i))

    wxne = sm.var('wxne' + str(i))
    txne = sm.var('txne' + str(i))
    Mxne = sm.var('Mxne' + str(i))
    Vxne = sm.var('Vxne' + str(i))
    wyme = sm.var('wyme' + str(i))
    tyme = sm.var('tyme' + str(i))
    Myme = sm.var('Myme' + str(i))
    Vyme = sm.var('Vyme' + str(i))

    wxna_slv = slv.args[0][i]
    txna_slv = slv.args[0][i + M]
    Mxna_slv = slv.args[0][i + 2 * M]
    Vxna_slv = slv.args[0][i + 3 * M]
    wyma_slv = slv.args[0][i + 4 * M]
    tyma_slv = slv.args[0][i + 5 * M]
    Myma_slv = slv.args[0][i + 6 * M]
    Vyma_slv = slv.args[0][i + 7 * M]

    wxnb_slv = slv.args[0][i + 8 * M]
    txnb_slv = slv.args[0][i + 9 * M]
    Mxnb_slv = slv.args[0][i + 10 * M]
    Vxnb_slv = slv.args[0][i + 11 * M]
    wymb_slv = slv.args[0][i + 12 * M]
    tymb_slv = slv.args[0][i + 13 * M]
    Mymb_slv = slv.args[0][i + 14 * M]
    Vymb_slv = slv.args[0][i + 15 * M]

    wxnc_slv = slv.args[0][i + 16 * M]
    txnc_slv = slv.args[0][i + 17 * M]
    Mxnc_slv = slv.args[0][i + 18 * M]
    Vxnc_slv = slv.args[0][i + 19 * M]
    wymc_slv = slv.args[0][i + 20 * M]
    tymc_slv = slv.args[0][i + 21 * M]
    Mymc_slv = slv.args[0][i + 22 * M]
    Vymc_slv = slv.args[0][i + 23 * M]

    wxnd_slv = slv.args[0][i + 24 * M]
    txnd_slv = slv.args[0][i + 25 * M]
    Mxnd_slv = slv.args[0][i + 26 * M]
    Vxnd_slv = slv.args[0][i + 27 * M]
    wymd_slv = slv.args[0][i + 28 * M]
    tymd_slv = slv.args[0][i + 29 * M]
    Mymd_slv = slv.args[0][i + 30 * M]
    Vymd_slv = slv.args[0][i + 31 * M]

    wxne_slv = slv.args[0][i + 32 * M]
    txne_slv = slv.args[0][i + 33 * M]
    Mxne_slv = slv.args[0][i + 34 * M]
    Vxne_slv = slv.args[0][i + 35 * M]
    wyme_slv = slv.args[0][i + 36 * M]
    tyme_slv = slv.args[0][i + 37 * M]
    Myme_slv = slv.args[0][i + 38 * M]
    Vyme_slv = slv.args[0][i + 39 * M]

    wa = wa.subs([
        (wxna, wxna_slv),
        (txna, txna_slv),
        (Mxna, Mxna_slv),
        (Vxna, Vxna_slv),
        (wyma, wyma_slv),
        (tyma, tyma_slv),
        (Myma, Myma_slv),
        (Vyma, Vyma_slv)
    ])
    Mxa = Mxa.subs([
        (wxna, wxna_slv),
        (txna, txna_slv),
        (Mxna, Mxna_slv),
        (Vxna, Vxna_slv),
        (wyma, wyma_slv),
        (tyma, tyma_slv),
        (Myma, Myma_slv),
        (Vyma, Vyma_slv)
    ])
    Mya = Mya.subs([
        (wxna, wxna_slv),
        (txna, txna_slv),
        (Mxna, Mxna_slv),
        (Vxna, Vxna_slv),
        (wyma, wyma_slv),
        (tyma, tyma_slv),
        (Myma, Myma_slv),
        (Vyma, Vyma_slv)
    ])
    Vxa = Vxa.subs([
        (wxna, wxna_slv),
        (txna, txna_slv),
        (Mxna, Mxna_slv),
        (Vxna, Vxna_slv),
        (wyma, wyma_slv),
        (tyma, tyma_slv),
        (Myma, Myma_slv),
        (Vyma, Vyma_slv)
    ])
    Vya = Vya.subs([
        (wxna, wxna_slv),
        (txna, txna_slv),
        (Mxna, Mxna_slv),
        (Vxna, Vxna_slv),
        (wyma, wyma_slv),
        (tyma, tyma_slv),
        (Myma, Myma_slv),
        (Vyma, Vyma_slv)
    ])
    wb = wb.subs([
        (wxnb, wxnb_slv),
        (txnb, txnb_slv),
        (Mxnb, Mxnb_slv),
        (Vxnb, Vxnb_slv),
        (wymb, wymb_slv),
        (tymb, tymb_slv),
        (Mymb, Mymb_slv),
        (Vymb, Vymb_slv)
    ])
    txb = txb.subs([
        (wxnb, wxnb_slv),
        (txnb, txnb_slv),
        (Mxnb, Mxnb_slv),
        (Vxnb, Vxnb_slv),
        (wymb, wymb_slv),
        (tymb, tymb_slv),
        (Mymb, Mymb_slv),
        (Vymb, Vymb_slv)
    ])
    tyb = tyb.subs([
        (wxnb, wxnb_slv),
        (txnb, txnb_slv),
        (Mxnb, Mxnb_slv),
        (Vxnb, Vxnb_slv),
        (wymb, wymb_slv),
        (tymb, tymb_slv),
        (Mymb, Mymb_slv),
        (Vymb, Vymb_slv)
    ])
    Mxb = Mxb.subs([
        (wxnb, wxnb_slv),
        (txnb, txnb_slv),
        (Mxnb, Mxnb_slv),
        (Vxnb, Vxnb_slv),
        (wymb, wymb_slv),
        (tymb, tymb_slv),
        (Mymb, Mymb_slv),
        (Vymb, Vymb_slv)
    ])
    Myb = Myb.subs([
        (wxnb, wxnb_slv),
        (txnb, txnb_slv),
        (Mxnb, Mxnb_slv),
        (Vxnb, Vxnb_slv),
        (wymb, wymb_slv),
        (tymb, tymb_slv),
        (Mymb, Mymb_slv),
        (Vymb, Vymb_slv)
    ])
    Vxb = Vxb.subs([
        (wxnb, wxnb_slv),
        (txnb, txnb_slv),
        (Mxnb, Mxnb_slv),
        (Vxnb, Vxnb_slv),
        (wymb, wymb_slv),
        (tymb, tymb_slv),
        (Mymb, Mymb_slv),
        (Vymb, Vymb_slv)
    ])
    Vyb = Vyb.subs([
        (wxnb, wxnb_slv),
        (txnb, txnb_slv),
        (Mxnb, Mxnb_slv),
        (Vxnb, Vxnb_slv),
        (wymb, wymb_slv),
        (tymb, tymb_slv),
        (Mymb, Mymb_slv),
        (Vymb, Vymb_slv)
    ])
    wc = wc.subs([
        (wxnc, wxnc_slv),
        (txnc, txnc_slv),
        (Mxnc, Mxnc_slv),
        (Vxnc, Vxnc_slv),
        (wymc, wymc_slv),
        (tymc, tymc_slv),
        (Mymc, Mymc_slv),
        (Vymc, Vymc_slv)
    ])
    txc = txc.subs([
        (wxnc, wxnc_slv),
        (txnc, txnc_slv),
        (Mxnc, Mxnc_slv),
        (Vxnc, Vxnc_slv),
        (wymc, wymc_slv),
        (tymc, tymc_slv),
        (Mymc, Mymc_slv),
        (Vymc, Vymc_slv)
    ])
    tyc = tyc.subs([
        (wxnc, wxnc_slv),
        (txnc, txnc_slv),
        (Mxnc, Mxnc_slv),
        (Vxnc, Vxnc_slv),
        (wymc, wymc_slv),
        (tymc, tymc_slv),
        (Mymc, Mymc_slv),
        (Vymc, Vymc_slv)
    ])
    Mxc = Mxc.subs([
        (wxnc, wxnc_slv),
        (txnc, txnc_slv),
        (Mxnc, Mxnc_slv),
        (Vxnc, Vxnc_slv),
        (wymc, wymc_slv),
        (tymc, tymc_slv),
        (Mymc, Mymc_slv),
        (Vymc, Vymc_slv)
    ])
    Myc = Myc.subs([
        (wxnc, wxnc_slv),
        (txnc, txnc_slv),
        (Mxnc, Mxnc_slv),
        (Vxnc, Vxnc_slv),
        (wymc, wymc_slv),
        (tymc, tymc_slv),
        (Mymc, Mymc_slv),
        (Vymc, Vymc_slv)
    ])
    Vxc = Vxc.subs([
        (wxnc, wxnc_slv),
        (txnc, txnc_slv),
        (Mxnc, Mxnc_slv),
        (Vxnc, Vxnc_slv),
        (wymc, wymc_slv),
        (tymc, tymc_slv),
        (Mymc, Mymc_slv),
        (Vymc, Vymc_slv)
    ])
    Vyc = Vyc.subs([
        (wxnc, wxnc_slv),
        (txnc, txnc_slv),
        (Mxnc, Mxnc_slv),
        (Vxnc, Vxnc_slv),
        (wymc, wymc_slv),
        (tymc, tymc_slv),
        (Mymc, Mymc_slv),
        (Vymc, Vymc_slv)
    ])
    wd = wd.subs([
        (wxnd, wxnd_slv),
        (txnd, txnd_slv),
        (Mxnd, Mxnd_slv),
        (Vxnd, Vxnd_slv),
        (wymd, wymd_slv),
        (tymd, tymd_slv),
        (Mymd, Mymd_slv),
        (Vymd, Vymd_slv)
    ])
    txd = txd.subs([
        (wxnd, wxnd_slv),
        (txnd, txnd_slv),
        (Mxnd, Mxnd_slv),
        (Vxnd, Vxnd_slv),
        (wymd, wymd_slv),
        (tymd, tymd_slv),
        (Mymd, Mymd_slv),
        (Vymd, Vymd_slv)
    ])
    tyd = tyd.subs([
        (wxnd, wxnd_slv),
        (txnd, txnd_slv),
        (Mxnd, Mxnd_slv),
        (Vxnd, Vxnd_slv),
        (wymd, wymd_slv),
        (tymd, tymd_slv),
        (Mymd, Mymd_slv),
        (Vymd, Vymd_slv)
    ])
    Mxd = Mxd.subs([
        (wxnd, wxnd_slv),
        (txnd, txnd_slv),
        (Mxnd, Mxnd_slv),
        (Vxnd, Vxnd_slv),
        (wymd, wymd_slv),
        (tymd, tymd_slv),
        (Mymd, Mymd_slv),
        (Vymd, Vymd_slv)
    ])
    Myd = Myd.subs([
        (wxnd, wxnd_slv),
        (txnd, txnd_slv),
        (Mxnd, Mxnd_slv),
        (Vxnd, Vxnd_slv),
        (wymd, wymd_slv),
        (tymd, tymd_slv),
        (Mymd, Mymd_slv),
        (Vymd, Vymd_slv)
    ])
    Vxd = Vxd.subs([
        (wxnd, wxnd_slv),
        (txnd, txnd_slv),
        (Mxnd, Mxnd_slv),
        (Vxnd, Vxnd_slv),
        (wymd, wymd_slv),
        (tymd, tymd_slv),
        (Mymd, Mymd_slv),
        (Vymd, Vymd_slv)
    ])
    Vyd = Vyd.subs([
        (wxnd, wxnd_slv),
        (txnd, txnd_slv),
        (Mxnd, Mxnd_slv),
        (Vxnd, Vxnd_slv),
        (wymd, wymd_slv),
        (tymd, tymd_slv),
        (Mymd, Mymd_slv),
        (Vymd, Vymd_slv)
    ])
    we = we.subs([
        (wxne, wxne_slv),
        (txne, txne_slv),
        (Mxne, Mxne_slv),
        (Vxne, Vxne_slv),
        (wyme, wyme_slv),
        (tyme, tyme_slv),
        (Myme, Myme_slv),
        (Vyme, Vyme_slv)
    ])
    txe = txe.subs([
        (wxne, wxne_slv),
        (txne, txne_slv),
        (Mxne, Mxne_slv),
        (Vxne, Vxne_slv),
        (wyme, wyme_slv),
        (tyme, tyme_slv),
        (Myme, Myme_slv),
        (Vyme, Vyme_slv)
    ])
    tye = tye.subs([
        (wxne, wxne_slv),
        (txne, txne_slv),
        (Mxne, Mxne_slv),
        (Vxne, Vxne_slv),
        (wyme, wyme_slv),
        (tyme, tyme_slv),
        (Myme, Myme_slv),
        (Vyme, Vyme_slv)
    ])
    Mxe = Mxe.subs([
        (wxne, wxne_slv),
        (txne, txne_slv),
        (Mxne, Mxne_slv),
        (Vxne, Vxne_slv),
        (wyme, wyme_slv),
        (tyme, tyme_slv),
        (Myme, Myme_slv),
        (Vyme, Vyme_slv)
    ])
    Mye = Mye.subs([
        (wxne, wxne_slv),
        (txne, txne_slv),
        (Mxne, Mxne_slv),
        (Vxne, Vxne_slv),
        (wyme, wyme_slv),
        (tyme, tyme_slv),
        (Myme, Myme_slv),
        (Vyme, Vyme_slv)
    ])
    Vxe = Vxe.subs([
        (wxne, wxne_slv),
        (txne, txne_slv),
        (Mxne, Mxne_slv),
        (Vxne, Vxne_slv),
        (wyme, wyme_slv),
        (tyme, tyme_slv),
        (Myme, Myme_slv),
        (Vyme, Vyme_slv)
    ])
    Vye = Vye.subs([
        (wxne, wxne_slv),
        (txne, txne_slv),
        (Mxne, Mxne_slv),
        (Vxne, Vxne_slv),
        (wyme, wyme_slv),
        (tyme, tyme_slv),
        (Myme, Myme_slv),
        (Vyme, Vyme_slv)
    ])

print('Готово\n')

wwa = (wa * D) / (q * a1 ** 4)
wwb = (wb * D) / (q * a2 ** 4)
wwc = (wc * D) / (q * a3 ** 4)
wwd = (wd * D) / (q * a4 ** 4)
wwe = (we * D) / (q * a5 ** 4)

MMxa = Mxa / (q * a1 ** 2)
MMxb = Mxb / (q * a2 ** 2)
MMxc = Mxc / (q * a2 ** 2)
MMxd = Mxd / (q * a2 ** 2)
MMxe = Mxe / (q * a2 ** 2)

VVxa = Vxa / (q * a1)
VVxb = Vxb / (q * a2)
VVxc = Vxc / (q * a2)
VVxd = Vxd / (q * a2)
VVxe = Vxe / (q * a2)

wwwa = sm.zeros(M + 2, 5)
wwwb = sm.zeros(M + 2, 5)
wwwc = sm.zeros(M + 2, 5)
wwwd = sm.zeros(M + 2, 5)
wwwe = sm.zeros(M + 2, 5)
MMMxa = sm.zeros(M + 2, 1)
MMMxb = sm.zeros(M + 2, 1)
MMMxc = sm.zeros(M + 2, 1)
MMMxd = sm.zeros(M + 2, 1)
MMMxe = sm.zeros(M + 2, 1)
VVVxa = sm.zeros(M + 2, 1)
VVVxb = sm.zeros(M + 2, 1)
VVVxc = sm.zeros(M + 2, 1)
VVVxd = sm.zeros(M + 2, 1)
VVVxe = sm.zeros(M + 2, 1)
xai = sm.zeros(M + 2, 1)
xbi = sm.zeros(M + 2, 1)
xci = sm.zeros(M + 2, 1)
xdi = sm.zeros(M + 2, 1)
xei = sm.zeros(M + 2, 1)
yai = sm.zeros(M + 2, 5)
ybi = sm.zeros(M + 2, 5)
yci = sm.zeros(M + 2, 5)
ydi = sm.zeros(M + 2, 5)
yei = sm.zeros(M + 2, 5)

# Вывод значений прогиба на серединах пластин
print(wwa.subs([(x, a1 / 2), (y, b1 / 2)]))
print(wwb.subs([(x, a2 / 2), (y, b2 / 2)]))
print(wwc.subs([(x, a3 / 2), (y, b3 / 2)]))
print(wwd.subs([(x, a4 / 2), (y, b4 / 2)]))
print(wwe.subs([(x, a5 / 2), (y, b5 / 2)]))

print('Вывод графиков')

# Вывод прогиба
for i in range(M + 2):
    xai[i] = (a1 / (M + 1)) * i
    xbi[i] = (a2 / (M + 1)) * i
    xci[i] = (a3 / (M + 1)) * i
    xdi[i] = (a4 / (M + 1)) * i
    xei[i] = (a5 / (M + 1)) * i

    yai[i, 0] = 0
    ybi[i, 0] = 0
    yci[i, 0] = 0
    ydi[i, 0] = 0
    yei[i, 0] = 0

    yai[i, 1] = b1 / 4
    ybi[i, 1] = b2 / 4
    yci[i, 1] = b3 / 4
    ydi[i, 1] = b4 / 4
    yei[i, 1] = b5 / 4

    yai[i, 2] = b1 / 2
    ybi[i, 2] = b2 / 2
    yci[i, 2] = b3 / 2
    ydi[i, 2] = b4 / 2
    yei[i, 2] = b5 / 2

    yai[i, 3] = (b1 * 3) / 4
    ybi[i, 3] = (b2 * 3) / 4
    yci[i, 3] = (b3 * 3) / 4
    ydi[i, 3] = (b4 * 3) / 4
    yei[i, 3] = (b5 * 3) / 4

    yai[i, 4] = b1
    ybi[i, 4] = b2
    yci[i, 4] = b3
    ydi[i, 4] = b4
    yei[i, 4] = b5

    wwwa[i, 0] = wwa.subs([(x, xai[i]), (y, yai[i, 0])])
    wwwb[i, 0] = wwb.subs([(x, xbi[i]), (y, ybi[i, 0])])
    wwwc[i, 0] = wwc.subs([(x, xci[i]), (y, yci[i, 0])])
    wwwd[i, 0] = wwd.subs([(x, xdi[i]), (y, ydi[i, 0])])
    wwwe[i, 0] = wwe.subs([(x, xei[i]), (y, yei[i, 0])])

    wwwa[i, 1] = wwa.subs([(x, xai[i]), (y, yai[i, 1])])
    wwwb[i, 1] = wwb.subs([(x, xbi[i]), (y, ybi[i, 1])])
    wwwc[i, 1] = wwc.subs([(x, xci[i]), (y, yci[i, 1])])
    wwwd[i, 1] = wwd.subs([(x, xdi[i]), (y, ydi[i, 1])])
    wwwe[i, 1] = wwe.subs([(x, xei[i]), (y, yei[i, 1])])

    wwwa[i, 2] = wwa.subs([(x, xai[i]), (y, yai[i, 2])])
    wwwb[i, 2] = wwb.subs([(x, xbi[i]), (y, ybi[i, 2])])
    wwwc[i, 2] = wwc.subs([(x, xci[i]), (y, yci[i, 2])])
    wwwd[i, 2] = wwd.subs([(x, xdi[i]), (y, ydi[i, 2])])
    wwwe[i, 2] = wwe.subs([(x, xei[i]), (y, yei[i, 2])])

    wwwa[i, 3] = wwa.subs([(x, xai[i]), (y, yai[i, 3])])
    wwwb[i, 3] = wwb.subs([(x, xbi[i]), (y, ybi[i, 3])])
    wwwc[i, 3] = wwc.subs([(x, xci[i]), (y, yci[i, 3])])
    wwwd[i, 3] = wwd.subs([(x, xdi[i]), (y, ydi[i, 3])])
    wwwe[i, 3] = wwe.subs([(x, xei[i]), (y, yei[i, 3])])

    wwwa[i, 4] = wwa.subs([(x, xai[i]), (y, yai[i, 4])])
    wwwb[i, 4] = wwb.subs([(x, xbi[i]), (y, ybi[i, 4])])
    wwwc[i, 4] = wwc.subs([(x, xci[i]), (y, yci[i, 4])])
    wwwd[i, 4] = wwd.subs([(x, xdi[i]), (y, ydi[i, 4])])
    wwwe[i, 4] = wwe.subs([(x, xei[i]), (y, yei[i, 4])])
    

# Смещение графиков для совместного вывода 1 + 2 + 4
for i in range(M+2):
    xbi[i] += a1
    xdi[i] += a1 + a2

# График 1 + 2 + 4, 3, 5
plt.subplot(211)
plt.plot(xai, wwwa[:, 0], '-k', label='y = 0')
plt.plot(xbi, wwwb[:, 0], '-k')
plt.plot(xdi, wwwd[:, 0], '-k')
plt.plot(xai, wwwa[:, 1], '--k', label='y = b/4')
plt.plot(xbi, wwwb[:, 1], '--k')
plt.plot(xdi, wwwd[:, 1], '--k')
plt.plot(xai, wwwa[:, 2], ':k', label='y = b/2')
plt.plot(xbi, wwwb[:, 2], ':k')
plt.plot(xdi, wwwd[:, 2], ':k')
plt.plot(xai, wwwa[:, 3], '-.k', label='y = 3b/4')
plt.plot(xbi, wwwb[:, 3], '-.k')
plt.plot(xdi, wwwd[:, 3], '-.k')
plt.plot(xai, wwwa[:, 4], '-k', label='y = b')
plt.plot(xbi, wwwb[:, 4], '-k')
plt.plot(xdi, wwwd[:, 4], '-k')
plt.title("1 + 2 + 4")
plt.legend()
plt.subplot(234)
plt.plot(xci, wwwc[:, 0], '-k', label='y = 0')
plt.plot(xci, wwwc[:, 1], '--k', label='y = b/4')
plt.plot(xci, wwwc[:, 2], ':k', label='y = b/2')
plt.plot(xci, wwwc[:, 3], '-.k', label='y = 3b/4')
plt.plot(xci, wwwc[:, 4], '-k', label='y = b')
plt.title("3")
plt.legend()
plt.subplot(236)
plt.plot(xci, wwwe[:, 0], '-k', label='y = 0')
plt.plot(xci, wwwe[:, 1], '--k', label='y = b/4')
plt.plot(xci, wwwe[:, 2], ':k', label='y = b/2')
plt.plot(xci, wwwe[:, 3], '-.k', label='y = 3b/4')
plt.plot(xci, wwwe[:, 4], '-k', label='y = b')
plt.title("5")
plt.legend()
plt.show()


"""
# Вывод момента
for i in range(M + 2):

    xai[i, 0] = 0
    yai[i, 0] = (b1 / (M + 1)) * i
    xci[i, 0] = 0
    yci[i, 0] = (b3 / (M + 1)) * i

    xdi[i, 0] = a4
    ydi[i, 0] = (b4 / (M + 1)) * i
    xei[i, 0] = a5
    yei[i, 0] = (b5 / (M + 1)) * i

    MMMxa[i] = MMxa.subs([(x, xai[i, 0]), (y, yai[i, 0])])
    MMMxc[i] = MMxc.subs([(x, xci[i, 0]), (y, yci[i, 0])])

    MMMxd[i] = MMxd.subs([(x, xdi[i, 0]), (y, ydi[i, 0])])
    MMMxe[i] = MMxe.subs([(x, xei[i, 0]), (y, yei[i, 0])])

# Смещение графиков для совместного вывода 1 + 3 и 4 + 5
for i in range(M+2):
    yci[i] += b1
    yei[i] += b4

# График 1 + 3, 4 + 5
plt.subplot(221)
plt.plot(yai[:, 0], MMMxa)
plt.subplot(222)
plt.plot(yci[:, 0], MMMxc)
#plt.title("1 + 3")
plt.subplot(223)
plt.plot(ydi[:, 0], MMMxd)
plt.subplot(224)
plt.plot(yei[:, 0], MMMxe)
#plt.title("4 + 5")
plt.show()
"""
"""
# Вывод перерезывающей силы
for i in range(M + 2):
    xai[i] = 0
    yai[i] = (b1 / (M + 1)) * i
    xci[i] = 0
    yci[i] = (b3 / (M + 1)) * i

    xdi[i] = a4
    ydi[i] = (b4 / (M + 1)) * i
    xei[i] = a5
    yei[i] = (b5 / (M + 1)) * i

    VVVxa[i] = VVxa.subs([(x, xai[i]), (y, yai[i])])
    VVVxc[i] = VVxc.subs([(x, xci[i]), (y, yci[i])])

    VVVxd[i] = VVxd.subs([(x, xdi[i]), (y, ydi[i])])
    VVVxe[i] = VVxe.subs([(x, xei[i]), (y, yei[i])])

# Смещение графиков для совместного вывода 1 + 3 и 4 + 5
for i in range(M+2):
    yci[i] += b1
    yei[i] += b4

# График 1 + 3, 4 + 5
plt.subplot(121)
plt.plot(yai, VVVxa)
plt.plot(yci, VVVxc)
plt.title("1 + 3")
plt.subplot(122)
plt.plot(ydi, VVVxd)
plt.plot(yei, VVVxe)
plt.title("4 + 5")
plt.show()
"""