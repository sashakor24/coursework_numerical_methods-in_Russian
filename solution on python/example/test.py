from math import pi, sin, cos
import numpy as np
from scipy.integrate import RK45
import matplotlib.pyplot as plt

def solution(t): #Точное решение ДУ
    return sin(pi*t)

def f(t,x):
    return pi*cos(pi*t)+x-sin(pi*t)


def Euler(t,x,h,f):
    for i in range(x.size - 1):
        x[i+1] = x[i] + h*f(t[i],x[i])

def RK4(t,x,h,f):
    for i in range(x.size - 1):
        k1=f(t[i],x[i])
        k2=f(t[i]+h/2,x[i]+k1*h/2)
        k3=f(t[i]+h/2,x[i]+k2*h/2)
        k4=f(t[i]+h,x[i]+k3*h)
        x[i+1] = x[i] + (h/6)*(k1+2*k2+2*k3+k4)

def CError(x1,x2):
    return max([abs(x1[i]-x2[i]) for i in range(x1.size)])/max([abs(x) for x in x2])

def LError(x1,x2):
    return sum([abs(x1[i]-x2[i]) for i in range(x1.size)])/sum([abs(x) for x in x2])

def do(N):
    h = 1/N
    t = np.array([i*h for i in range(N+1)]) #сетка
    x_eu = np.zeros(N+1,float)
    x_rk = np.zeros(N+1,float)
    x_lib = np.zeros(N+1,float)
    Euler(t,x_eu,h,f)
    RK4(t,x_rk,h,f)

    rks = RK45(f,0,[0],1)#(f,x[0],[y0],x[N])
    rks.first_step = h
    rks.max_step = h
    rks.min_step = h
    rks.h_abs = h
    rks.atol = 1
    rks.rtol = 1
    x_lib[0]= rks.y[0]
    i = 1

    while(True):
        rks.step()
        if 1 - rks.t >= h/2: # t[N]- текущая позиция
            x_lib[i] = rks.y[0]
            i += 1
        elif rks.status != 'finished':
            rks.step()
            x_lib[i] = rks.y[0]
            i += 1
            break
        else:
            x_lib[i] = rks.y[0]
            i += 1
            break

    s = np.array([solution(t[i]) for i in range(N+1)]) #точное решение на узлах сетки
    with open('test.txt', 'a') as file:
        file.write('N = '+str(N)+', ')
        file.write('C_eu: '+str(CError(s,x_eu))+', ')
        file.write('L_eu: '+str(LError(s,x_eu))+', ')
        file.write('C_rk: '+str(CError(s,x_rk))+', ')
        file.write('L_rk: '+str(LError(s,x_rk))+', ')
        file.write('C_lib: '+str(CError(s,x_lib))+', ')
        file.write('L_lib: '+str(LError(s,x_lib))+';\n')
    plt.plot(t,x_eu, 'g') #зелёный график - график численного решения методом Эйлера
    plt.plot(t,x_rk, 'b') #синий график - график численного решения методом Рунге-Кутты
    plt.plot(t,x_lib, 'r') #красный график - график численного решения методом Рунге-Кутты с использованием библиотеки
    if N<100: #чтобы точное решение более гладко выглядело, построим его на 101 узле, если N<100
        t = np.array([i/100 for i in range(101)])
        s = np.array([solution(t[i]) for i in range(101)])
    plt.plot(t,s, 'brown') #коричневый график - график, построенный соединением точных значений в точках сетки
    plt.savefig('тест_'+str(N)+'.png', dpi = 300)
    plt.close()
