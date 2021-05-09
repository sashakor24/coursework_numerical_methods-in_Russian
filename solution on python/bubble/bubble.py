from math import pi, sin, cos
import numpy as np
from scipy.integrate import RK45
import matplotlib.pyplot as plt

p_inf = 1 #атм
gamma = 1.4
C0 = 1.5e5 #см/сек
p0 = 0.01 #атм
R0 = 1 #см эта строчка нужна только, чтобы создать глобальную переменную,
# R0 инициализируется в функции "do"
n = 7
B = 3214 #атм
rho0 = 1 #kg/см^3

def p(R):
    return p0*(R0/R)**(3*gamma)
    
def F(R):
    return 1 + (p(R)-p_inf)/B

def f2(t,xi):
    return ( n*B/(2*xi[0]*rho0*(n-1)) * (F(xi[0])**(1-1/n)-1) - (xi[1]**2)/xi[0] -
    (3*gamma*R0**(3*gamma)*p0)/(C0*rho0) * F(xi[0])**(-1/n) * xi[1]/(xi[0]**(3*gamma-1)) )

def f(t,xi):#f(t - float,xi[0] - float, xi[1] - float)
    return np.array([xi[1],f2(t,xi)])

def Euler(t,x,h,f):
    for i in range(x.shape[0] - 1):
        x[i+1] = x[i] + h*f(t[i],x[i])
        
def RK4(t,x,h,f):
    for i in range(x.shape[0] - 1):
        k1=f(t[i],x[i])
        k2=f(t[i]+h/2,x[i]+k1*h/2)
        k3=f(t[i]+h/2,x[i]+k2*h/2)
        k4=f(t[i]+h,x[i]+k3*h)
        x[i+1] = x[i] + (h/6)*(k1+2*k2+2*k3+k4)

def do(N, R, bEu, bRk, bLib):# eu - зелёный, rk - синий, lib - красный
    global R0
    R0 = R
    h = 1/N
    t = np.array([i*h for i in range(N+1)]) #сетка
    # x[i]=[R[i],M[i]]=[x[i][0],x[i][1]]
    if bEu:
        x_eu = np.zeros((N+1,2),float) # инициализируем нулями
        x_eu[0][0] = R0 # задача Коши: R[0]=R0, M[0] = 0
        Euler(t,x_eu,h,f)
        plt.plot(t,x_eu[:,0], 'g') # строим график
    if bRk:
        x_rk = np.zeros((N+1,2),float)
        x_rk[0][0] = R0
        RK4(t,x_rk,h,f)
        plt.plot(t,x_rk[:,0], 'b')
    if bLib:
        x_lib = np.zeros((N+1,2),float)
        x_lib[0][0] = R0
        rks = RK45(f,t[0],x_lib[0],t[N])#(f,t[0],[x[0][1],x[0][2]],t[N])
        rks.first_step = h
        rks.max_step = h
        rks.min_step = h
        rks.h_abs = h
        rks.atol = 1
        rks.rtol = 1
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

        plt.plot(t, x_lib[:,0], 'r')
        
    plt.savefig('bubble_N='+str(N)+'_R0='+str(R0)+'.png', dpi = 300) # сохраняем в файл
    plt.close()
