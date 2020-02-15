
import math
import sympy
from sympy import *

class DensityError():

    def __init__(self,ps,r,t,ps_dsa,ps_dsa_error,t_error):
        self.ps = ps
        self.r = r
        self.t = t + 273.15
        self.ps_dsa = ps_dsa * 703.0695796
        self.ps_dsa_error = ps_dsa_error/100
        self.t_error = t_error
        
    def density(self):
        density = self.ps / self.r / (self.t)
        return density

    def ps_error_abs(self):
        ps_error_abs = (self.ps_dsa) * (self.ps_dsa_error)
        return ps_error_abs

    def diff(self):
        P = sympy.Symbol('P')
        R = sympy.Symbol('R')
        T = sympy.Symbol('T')
        ρ = P / R / T 
        dρdP = sympy.diff(ρ,P)
        dρdR = sympy.diff(ρ,R)
        dρdT = sympy.diff(ρ,T)
        print('dρ/dP= '+ str(dρdP))
        print('dρ/dR= '+ str(dρdR))
        print('dρ/dT= '+ str(dρdT))



class AreaError():
    def __init__(self,diameter,diameter_error):
        self.diameter = diameter /1000
        self.diameter_error = diameter_error /1000

    def area(self):
        diameter = self.diameter
        area = math.pi * diameter**2 /4
        return area

    def diameter_error_abs(self):
        diameter_error = self.diameter_error
        return diameter_error



class AlphaError():
    def __init__(self,alpha,alpha_error):
        self.alpha = alpha
        self.alpha_error = alpha_error/100

    def alpha_error_abs(self):
        alpha_error_abs = self.alpha * self.alpha_error
        return alpha_error_abs


class VelocityError():
    def __init__(self,dynamicpressure,pv_dsa,pv_dsa_error,density,density_error):
        self.dynamicpressure = dynamicpressure
        self.pv_dsa = pv_dsa*703.0695796
        self.pv_dsa_error =pv_dsa_error /100
        self.density = density
        self.density_error =density_error

    def pv_error_abs(self):
        pv_error_abs = self.pv_dsa * self.pv_dsa_error
        return pv_error_abs



class MassError():
    def mass_diff(self):
        α = sympy.Symbol('α')
        ρ = sympy.Symbol('ρ')
        A = sympy.Symbol('A')
        Pv = sympy.Symbol('Pv')
        G = α * ρ * A * sympy.sqrt(2 * Pv / ρ) 
        dGdα = sympy.diff(G,α)
        dGdρ = sympy.diff(G,ρ)
        dGdA = sympy.diff(G,A)
        dGdPv = sympy.diff(G,Pv)
        print('dG/dα= '+ str(dGdα))
        print('dG/dρ= '+ str(dGdρ))
        print('dG/dA= '+ str(dGdA))
        print('dG/dPv= '+ str(dGdPv))


if __name__ == '__main__':
    #静圧Ps (mmAq)
    #ガス定数R　(m/k)
    #温度T (℃)
    #DSAレンジ(psi)
    #DSA誤差 (%)
    #温度誤差(℃)


    ps = 11050
    r = 29.3
    t = 15
    ps_dsa = 50
    ps_dsa_error = 0.12
    t_error = 0.3

    print()
    z = DensityError(ps, r, t, ps_dsa, ps_dsa_error ,t_error) 
    density_error = math.sqrt(z.ps_error_abs()**2 * (1/(z.r*z.t))**2 + \
        z.t_error**2 * ((-1) * z.ps/(z.r * z.t**2))**2)
    print('密度 '+ str(round(z.density(),3)) + '[kg/s]')
    print('壁圧DSAのレンジ ' +str(ps_dsa) +'[psi]')
    print('DSAの誤差 ' +str(ps_dsa_error) +'[%F.S]')
    print('DSAの絶対誤差 ' +str(round(z.ps_error_abs(),3)) +'[mmAq]')
    z.diff()
    print('密度誤差は ' + str(round(density_error,4)) + '[kg/m3]')


    #配管直径 mm
    #直径誤差　mm

    d = 500
    d_error = 5
    print()
    c = AreaError(d ,d_error)
    print('配管直径 '+ str(d) + '[m]')
    print('配管直径誤差 '+ str(round(c.diameter_error_abs(),3)) + '[m]')


    #流量係数α
    #流量計数誤差[%]
    a = 1.0
    a_error = 0.5
    print()
    a = AlphaError(a ,a_error)
    print('流量係数α '+ str(a))
    print('流量係数誤差 '+ str(round(a.alpha_error_abs(),3)))


    
    #動圧(mmAq)
    #DSAレンジ (psi)
    #DSA誤差　(%)
    pv = 15
    pv_dsa = 1
    pv_dsa_error = 0.12
    print()
    v = VelocityError(pv, pv_dsa_error, pv_dsa_error, z.density(), density_error)
    print('差圧DSAのレンジ '+ str(pv_dsa)+'[psi]')
    print('DSAの誤差 '+ str(pv_dsa_error)+'[%F.S]')

    print()
    print('流量の計算式を微分')
    e = MassError()
    e.mass_diff()

    a_sensitivity = a.alpha_error_abs()**2 \
        * (math.sqrt(2)*c.area()*math.sqrt(pv/10000*98066.5/z.density()) )**2

    density_sensitivity = density_error()**2 \
        * (math.sqrt(2)*c.area()*a)

    mass_error = math.sqrt()
    print(mass_error)
