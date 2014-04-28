# -*- coding:utf-8 -*-

from sympy import *
from sympy.abc import x,L,theta

init_printing()

import pyperclip as pp

Fa = symbols('F_a')
Da = symbols('D_a')
D = Matrix([Da])

B0 = Matrix([1/L,0]).T
Temp = B0.T*D*B0
K0 = Temp.integrate((x,0,L))

BL = theta*Matrix([0,1/L]).T
Temp = B0.T*D*BL+BL.T*D*BL+BL.T*D*B0
KL = Temp.integrate((x,0,L))

dBL = Matrix([0,1/L])*Matrix([0,1/L]).T
Temp = dBL*Fa
KG = Temp.integrate((x,0,L))

K = K0+KL+KG

B0 = Matrix([1/L,0,0]).T
Temp = B0.T*D*B0
K0 = Temp.integrate((x,0,L))

thetav = symbols('theta_v')
thetaw = symbols('theta_w')
BL = thetav*Matrix([0,1/L,0]).T+thetaw*Matrix([0,0,1/L]).T
Temp = B0.T*D*BL+BL.T*D*BL+BL.T*D*B0
KL = Temp.integrate((x,0,L))

dBL = Matrix([0,1/L,0])*Matrix([0,1/L,0]).T+Matrix([0,0,1/L])*Matrix([0,0,1/L]).T
Temp = dBL*Fa
KG = Temp.integrate((x,0,L))

K = K0+KL+KG

