# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis

'''
from nsapy import Domain,Analysis,NodeResult
from nsapy.Material import Material1D
from nsapy import Section,Element
import numpy as np
from pylab import *
from numba import autojit

def main():

    truss = Domain(2,2,'test04')

    truss.add_mat(Material1D.Bilinear,1,2.e5,200.,0.02)
    truss.add_sec(Section.Truss,1,truss.mat[1],100.)

    L = 1000.0
    mass = [200.,200.,200.]
    
    truss.add_node(1,(0.0,0.,0.),mass)
    truss.add_node(2,(L  ,0.,0.),mass)
    
    truss.add_ele(Element.Truss2D,1,truss.sec[1],(truss.node[1],truss.node[2]))
    dt = 0.01
    t = np.linspace(dt,1.,int(1./dt))
    F = -2.e4*1.2*np.sin(2.*np.pi*t)
    truss.add_load(1,truss.node[2],[1],[F],t)

    truss.add_cons(1,truss.node[1],[1,2],[0.0,0.0])
    truss.add_cons(2,truss.node[2],[2],[0.0])

    truss.build_model()
    # truss.apply_cons(0)

    # analyse = Analysis.Analysis_Eigen(truss)
    # analyse.execute(1)
    
    analyse = Analysis.Analysis_NonLinear_Static(truss)
    disp = NodeResult(truss,2,'deformation',1)
    react = NodeResult(truss,1,'reaction',1)
    analyse.results.append(disp)
    analyse.results.append(react)

    analyse.execute(len(F))

    # print 'T = ',truss.T[0],'s.'
    plot(disp.values,react.values,'-o');grid(1)
    show()

    return truss

if __name__ == '__main__':

    # main = autojit(main)
    truss = main()