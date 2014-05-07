# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis

'''
from nsapy import Domain,Analysis,NodeResult
from nsapy.Material import Material1D
from nsapy import Section,Element
import numpy as np
import pylab as pl

def main():

    frame = Domain(2,3,'test07')

    frame.add_mat(Material1D.Bilinear,1,2.e5,200.,0.01)

    b = 200.
    h = 400.
    A = b*h
    I = b*h**3/12.
    frame.add_sec(Section.FiberSection2D,1)
    frame.sec[1].add_rect(b,h,nb=1,nh=10,mat=frame.mat[1])

    L = 3000.0
    nmass = [200.0,200.0,1.e-9]
    
    frame.add_node(1,(0.0,0.,0.),nmass)
    frame.add_node(2,(L  ,0.,0.),nmass)
    nintp = 10

    frame.add_ele(Element.Frame2D,1,frame.sec[1],(frame.node[1],frame.node[2]),nintp)

    N = 10
    t = np.linspace(1./N,1.,N)
    F1 = -700.e3*t
    
    # F1 = -8.e3*np.sin(2.*np.pi/4.*t)
    # F2 = -6.e3*np.sin(2.*np.pi/4.*t)

    frame.add_load(1,frame.node[2],[2],[F1],t)
    frame.add_cons(1,frame.node[1],[1,2,3],[0.0,0.0,0.0])

    frame.build_model()
    # frame.apply_cons(0)

    # analyse = Analysis.Analysis_Linear_Eigen(frame)
    # analyse.execute(1)

    # analyse = Analysis.Analysis_Linear_Static(frame)
    analyse = Analysis.Analysis_NonLinear_Static(frame)
    analyse.tol = 1e-6
    analyse.maxiter = 10

    disp = NodeResult(frame,2,'deformation',2)
    react = NodeResult(frame,1,'reaction',2)
    analyse.results.append(disp)
    analyse.results.append(react)

    # analyse.execute(6)
    analyse.execute(len(t))
    
    u = F1*L**3/2.e5/I/3.
    
    # print '            Disp          Moment'
    # print 'nsapy    %.6f    %.6f'%(frame.U[-2], -frame.UF[2])
    # print 'Theory   %.6f    %.6f'%( u[-1], F1[-1]*L)

    for i in xrange(len(t)):
        print t[i],disp.values[i]

    pl.plot(disp.values,react.values,'-')
    pl.plot(disp.values,-F1,'-')
    pl.grid(1)
    pl.show()
    return frame

if __name__ == '__main__':

    frame = main()