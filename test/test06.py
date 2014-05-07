# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis

'''
from nsapy import Domain,Analysis
from nsapy.Material import Material1D
from nsapy import Section,Element
import numpy as np

def main():

    frame = Domain(2,3,'test06')

    frame.add_mat(Material1D.Material_Elastic,1,2.e5)

    b = 200.
    h = 400.
    A = b*h
    I = b*h**3/12.
    frame.add_sec(Section.Elastic_Frame2D,1,frame.mat[1],A,I)

    L = 3000.0
    nmass = [200.0,200.0,0.0]
    
    frame.add_node(1,(0.0,0.,0.),nmass)
    frame.add_node(2,(L  ,0.,0.),nmass)

    frame.add_ele(Element.Frame2D_Elastic,1,frame.sec[1],(frame.node[1],frame.node[2]))

    t = np.linspace(0.1,1,10)
    F1 = -10.e3*t
    
    # F1 = -8.e3*np.sin(2.*np.pi/4.*t)
    # F2 = -6.e3*np.sin(2.*np.pi/4.*t)

    frame.add_load(1,frame.node[2],[2],[F1],t-.1)
    frame.add_cons(1,frame.node[1],[1,2,3],[0.0,0.0,0.0])

    frame.build_model()
    # frame.apply_cons(0)

    # analyse = Analysis.Analysis_Linear_Eigen(frame)
    # analyse.execute(1)

    analyse = Analysis.Analysis_Linear_Static(frame)
    # analyse = Analysis.Analysis_NonLinear_Static(frame)
    analyse.execute(len(t))
    
    u = F1*L**3/2.e5/I/3.
    
    print '            Disp          Moment'
    print 'nsapy    %.6f    %.6f'%(frame.U[-2], -frame.UF[2])
    print 'Theory   %.6f    %.6f'%( u[-1], F1[-1]*L)
    
    return frame

if __name__ == '__main__':

    frame = main()