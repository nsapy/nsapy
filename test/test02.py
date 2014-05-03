# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis

'''
from nsapy import Domain,Analysis
from nsapy.Material import Material1D
from nsapy import Section,Element
import numpy as np

def main():

    truss = Domain(2,2,'test02')

    truss.add_mat(Material1D.Material_Elastic,1,2.e5)
    # truss.add_mat(Material1D.Bilinear,1,2.e5,200.,0.02)
    truss.add_sec(Section.Elastic_Truss,1,truss.mat[1],100.)

    L = 4500.0
    l1 = 1500.0
    l2 = 750.0
    h = 500.0
    nmass = [200.0,200.0,200.0]
    
    truss.add_node(1,(0.0,0.,0.),nmass)
    truss.add_node(2,(L,0.,0.),nmass)
    truss.add_node(3,(l1,0.,0.),nmass)
    truss.add_node(4,(l1+l2,0.,0.),nmass)
    truss.add_node(5,(l1+2*l2,0.,0.),nmass)
    truss.add_node(6,(l1,-h,0.),nmass)
    truss.add_node(7,(l1+2*l2,-h,0.),nmass)

    truss.add_ele(Element.Truss2D,1,truss.sec[1],(truss.node[1],truss.node[3]))
    truss.add_ele(Element.Truss2D,2,truss.sec[1],(truss.node[3],truss.node[4]))
    truss.add_ele(Element.Truss2D,3,truss.sec[1],(truss.node[4],truss.node[5]))
    truss.add_ele(Element.Truss2D,4,truss.sec[1],(truss.node[5],truss.node[2]))
    truss.add_ele(Element.Truss2D,5,truss.sec[1],(truss.node[1],truss.node[6]))
    truss.add_ele(Element.Truss2D,6,truss.sec[1],(truss.node[6],truss.node[7]))
    truss.add_ele(Element.Truss2D,7,truss.sec[1],(truss.node[7],truss.node[2]))
    truss.add_ele(Element.Truss2D,8,truss.sec[1],(truss.node[3],truss.node[6]))
    truss.add_ele(Element.Truss2D,9,truss.sec[1],(truss.node[6],truss.node[4]))
    truss.add_ele(Element.Truss2D,10,truss.sec[1],(truss.node[4],truss.node[7]))
    truss.add_ele(Element.Truss2D,11,truss.sec[1],(truss.node[7],truss.node[5]))
    

    t = np.linspace(0.01,1,100)
    F1 = -8.e3*t
    F2 = -6.e3*t
    
    # F1 = -8.e3*np.sin(2.*np.pi/4.*t)
    # F2 = -6.e3*np.sin(2.*np.pi/4.*t)

    truss.add_load(1,truss.node[3],[2],[F1],t-0.01)
    truss.add_load(2,truss.node[4],[2],[F2],t-0.01)
    truss.add_load(3,truss.node[5],[2],[F1],t-0.01)

    truss.add_cons(1,truss.node[1],[1,2],[0.0,0.0])
    truss.add_cons(2,truss.node[2],[2],[0.0,])

    truss.build_model()
    # truss.apply_cons(0)

    # analyse = Analysis.Analysis_Linear_Eigen(truss)
    # analyse.execute(1)
    analyse = Analysis.Analysis_Linear_Static(truss)
    # analyse = Analysis.Analysis_NonLinear_Static(truss)
    analyse.execute(len(t))
    
    return truss

if __name__ == '__main__':

    # main = autojit(main)
    truss = main()

    # print 'T = ',truss.T[0],'s.'
    # print 'Element Disp = ',truss.ele[1].local_deformation[0],'mm.'