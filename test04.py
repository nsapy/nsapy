# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis

'''
from nsapy import Domain,Analysis,PostProcessor
from nsapy.Material import Material1D
from nsapy import Section,Element
import numpy as np

def main():

    truss = Domain(2,2,'test04')

    truss.add_mat(Material1D.Bilinear,1,2.e5,200.,0.01)
    truss.add_sec(Section.Section_Elastic_Truss,1,truss.mat[1],100.)

    L = 1000.0
    
    truss.add_node(1,(0.0,0.,0.),200.)
    truss.add_node(2,(L  ,0.,0.),200.)
        
    truss.add_ele(Element.Truss2D,1,truss.sec[1],(truss.node[1],truss.node[2]))
    F = -2.1e4*np.linspace(0.1,1,10)
    truss.add_load(1,truss.node[2],[1],[F])

    truss.add_cons(1,truss.node[1],[1,2],[0.0,0.0])
    truss.add_cons(2,truss.node[2],[2],[0.0])

    truss.build_model()
    truss.apply_cons()
    truss.write_gmsh()

    analyse = Analysis.Analysis_Eigen(truss)
    analyse.execute(1)
    analyse = Analysis.Analysis_Linear_Static(truss)
    analyse.execute(10)
    
    return truss

if __name__ == '__main__':

    truss = main()

    print 'T = ',truss.T[0],'s.'
    print 'Element Disp = ',truss.ele[1].local_deformation[0],'mm.'