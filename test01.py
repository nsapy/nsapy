# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis

'''
from nsapy import Domain,Analysis,PostProcessor
from nsapy.Material import Material1D
from nsapy import Section,Element

def main():

    truss = Domain(2,2)

    truss.add_mat(Material1D.Material_Elastic,1,2.e5)
    truss.add_sec(Section.Section_Elastic,1,truss.mat[1],100.)

    L = 1000.0
    
    truss.add_node(1,(0.0,0.,0.),200.)
    truss.add_node(2,(L  ,0.,0.),200.)
        
    truss.add_ele(Element.Truss2D,1,truss.sec[1],(truss.node[1],truss.node[2]))
    F = -1.e4
    truss.add_load(1,truss.node[2],[1],[F])

    truss.add_cons(1,truss.node[1],[1,2],[0.0,0.0])
    truss.add_cons(2,truss.node[2],[2],[0.0])

    truss.build_model()
    # truss.export_gmsh_scr_2d('truss')

    analyse = Analysis.Analysis_Linear_Eigen(truss)
    analyse.execute(1)
    analyse = Analysis.Analysis_Linear_Static(truss)
    analyse.execute()

    post = PostProcessor(truss)
    post.get_all()
    
    return truss

if __name__ == '__main__':

    # main = autojit(main)
    truss = main()

    print 'T = ',truss.T[0],'s.'
    print 'Element Disp = ',truss.ele[1].local_deformation[0],'mm.'