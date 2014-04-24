# -*- coding:utf-8 -*-
'''
**Nsapy ---- Python Module for Nonlinear Structural Analysis

This is a testing script.

'''
from nsapy import Domain,Analysis,PostProcessor
from nsapy.Material import Material1D
from nsapy import Section,Element

def main():

    truss = Domain(2,2)

    truss.add_mat(Material1D.Material_Elastic,1,1.)
    truss.add_sec(Section.Section_Elastic,1,truss.mat[1],1.)

    d = 1000.0

    for i in range(7):
        truss.add_node(11+i,(i*d,0,0),1.)
        truss.add_node(21+i,(i*d,d,0),1.)
        truss.add_ele(Element.Truss2D,'V-%d'%(1+i),truss.sec[1],(truss.node[11+i],truss.node[21+i]))
        truss.add_load(1+i,truss.node[21+i],[2],[-1.])

    truss.load[1].dof_value[0] = -.5
    truss.load[7].dof_value[0] = -.5

    truss.add_cons(1,truss.node[11],[1,2],[0.0,0.0])
    truss.add_cons(2,truss.node[17],[1,2],[0.0,0.0])

    for i in range(6):
        truss.add_ele(Element.Truss2D,'H1-%d'%(1+i),truss.sec[1],(truss.node[11+i],truss.node[12+i]))
        truss.add_ele(Element.Truss2D,'H2-%d'%(1+i),truss.sec[1],(truss.node[21+i],truss.node[22+i]))
        if i < 3:
            truss.add_ele(Element.Truss2D,'DL-%d'%(1+i),truss.sec[1],(truss.node[12+i],truss.node[21+i]))
        else:
            truss.add_ele(Element.Truss2D,'DR-%d'%(1+i),truss.sec[1],(truss.node[11+i],truss.node[22+i]))     

    truss.build_model()
    truss.export_gmsh_scr_2d('truss')

    analyse = Analysis.Analysis_Linear_Static(truss)
    analyse.execute()

    post = PostProcessor(truss)
    post.get_all()
    
    return truss

if __name__ == '__main__':

    # main = autojit(main)
    truss = main()
