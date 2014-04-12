# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis
**用于结构非线性分析的Python模块

'''
from numba import jit,autojit
import numpy as np
import numpy.linalg as npl
import scipy.linalg as spl
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import matplotlib.pylab as pl

VECX = np.asarray([1.,0.,0.])
VECY = np.asarray([0.,1.,0.])
VECZ = np.asarray([0.,0.,1.])


class NsaBase(object):
    """NsaBase"""
    def __init__(self, *arg):
        self.tag = arg[0]
        
class Material(NsaBase):
    """Material"""
    def __init__(self, *arg):
        super(Material,self).__init__(*arg)

class Material_Elastic(Material):
    """Material_Elastics"""
    def __init__(self, *arg):
        super(Material_Elastic,self).__init__(*arg)
        self.E0 = float(arg[1])
        self.v = float(arg[2]) if len(arg)>=3 else 0.3
        self.G0 = self.E0/2./(1.-self.v)
        
class GeometricTopology(NsaBase):
    """GeometricTopology"""
    def __init__(self, *arg):
        super(GeometricTopology,self).__init__(*arg)
    
class Node(GeometricTopology):
    """Node"""
    def __init__(self, *arg):

        super(Node, self).__init__(*arg)
        self.coord = np.asarray(arg[1],dtype=float)
        self.mass = float(arg[2]) if len(arg)>=3 else 0.0
        self.dim = 3
        self.dof = 3
        self.index = 1

    def get_global_dof_index(self):

        dof_index = np.arange(self.dof,dtype=int)
        dof_index_start = self.index*self.dof
        self.dof_index = dof_index+dof_index_start

    def get_deformation(self,U):

        self.deformation = np.asarray([U[self.dof_index[i]] for i in range(self.dof)])


class Section(GeometricTopology):
    """Section"""
    def __init__(self, *arg):
        super(Section, self).__init__(*arg)

class Section_Elastic(Section):
    """Section_Elastic"""
    def __init__(self, *arg):
        super(Section_Elastic, self).__init__(*arg)
        self.mat = arg[1]
        self.E = self.mat.E0
        self.G = self.mat.G0
        self.A = arg[2] if len(arg)>=3 else 0.0
        self.I = arg[3] if len(arg)>=4 else np.zeros(3)
        self.J = arg[4] if len(arg)>=5 else np.zeros(3)
        self.shear_factor = arg[5] if len(arg)>=6 else 1.0

        self.ka = self.E*self.A
        self.ks = self.G*self.A*self.shear_factor
        self.kf = self.E*self.I
        self.kt = self.E*self.J
        
class Element(GeometricTopology):
    """Element"""
    def __init__(self, *arg):

        super(Element, self).__init__(*arg)
        self.section = arg[1]
        self.node = arg[2]

class Line(Element):
    """Line"""
    def __init__(self, *arg):

        super(Line, self).__init__(*arg)
        self.vecx = self.node[1].coord - self.node[0].coord
        self.length = npl.norm(self.vecx)
        self.nvecx = self.vecx/self.length
        self.dof = self.node[0].dof+self.node[1].dof
        self.dim = self.node[0].dim
        self.dcosx = np.asarray([np.dot(self.nvecx,VECX[:self.dim]),np.dot(self.nvecx,VECY[:self.dim])])

    def assemble_stiffness_matrix(self):
    
        T0 = self.dcosx
        self.T = spl.block_diag(T0,T0)

        self.K = np.dot(np.dot(self.T.T,self.Kl),self.T)

        self.dof_index = np.hstack((self.node[0].dof_index
                            ,self.node[1].dof_index))

        self.rows = np.repeat(self.dof_index,self.dof)
        self.cols = np.tile(self.dof_index,self.dof)
        self.vals = self.K.flatten()
        self.rows = list(self.rows)
        self.cols = list(self.cols)
        self.vals = list(self.vals)

    def get_deformation_force(self):

        self.global_deformation = np.array([self.node[0].deformation,self.node[1].deformation]).flatten().T
        self.global_force = np.dot(self.K,self.global_deformation)
        self.local_deformation = np.dot(self.T,self.global_deformation)
        self.local_force = np.dot(self.Kl,self.local_deformation)

class Truss2D(Line):
    """Truss2D"""
    def __init__(self, *arg):
        super(Truss2D, self).__init__(*arg)
        self.k = self.section.ka/self.length
        self.Kl = self.k*np.array([[1.0,-1.0],[-1.0,1.0]])
        self.dim = 2
        self.dof = 4
        
class Truss3D(Truss2D):
    """Truss3D"""
    def __init__(self, *arg):
        super(Truss3D, self).__init__(*arg)
        self.dim = 3
        self.dof = 6

class Constraint(GeometricTopology):
    """Constraint"""
    def __init__(self, *arg):
        super(Constraint, self).__init__(*arg)
        self.node = arg[1]
        self.dof_tag = np.asarray(arg[2],dtype=int)-1
        self.dof_value = np.asarray(arg[3],dtype=float)
        
class Load(NsaBase):
    """Load"""
    def __init__(self, *arg):
        super(Load, self).__init__(*arg)
        self.node = arg[1]
        self.dof_tag = np.asarray(arg[2],dtype=int)-1
        self.dof_value = np.asarray(arg[3],dtype=float)

class Analysis(NsaBase):
    """Analysis"""
    def __init__(self, domain):
        self.domain = domain

class Analysis_Linear(Analysis):
    """Analysis_Linear"""
    def __init__(self, domain):
        super(Analysis_Linear, self).__init__(domain)

class Analysis_Linear_Eigen(Analysis_Linear):
    """Analysis_Linear"""
    def __init__(self, domain):
        super(Analysis_Linear_Eigen, self).__init__(domain)

    def execute(self):
        self.domain.Omega_2,self.domain.mode = spsl.eig(self.domain.K,self.domain.M)

class Analysis_Linear_Static(Analysis):
    """Analysis_Linear_Static"""
    def __init__(self, domain):
        super(Analysis_Linear_Static, self).__init__(domain)

    def execute(self):
        self.domain.U = spsl.spsolve(self.domain.K,self.domain.F)
        
class Domain(NsaBase):
    """Domain"""
    def __init__(self,*arg):

        self.DIM, self.DOF = int(arg[0]),int(arg[1])
        self.node = {}
        self.nodeindex = 0
        self.mat = {}
        self.sec = {}
        self.ele = {}
        self.cons = {}
        self.load = {}
        self.alength = 0.0

        self.K_rows = []
        self.K_cols = []
        self.K_vals = []
        self.M_vals = []

    def add_node(self,*arg):
        n = Node(*arg)
        n.dof = self.DOF
        n.dim = self.DIM
        n.index = self.nodeindex
        n.get_global_dof_index()
        n.coord = n.coord[:n.dim]
        self.M_vals += [n.mass for i in range(n.dof)]
        self.node[n.tag] = n
        self.nodeindex += 1


    def add_mat(self,mattype,*arg):
        mat = mattype(*arg)
        self.mat[mat.tag] = mat

    def add_sec(self,sectype,*arg):
        sec = sectype(*arg)
        self.sec[sec.tag] = sec

    def add_ele(self,eletype,*arg):
        ele = eletype(*arg)
        self.alength += ele.length
        ele.assemble_stiffness_matrix()
        self.K_rows += ele.rows
        self.K_cols += ele.cols
        self.K_vals += ele.vals
        self.ele[ele.tag] = ele
        
    def add_cons(self,*arg):
        cons = Constraint(*arg)
        self.cons[cons.tag] = cons

    def add_load(self,*arg):
        load = Load(*arg)
        self.load[load.tag] = load

    def build_model(self):

        self.nnode = len(self.node)
        self.nele = len(self.ele)
        self.alength /= self.nele
        self.ndof = self.nnode*self.DOF

        self.F = np.zeros(self.ndof)
        self.U = np.zeros(self.ndof)
        self.UF = np.zeros(self.ndof)

        self.assemble()
        self.apply_load()
        self.apply_cons()

    def assemble(self):
        self.K = sps.coo_matrix((self.K_vals, (self.K_rows,self.K_cols)), shape=(self.ndof,self.ndof))
        self.M = sps.diags(self.M_vals,0)

    def apply_load(self):
        for loadi in self.load.values():
            for i in range(len(loadi.dof_tag)):
                dof_index = loadi.node.dof_index[loadi.dof_tag[i]]
                self.F[dof_index] = loadi.dof_value[i]

    def apply_cons(self):
        K = self.K.toarray()
        for consi in self.cons.values():
            for i in range(len(consi.dof_tag)):
                dof_index = consi.node.dof_index[consi.dof_tag[i]]
                self.F = self.F-consi.dof_value[i]*K[:,dof_index]
                K[:,dof_index] = 0.0
                K[dof_index,:] = 0.0
                K[dof_index,dof_index] = 1.0
                self.F[dof_index] = consi.dof_value[i]
        self.K = sps.csr_matrix(K)

class PostProcessor(NsaBase):
    """PostProcessor"""
    def __init__(self, domain):
        self.domain = domain

    def get_all(self):
        self.get_node_deformation()
        self.get_element_deformation_force()
        self.get_unbalancedforce()

    def get_node_deformation(self):
        U = self.domain.U
        for ni in self.domain.node.values():
            ni.get_deformation(U)

    def get_element_deformation_force(self):
        for ei in self.domain.ele.values():
            ei.get_deformation_force()

    def get_unbalancedforce(self):
        for ei in self.domain.ele.values():
            for i in range(len(ei.dof_index)):
                self.domain.UF[ei.dof_index[i]] += -ei.global_force[i]

    def export_autocad_scr_2d(self,*fn):

        scrf = open('%s.scr'%fn,'w')
        scale = 0.1*self.domain.alength/max(abs(self.domain.U))
        for ni in self.domain.node.values():
            x,y = ni.coord
            print >>scrf,'POINT %.2f,%.2f,0.0'%(x,y)
            x,y = ni.coord + ni.deformation*scale
            print >>scrf,'POINT %.2f,%.2f,0.0'%(x,y)
        for ei in self.domain.ele.values():
            x1,y1 = ei.node[0].coord
            x2,y2 = ei.node[1].coord
            print >>scrf,'LINE %.2f,%.2f,0.0 %.2f,%.2f,0.0 '%(x1,y1,x2,y2)
            x1,y1 = ei.node[0].coord + ei.node[0].deformation*scale
            x2,y2 = ei.node[1].coord + ei.node[1].deformation*scale
            print >>scrf,'LINE %.2f,%.2f,0.0 %.2f,%.2f,0.0 '%(x1,y1,x2,y2)
        print >>scrf,'ZOOM E '
        scrf.close()
        # print 'AutoCAD Script Exported. '



def main():

    truss = Domain(2,2)

    truss.add_mat(Material_Elastic,1,1.)
    truss.add_sec(Section_Elastic,1,truss.mat[1],1.)

    d = 1000.0

    for i in range(7):
        truss.add_node(11+i,(i*d,0,0),1.)
        truss.add_node(21+i,(i*d,d,0),1.)
        truss.add_ele(Truss2D,'V-%d'%(1+i),truss.sec[1],(truss.node[11+i],truss.node[21+i]))
        truss.add_load(1+i,truss.node[21+i],[2],[-1.])

    truss.load[1].dof_value[0] = -.5
    truss.load[7].dof_value[0] = -.5

    truss.add_cons(1,truss.node[11],[1,2],[0.0,0.0])
    truss.add_cons(2,truss.node[17],[1,2],[0.0,0.0])

    for i in range(6):
        truss.add_ele(Truss2D,'H1-%d'%(1+i),truss.sec[1],(truss.node[11+i],truss.node[12+i]))
        truss.add_ele(Truss2D,'H2-%d'%(1+i),truss.sec[1],(truss.node[21+i],truss.node[22+i]))
        if i < 3:
            truss.add_ele(Truss2D,'DL-%d'%(1+i),truss.sec[1],(truss.node[12+i],truss.node[21+i]))
        else:
            truss.add_ele(Truss2D,'DR-%d'%(1+i),truss.sec[1],(truss.node[11+i],truss.node[22+i]))     

    truss.build_model()

    analyse = Analysis_Linear_Static(truss)
    analyse.execute()

    post = PostProcessor(truss)
    post.get_all()
    post.export_autocad_scr_2d('truss')

    return truss



if __name__ == '__main__':

    # main = autojit(main)
    truss = main()