# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis
**用于结构非线性分析的Python模块

'''
import numpy as np
import scipy.sparse as sps

VECX = np.array([1.,0.,0.])
VECY = np.array([0.,1.,0.])
VECZ = np.array([0.,0.,1.])

class NsaBase(object):
    """NsaBase"""
    def __init__(self, *arg):
        self.tag = arg[0]
        
class Material(NsaBase):
    """Material"""
    def __init__(self, *arg):
        super(Material,self).__init__(*arg)

class Material1D(Material):
    """Material"""
    def __init__(self, *arg):
        super(Material1D,self).__init__(*arg)
        
class GeometricTopology(NsaBase):
    """GeometricTopology"""
    def __init__(self, *arg):
        super(GeometricTopology,self).__init__(*arg)
    
class Node(GeometricTopology):
    """Node"""
    def __init__(self, *arg):

        super(Node, self).__init__(*arg)
        self.coord = np.array(arg[1],dtype=float)
        self.mass = float(arg[2]) if len(arg)>=3 else 0.0
        self.dim = 3
        self.dof = 3
        self.index = 1
        self.coord0 = self.coord

    def get_global_dof_index(self):

        dof_index = np.arange(self.dof,dtype=int)
        dof_index_start = self.index*self.dof
        self.dof_index = dof_index+dof_index_start

    def get_deformation(self,U):

        self.deformation = np.array([U[self.dof_index[i]] for i in range(self.dof)])

    def update_coord(self):

        self.coord += self.deformation


class Section(GeometricTopology):
    """Section"""
    def __init__(self, *arg):
        super(Section, self).__init__(*arg)
        
class Element(GeometricTopology):
    """Element"""
    def __init__(self, *arg):

        super(Element, self).__init__(*arg)
        self.section = arg[1]
        self.node = arg[2]
        self.index = 1

class Constraint(NsaBase):
    """Constraint"""
    def __init__(self, *arg):
        super(Constraint, self).__init__(*arg)
        self.node = arg[1]
        self.dof_tag = np.array(arg[2],dtype=int)-1
        self.dof_value = np.array(arg[3],dtype=float)
        
        
class Load(NsaBase):
    """Load"""
    def __init__(self, *arg):
        super(Load, self).__init__(*arg)
        self.node = arg[1]
        self.dof_tag = np.array(arg[2],dtype=int)-1
        self.dof_value = np.array(arg[3],dtype=float)
        self.nsteps = len(self.dof_value[0])
        self.time = np.array(arg[4],dtype=float) if len(arg)>5 else np.arange(1,self.nsteps)
        
class Analysis(NsaBase):
    """Analysis"""
    def __init__(self, domain):
        self.domain = domain
        self.nsteps = 1
        
        
class Domain(NsaBase):
    """Domain"""
    def __init__(self,*arg):

        self.DIM, self.DOF = int(arg[0]),int(arg[1])
        self.node = {}
        self.nodeindex = 0
        self.mat = {}
        self.sec = {}
        self.ele = {}
        self.eleindex = 0
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
        ele.index = self.eleindex
        self.alength += ele.length
        ele.assemble_stiffness_matrix()
        self.K_rows += ele.rows
        self.K_cols += ele.cols
        self.K_vals += ele.vals
        self.ele[ele.tag] = ele
        self.eleindex += 1
        
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

    def assemble(self):
        self.K = sps.coo_matrix((self.K_vals, (self.K_rows,self.K_cols)), shape=(self.ndof,self.ndof))
        self.M = sps.diags(self.M_vals,0)

    def apply_load(self,step):
        '''施加节点荷载'''
        for loadi in self.load.values():
            for i in range(len(loadi.dof_tag)):
                dof_index = loadi.node.dof_index[loadi.dof_tag[i]]
                self.F[dof_index] = loadi.dof_value[i][step]-loadi.dof_value[i][step-1] if step>=1 else loadi.dof_value[i][step]

    def apply_cons(self):
        '''施加支座约束'''
        K = self.K.toarray()
        for consi in self.cons.values():
            for i in range(len(consi.dof_tag)):
                dof_index = consi.node.dof_index[consi.dof_tag[i]]
                self.F = self.F-consi.dof_value[i]*K[:,dof_index]
                K[:,dof_index] = 0.0
                K[dof_index,:] = 0.0
                K[dof_index,dof_index] = 1.0
                self.F[dof_index] = consi.dof_value[i]
        self.K = sps.csc_matrix(K)

    def read_gmsh(self, fn, *para):
        gmshfile = open('%s.msh'%fn,'r')
        nodeon,eleon = 0,0
        eletype = para[0]
        line = 'read_gmsh'
        while line!='':
            if line == '$Nodes': 
                nodeon = 1
                nnode = int(gmshfile.readline())
            elif line == '$EndNodes': nodeon = 0
            if line == '$Elements': 
                eleon = 1
                nele = int(gmshfile.readline())
            elif line == '$EndElements': eleon = 0
            if nodeon == 1:
                line = gmshfile.readline().strip('\n')
                line = line.replace(' ',',')
                exec 'tag,x,y,z = %s'%line
                self.add_node(tag,[x,y,z])
            elif eleon == 1:
                line = gmshfile.readline().strip('\n')
                line = line.replace(' ',',')
                exec 'tag,x,y,z = %s'%line
                self.add_ele(tag,eletype,[x,y,z])
            else:
                line = gmshfile.readline().strip('\n')
    
    
    def export_gmsh_scr_2d(self,*fn):

        scrf = open('%s.geo'%fn,'w')
        ps = 0.02*self.alength
        for ni in self.node.values():
            x,y = ni.coord
            i = ni.index+1
            print >>scrf,'Point(%d) = {%.2f,%.2f,0.0,%.2e};'%(i,x,y,ps)
        
        for ei in self.ele.values():
            n1 = ei.node[0].index+1
            n2 = ei.node[1].index+1
            print >>scrf,'Line(%d) = {%d,%d};'%((ei.index+1),n1,n2)
            i+=1
        scrf.close()
        print 'Gmsh model file \"%s.geo\" has been exported.'%fn

class Results(NsaBase):
    """Results"""
    def __init__(self):
        super(Results, self).__init__(*arg)
    
class Results_Node(Results):
    """Results"""
    def __init__(self,domain):
        super(Results_Node, self).__init__(*arg)
        self.domain = domain
        self.node = {}

    def add_node(self,tag):
        self.node[tag] = self.domain.node[tag]

    def save_result(self):
        pass




        

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

    def export_gmsh_scr_2d(self,*fn):

        scrf = open('%s.geo'%fn,'w')
        ps = 0.01*self.domain.alength
        for ni in self.domain.node.values():
            x,y = ni.coord
            i = ni.index+1
            print >>scrf,'Point(%d) = {%.2f,%.2f,0.0,%.2e};'%(i,x,y,ps)
        
        i = 1
        for ei in self.domain.ele.values():
            n1 = ei.node[0].index+1
            n2 = ei.node[1].index+1
            print >>scrf,'Line(%d) = {%d,%d};'%(i,n1,n2)
            i+=1
        scrf.close()
        
        print 'Gmsh model file \"%s\" has been exported.'%fn

if __name__ == '__main__':

    pass