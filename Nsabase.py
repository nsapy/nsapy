# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis
**用于结构非线性分析的Python模块

'''
import numpy as np
import scipy.sparse as sps
from copy import deepcopy

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
    
class RefPoint(GeometricTopology):
    """RefPoint"""
    def __init__(self, *arg):

        super(RefPoint, self).__init__(*arg)
        self.coord = np.array(arg[1],dtype=float)

class Node(GeometricTopology):
    """Node"""
    def __init__(self, *arg):

        super(Node, self).__init__(*arg)
        self.coord = np.array(arg[1],dtype=float)
        self.mass = np.array(arg[2],dtype=float) if len(arg)>=3 else None
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

    def get_unbalanced_force(self,UF):

        self.UF = np.array([UF[self.dof_index[i]] for i in range(self.dof)])
        self.reaction = -self.UF

class Section(GeometricTopology):
    """Section"""
    def __init__(self, *arg):
        super(Section, self).__init__(*arg)
        self.mat = deepcopy(arg[1])
        
class Element(GeometricTopology):
    """Element"""
    def __init__(self, *arg):

        super(Element, self).__init__(*arg)
        self.section = deepcopy(arg[1])
        self.node = arg[2]

class Constraint(NsaBase):
    """Constraint"""
    def __init__(self, *arg):
        super(Constraint, self).__init__(*arg)
        self.node = arg[1]
        self.dof_tag = np.array(arg[2],dtype=int)-1
        self.dof_value = np.array(arg[3],dtype=float)
        self.dof_index = []
        
        
class Load(NsaBase):
    """Load"""
    def __init__(self, *arg):
        super(Load, self).__init__(*arg)
        self.node = arg[1]
        self.dof_tag = np.array(arg[2],dtype=int)-1
        self.dof_value = np.array(arg[3],dtype=float)
        self.dof_value_der = np.array(arg[3],dtype=float)
        self.dof_value_der[:,1:] = np.diff(self.dof_value_der)
        self.nsteps = len(self.dof_value[0])
        self.time = np.array(arg[4],dtype=float) if len(arg)>4 else np.arange(self.nsteps)
        
class Analysis(NsaBase):
    """Analysis"""
    def __init__(self, domain):
        self.domain = domain
        self.create_data_file()

    def create_data_file(self):
        meshfile = open(self.domain.name+'.msh','r')
        meshdata = meshfile.read()
        meshfile.close()

        self.nodedatafile = open(self.domain.name+'_nodedata.msh','w')
        self.nodedatafile.write(meshdata)
        
        self.eledatafile = open(self.domain.name+'_eledata.msh','w')
        self.eledatafile.write(meshdata)

    def save_data_file(self):
        self.nodedatafile.close()
        self.eledatafile.close()

    def write_data_file(self):

        self.nodedatafile.write('$NodeData\n')
        self.nodedatafile.write('1\n\"%s\"\n'%(self.domain.name))
        self.nodedatafile.write('1\n%f\n'%self.ctime)
        self.nodedatafile.write('3\n%d\n%d\n%d\n'%(0,3,self.domain.nnode))
        for ni in self.domain.node.values():
            u = np.zeros(3)
            if ni.dof in [2,3]:
                dof = 2
            elif ni.dof in [4,6]:
                dof = 3
            u[:dof] = ni.deformation[:dof]
            u1,u2,u3 = tuple(u)
            self.nodedatafile.write('%d %f %f %f\n'%(ni.tag,u1,u2,u3))

        self.nodedatafile.write('$EndNodeData\n')

        self.eledatafile.write('$ElementData\n')
        self.eledatafile.write('1\n\"Load Step %d\"\n'%self.cstep)
        self.eledatafile.write('1\n%f\n'%self.ctime)
        self.eledatafile.write('3\n%d\n%d\n%d\n'%(self.cstep,3,self.domain.nele))

        for ei in self.domain.ele.values():
            f = np.zeros(3)
            f[:ei.local_dof] = ei.local_force
            f1,f2,f3 = tuple(f)
            self.eledatafile.write('%d %f %f %f\n'%(ei.tag,f1,f2,f3))
                    
        self.eledatafile.write('$EndElementData\n')

    def execute(self,nsteps):
        pass
        
class Domain(NsaBase):
    """Domain"""
    def __init__(self,*arg):

        self.DIM, self.DOF = int(arg[0]),int(arg[1])
        self.name = arg[2] if len(arg)>2 else 'domain'
        self.refp = {}
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

        self.constraint_dof = {}
        self.constraint_value = []
        self.load_dof = []
        self.load_value = []
        
    def add_refpoint(self,*arg):
        rp = RefPoint(*arg)
        self.refp[rp.tag] = rp

    def add_node(self,*arg):
        n = Node(*arg)
        n.dof = self.DOF
        n.dim = self.DIM
        n.index = self.nodeindex
        n.get_global_dof_index()
        n.coord = n.coord[:n.dim]
        self.M_vals += [n.mass[i] for i in range(n.dof)]
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
        for i in range(len(cons.dof_tag)):
            dof_index = cons.node.dof_index[cons.dof_tag[i]]
            cons.dof_index.append(dof_index)
            self.constraint_dof[dof_index] = cons.dof_value[i]
        self.cons[cons.tag] = cons

    def add_load(self,*arg):
        load = Load(*arg)
        for i in range(len(load.dof_tag)):
            dof_index = load.node.dof_index[load.dof_tag[i]]
            self.load_dof.append(dof_index)
        self.load[load.tag] = load

    def build_model(self):

        self.nnode = len(self.node)
        self.nele = len(self.ele)
        self.alength /= self.nele
        self.ndof = self.nnode*self.DOF

        self.dof_index_list = range(self.ndof)
        self.free_dof = list(set(self.dof_index_list) - set(self.constraint_dof.keys()))
        
        self.constraint_ndof = len(self.constraint_dof)
        self.free_ndof = self.ndof - self.constraint_ndof

        self.F = np.zeros(self.ndof)
        self.dF = np.zeros(self.ndof)
        self.U = np.zeros(self.ndof)
        self.dU = np.zeros(self.ndof)
        self.UF = np.zeros(self.ndof)
        self.Ft = np.zeros(self.ndof)

        # self.Reduced_F  = np.zeros(self.free_ndof)
        # self.Reduced_dF = np.zeros(self.free_ndof)
        # self.Reduced_U  = np.zeros(self.free_ndof)
        # self.Reduced_dU = np.zeros(self.free_ndof)
        # self.Reduced_UF = np.zeros(self.free_ndof)

        self.assemble()
        self.write_gmsh()

    def assemble(self):

        self.K = sps.coo_matrix((self.K_vals, (self.K_rows,self.K_cols)), shape=(self.ndof,self.ndof))
        self.M = sps.diags(self.M_vals,0)

    def reassemble_stiffness_matrix(self):
        self.K_rows = []
        self.K_cols = []
        self.K_vals = []
        
        for ei in self.ele.values():
            ei.update_stiffness()
            self.K_rows += ei.rows
            self.K_cols += ei.cols
            self.K_vals += ei.vals

        self.K = sps.coo_matrix((self.K_vals, (self.K_rows,self.K_cols)), shape=(self.ndof,self.ndof))

    def apply_load(self,step):
        '''施加节点荷载'''
        for loadi in self.load.values():
            for i in range(len(loadi.dof_tag)):
                dof_index = loadi.node.dof_index[loadi.dof_tag[i]]
                self.F[dof_index] = loadi.dof_value[i][step]
                self.dF[dof_index] = loadi.dof_value_der[i][step]

    def apply_cons(self,step):
        '''施加支座约束'''
        K = self.K.toarray()

        for dof_index in self.constraint_dof.keys():
            dof_value = self.constraint_dof[dof_index]
            self.F = self.F-dof_value*K[:,dof_index]
            if step==0:
                self.dF = self.dF-dof_value*K[:,dof_index]

            K[:,dof_index] = 0.0
            K[dof_index,:] = 0.0
            K[dof_index,dof_index] = 1.0

        self.K = sps.csc_matrix(K)

    def get_result(self):

        '''为计算节点反力向量'''

        for ni in self.node.values():
            ni.get_deformation(self.U)
        
        self.Ft = np.zeros(self.ndof)

        for ei in self.ele.values():
            ei.get_deformation_force()
            for i in range(len(ei.dof_index)):
                self.Ft[ei.dof_index[i]] += ei.global_force[i]
        
        self.UF = -self.Ft

        # 约束自由度反力置零
        for i in self.constraint_dof.keys():
            self.Ft[i] = 0.0

        for ni in self.node.values():
            ni.get_unbalanced_force(self.UF)

    def update(self):
        
        for ni in self.node.values():
            ni.get_deformation(self.U)

        self.K_rows = []
        self.K_cols = []
        self.K_vals = []
        self.UF = np.zeros(self.ndof)
        
        for ei in self.ele.values():
            ei.get_deformation_force()
            # 更新应力状态
            ei.update()
            # 更新单元刚度矩阵
            ei.update_stiffness()
            self.K_rows += ei.rows
            self.K_cols += ei.cols
            self.K_vals += ei.vals
            # 更新节点不平衡力        
            for i in range(len(ei.dof_index)):
                self.UF[ei.dof_index[i]] += -ei.global_force[i]

        # 更新整体刚度矩阵
        self.K = sps.coo_matrix((self.K_vals, (self.K_rows,self.K_cols)), shape=(self.ndof,self.ndof))

        for ni in self.node.values():
            ni.get_unbalanced_force(self.UF)

    def write_gmsh(self):
        gmshfile = open('%s.msh'%self.name,'w')
        gmsh_head = '''$MeshFormat
2.8 0 8
$EndMeshFormat\n'''
        gmshfile.write(gmsh_head)
        gmshfile.write('$Nodes\n')
        gmshfile.write('%d\n'%self.nnode)
        for ni in self.node.values():
            i = ni.tag
            coord = np.zeros(3)
            coord[:ni.dim] = ni.coord
            x,y,z = tuple(coord)
            print >>gmshfile,'%d %f %f %f'%(i,x,y,z)
        gmshfile.write('$EndNodes\n')
        
        gmshfile.write('$Elements\n')
        gmshfile.write('%d\n'%self.nele)
        for ei in self.ele.values():
            i = ei.tag
            n1 = ei.node[0].tag
            n2 = ei.node[1].tag
            print >>gmshfile,'%d 1 1 1 %d %d'%(i,n1,n2)
        gmshfile.write('$EndElements\n')
        gmshfile.close()

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
    
    
    def export_gmsh_scr_2d(self):

        scrf = open('%s.geo'%self.name,'w')
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
        pass
    
class NodeResult(Results):
    """Results"""
    def __init__(self,domain,*tag):
        super(NodeResult, self).__init__()
        self.domain = domain
        self.node = self.domain.node[tag[0]]
        self.fieldtag = tag[1]
        self.doftag = tag[2]-1

        self.values = []

    def add_node(self,tag,fieldtag):
        pass

    def __add__(self,other):
        self.values = np.array(self.values)
        other.values = np.array(other.values)
        self.values += other.values
        return self

    def save_result(self):
        exec 'self.values.append(self.node.%s[%d])'%(self.fieldtag,self.doftag)
        

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

    def get_unbalanced_force(self):
        for ei in self.domain.ele.values():
            for i in range(len(ei.dof_index)):
                self.domain.UF[ei.dof_index[i]] += -ei.global_force[i]

if __name__ == '__main__':

    pass