# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis
**用于结构非线性分析的Python模块

'''

import numpy as np
import scipy.sparse as sps
import scipy.linalg as spl
from scipy.special import p_roots

from nsapy.Nsabase import Element,RefPoint
from nsapy.Nsabase import VECX,VECY,VECZ

from copy import deepcopy

def dcos(vec):
    '''计算向量的方向余弦'''

    dim = len(vec)
    res = np.zeros(dim)
    res[0] = vec.dot(VECX[:dim])
    res[1] = vec.dot(VECY[:dim])
    if dim>2:
        res[2] = vec.dot(VECZ[:dim])
    return res

class Line(Element):
    """Line"""
    def __init__(self, *arg):

        # 初始化单元
        super(Line, self).__init__(*arg)
        self.vecx = self.node[1].coord - self.node[0].coord # 单元局部x轴向量
        self.length = np.linalg.norm(self.vecx) # 单元长度
        self.nvecx = self.vecx/self.length # 单元局部x轴单位向量
        self.dof = self.node[0].dof+self.node[1].dof # 单元节点自由度总数
        self.dim = self.node[0].dim # 单元维度
        self.local_dof = 1 # 单元局部自由度
        # 局部x轴的方向余弦
        self.dcosx = dcos(self.nvecx)
        self.local_force = np.zeros(self.local_dof)
        
    def get_trans_matrix(self):
        '''几何转换矩阵'''
        pass

    def get_stiffness_matrix(self):
        self.K = self.T.T.dot(self.Ke).dot(self.T)

    def assemble_stiffness_matrix(self):
        
        '''将单元刚度矩阵定位到全局刚度矩阵'''
        # 读取节点自由度的全局编号
        self.dof_index = np.hstack((self.node[0].dof_index
                            ,self.node[1].dof_index))

        #　将局部刚度矩阵按行重新排列成一维向量
        self.vals = self.K.flatten() # 向量化的刚度矩阵
        self.rows = np.repeat(self.dof_index,self.dof) # 刚度元素在全局矩阵中的行号
        self.cols = np.tile(self.dof_index,self.dof) # 刚度元素在全局矩阵中的列号
        
        # 将上述向量转化为list对象以便于组合操作
        self.rows = list(self.rows)
        self.cols = list(self.cols)
        self.vals = list(self.vals)

    def get_deformation_force(self):
        ''' 获取节点变形及单元内力 '''
        self.get_deformation()
        self.get_force()

    def get_deformation(self):
        ''' 获取节点变形及应变 '''
        self.global_deformation = np.array([self.node[0].deformation,self.node[1].deformation]).flatten().T
        self.local_deformation = self.T.dot(self.global_deformation)
        self.strain = self.local_deformation/self.length

    def get_force(self):
        ''' 获取单元内力 '''
        self.section.set_strain(self.strain)
        self.section.get_force_stiffness()
        self.local_force = self.section.inner_force
        self.global_force = self.T.T.dot(self.local_force).flatten()
        
    def update_stiffness(self):

        self.get_stiffness()
        self.get_stiffness_matrix()
        self.assemble_stiffness_matrix()

    def update_history_para(self):
        '''非线性迭代收敛时将状态变量更新'''
        self.section.update_history_para()
        
            
class Truss2D(Line):
    """Truss2D"""
    def __init__(self, *arg):
        super(Truss2D, self).__init__(*arg)
        
        self.dim = 2
        self.dof = 4
        self.get_stiffness()
        self.get_trans_matrix()
        self.get_stiffness_matrix()
        
    def get_stiffness(self):
        self.Ke = np.array([[self.section.Da/self.length]])
        
    def get_trans_matrix(self):
        '''几何转换矩阵'''
        self.nvecx = self.nvecx[:self.dim]
        self.dcosx = self.dcosx[:self.dim]
        
        self.nvecy = np.array([-self.nvecx[1],self.nvecx[0]])
        self.dcosy = dcos(self.nvecy)
        self.dcos = np.vstack((self.dcosx,self.dcosy))

        self.T0 = np.array([[1.],[0.]]).T
        self.T1 = np.array([[-1.,0.,1.,0.],[0.,-1.,0.,1.]])
        self.T2 = spl.block_diag(self.dcos,self.dcos)
        self.T = self.T0.dot(self.T1).dot(self.T2)
        
class Truss3D(Line):
    """Truss3D"""
    def __init__(self, *arg):
        super(Truss3D, self).__init__(*arg)

        self.dim = 3
        self.dof = 6
        
        self.ref_point = arg[3] if len(arg)>3 else None    
        
        self.get_stiffness()
        self.get_trans_matrix()
        self.get_stiffness_matrix()

    def get_stiffness(self):
        self.Ke = np.array([[self.section.Da/self.length]])

    def get_stiffness_matrix(self):
        self.K = self.T.T.dot(self.Ke).dot(self.T)
        
    def get_trans_matrix(self):
        '''几何转换矩阵'''
        
        vecr = self.ref_point.coord - self.node[0].coord if self.ref_point is not None else VECY
        vecz = np.cross(self.nvecx,vecr)
        self.nvecz = vecz/np.linalg.norm(vecz)
        vecy = np.cross(self.nvecz,self.nvecx)
        self.nvecy = vecy/np.linalg.norm(vecy)
        
        self.dcosy = dcos(self.nvecy)
        self.dcosz = dcos(self.nvecz)
        self.dcos = np.vstack((self.dcosx,self.dcosy,self.dcosz))

        self.T0 = np.array([[1.],[0.],[0.]]).T
        self.T1 = np.array([[-1.,0.,0.,1.,0.,0.],[0.,-1.,0.,0.,1.,0.],[0.,0.,-1.,0.,0.,1.]])
        self.T2 = spl.block_diag(self.dcos,self.dcos)
        self.T = self.T0.dot(self.T1).dot(self.T2)

class Frame2D_Elastic(Line):
    """Frame2D"""
    def __init__(self, *arg):
        super(Frame2D_Elastic, self).__init__(*arg)
        self.dim = 2
        self.dof = 6

        self.local_dof = 3 # 单元局部自由度
        self.strain = np.zeros(self.local_dof)

        self.get_stiffness()
        self.get_trans_matrix()
        self.get_stiffness_matrix()
        
    def get_stiffness(self):

        Kea = [self.section.Da]
        Keb = self.section.Db/self.length*np.array([[4.,2.],[2.,4.]])
        self.Ke = spl.block_diag(Kea,Keb)
        
    def get_trans_matrix(self):
        '''几何转换矩阵'''
        self.nvecx = self.nvecx[:self.dim]
        self.dcosx = self.dcosx[:self.dim]
        
        self.nvecy = np.array([-self.nvecx[1],self.nvecx[0]])
        self.dcosy = dcos(self.nvecy)
        self.dcos = np.vstack((self.dcosx,self.dcosy))  

        oL = 1./self.length
        self.T0 = np.array([
                            [1.0,0.0,0.0,0.0],
                            [0.0,-oL,1.0,0.0],
                            [0.0,-oL,0.0,1.0],
                           ])
        self.T1 = np.array([
                            [-1.0,0.0,0.0,1.0,0.0,0.0],
                            [0.0,-1.0,0.0,0.0,1.0,0.0],
                            [ 0.0,0.0,1.0,0.0,0.0,0.0],
                            [ 0.0,0.0,0.0,0.0,0.0,1.0],
                           ])
        self.T2 = spl.block_diag(self.dcos,[1.0],self.dcos,[1.0])
        self.T = self.T0.dot(self.T1).dot(self.T2)

    def get_deformation(self):
        ''' 获取节点变形及应变 '''
        self.global_deformation = np.array([self.node[0].deformation,self.node[1].deformation]).flatten().T
        self.local_deformation = self.T.dot(self.global_deformation)
        self.strain[0] = self.local_deformation[0]/self.length
        self.strain[1] = self.local_deformation[1]*-4./self.length + self.local_deformation[2]*-2./self.length
        self.strain[2] = self.local_deformation[1]*2./self.length + self.local_deformation[2]*4./self.length

class Frame2D(Line):
    """Nonlinear Frame2D"""
    def __init__(self, *arg):
        super(Frame2D, self).__init__(*arg)
        self.n_intps = arg[3] if len(arg)>3 else 5
        self.dim = 2
        self.dof = 6

        self.local_dof = 3 # 单元局部自由度
        self.strain = np.zeros(self.local_dof)

        x,w = p_roots(self.n_intps)
        x = x.real
        self.loc_intps,self.wf_intps = 0.5*(x+1.0),0.5*w
        
        self.x = self.length*self.loc_intps
        self.intps = []
        self.Dax = np.ones(self.n_intps)
        self.Dbx = np.ones(self.n_intps)
        for i in xrange(self.n_intps):
            self.intps.append(deepcopy(self.section))
            self.Dax[i] = self.intps[i].Da
            self.Dbx[i] = self.intps[i].Db

        self.kappax = np.zeros(self.n_intps)
        
        self.get_stiffness()
        self.get_trans_matrix()
        self.get_stiffness_matrix()
             
    def get_stiffness(self):

        Kea = np.ones((1,1))
        Keb = np.ones((2,2))

        for i in xrange(self.n_intps):
            self.Dax[i] = self.intps[i].Da
            self.Dbx[i] = self.intps[i].Db
        
        # 采用高斯勒让德积分计算单元刚度
        ooL = 1./self.length
        Kea[0,0] = ooL*self.wf_intps.dot(self.Dax)
        Keb[0,0] = ooL*self.wf_intps.dot(self.Dbx*(6.*self.loc_intps-4.)**2)
        Keb[1,1] = ooL*self.wf_intps.dot(self.Dbx*(6.*self.loc_intps-2.)**2)
        Keb[0,1] = ooL*self.wf_intps.dot(self.Dbx*(6.*self.loc_intps-4.)*(6.*self.loc_intps-2.))
        Keb[1,0] = Keb[0,1]*1.0

        self.Ke = spl.block_diag(Kea,Keb)

    def get_trans_matrix(self):
        '''几何转换矩阵'''
        self.nvecx = self.nvecx[:self.dim]
        self.dcosx = self.dcosx[:self.dim]
        
        self.nvecy = np.array([-self.nvecx[1],self.nvecx[0]])
        self.dcosy = dcos(self.nvecy)
        self.dcos = np.vstack((self.dcosx,self.dcosy))  

        oL = 1./self.length
        self.T0 = np.array([
                            [1.0,0.0,0.0,0.0],
                            [0.0,-oL,1.0,0.0],
                            [0.0,-oL,0.0,1.0],
                           ])
        self.T1 = np.array([
                            [-1.0,0.0,0.0,1.0,0.0,0.0],
                            [0.0,-1.0,0.0,0.0,1.0,0.0],
                            [ 0.0,0.0,1.0,0.0,0.0,0.0],
                            [ 0.0,0.0,0.0,0.0,0.0,1.0],
                           ])
        self.T2 = spl.block_diag(self.dcos,[1.0],self.dcos,[1.0])
        self.T = self.T0.dot(self.T1).dot(self.T2)
    
    def get_deformation(self):
        ''' 获取节点变形及应变 '''
        self.global_deformation = np.array([self.node[0].deformation,self.node[1].deformation]).flatten().T
        self.local_deformation = self.T.dot(self.global_deformation)
        self.strain[0] = self.local_deformation[0]/self.length
        self.strain[1] = self.local_deformation[1]*-4./self.length + self.local_deformation[2]*-2./self.length
        self.strain[2] = self.local_deformation[1]*2./self.length + self.local_deformation[2]*4./self.length

        oL2 = 1./self.length**2
        self.kappax = -oL2*((4.*self.length-6.*self.x)*self.local_deformation[1]+(2.*self.length-6.*self.x)*self.local_deformation[2])            
    
    def get_force(self):
        ''' 获取单元内力 '''
        self.local_force = np.zeros(3)
        for i in xrange(self.n_intps):
            sec = self.intps[i]
            sec.set_strain([self.strain[0],self.kappax[i]])
            sec.get_force_stiffness()
            self.Dbx[i] = self.intps[i].Db
            self.local_force[0]+=self.wf_intps[i]*sec.inner_force[0]
            self.local_force[1]+=self.wf_intps[i]*(self.loc_intps[i]*6.0-4.0)*sec.inner_force[1]
            self.local_force[2]+=self.wf_intps[i]*(self.loc_intps[i]*6.0-2.0)*sec.inner_force[1]
        self.global_force = self.T.T.dot(self.local_force).flatten()

    def update_history_para(self):
        '''非线性迭代收敛时将状态变量更新'''
        for i in xrange(self.n_intps):
            sec = self.intps[i]
            sec.update_history_para()

class Spring2D(Truss2D):
    """docstring for Spring2D"""
    def __init__(self, *arg):
        self.node = arg[1]
        self.mat1 = deepcopy(arg[2])
        self.mat2 = deepcopy(arg[3])

    def get_trans_matrix(self):
        '''几何转换矩阵'''
        self.nvecx = self.nvecx[:self.dim]
        self.dcosx = self.dcosx[:self.dim]
        
        self.nvecy = np.array([-self.nvecx[1],self.nvecx[0]])
        self.dcosy = dcos(self.nvecy)
        self.dcos = np.vstack((self.dcosx,self.dcosy))

        self.T0 = np.diag([1.,1.])
        self.T1 = np.array([[-1.,0.,1.,0.],[0.,-1.,0.,1.]])
        self.T2 = spl.block_diag(self.dcos,self.dcos)
        self.T = self.T0.dot(self.T1).dot(self.T2)

    def get_stiffness(self):
        self.Ke = np.diag([self.mat1.tangent,self.mat2.tangent])

    def get_deformation(self):
        ''' 获取节点变形 '''
        self.global_deformation = np.array([self.node[0].deformation,self.node[1].deformation]).flatten().T
        self.local_deformation = self.T.dot(self.global_deformation)

    def get_force(self):
        ''' 获取单元内力 '''
        self.mat1.set_strain(self.local_deformation[0])
        self.mat1.get_stress_tangent()
        self.mat2.set_strain(self.local_deformation[1])
        self.mat2.get_stress_tangent()
        self.inner_force = np.array([[self.mat1.stress,self.mat2.stress]])

        self.local_force = self.section.inner_force
        self.global_force = self.T.T.dot(self.local_force.T).flatten()

    def update(self):
        self.mat1.update_history_para()
        self.mat2.update_history_para()