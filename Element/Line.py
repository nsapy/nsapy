# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis
**用于结构非线性分析的Python模块

'''

import numpy as np
import scipy.sparse as sps
import scipy.linalg as spl

from nsapy.Nsabase import Element
from nsapy.Nsabase import VECX,VECY,VECZ

class Line(Element):
    """Line"""
    def __init__(self, *arg):

        super(Line, self).__init__(*arg)
        self.vecx = self.node[1].coord - self.node[0].coord
        self.length = np.linalg.norm(self.vecx)
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
