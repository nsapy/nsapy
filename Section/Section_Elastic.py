# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis
**用于结构非线性分析的Python模块

'''

import numpy as np
import scipy.sparse as sps

from nsapy.Nsabase import Section

class Section_Elastic(Section):
    """Section_Elastic"""
    def __init__(self, *arg):
        super(Section_Elastic, self).__init__(*arg)
        self.E = self.mat.E0
        # self.G = self.mat.G0
        self.A = arg[2] if len(arg)>=3 else 0.0
        # self.I = arg[3] if len(arg)>=4 else np.zeros(2)
        # self.J = arg[4] if len(arg)>=5 else 0.0
        # self.shear_factor = arg[5] if len(arg)>=6 else 1.0

        self.Da = self.E*self.A
        # self.Ds = self.G*self.A*self.shear_factor
        # self.Db = self.E*self.I
        # self.Dt = self.E*self.J

        self.inner_force = np.zeros(5)

class Elastic_Truss(Section_Elastic):
    """Section_Elastic"""
    def __init__(self, *arg):
        super(Elastic_Truss, self).__init__(*arg)

    def get_force_stiffness(self,strain):
        self.inner_force = strain*self.Da

class Elastic_Frame2D(Section_Elastic):
    """Section_Elastic"""
    def __init__(self, *arg):
        super(Elastic_Frame2D, self).__init__(*arg)
        self.I = float(arg[3]) if len(arg)>=4 else 0.0
        self.Db = self.E*self.I
        self.D = np.array([self.Da,self.Db,self.Db])

    def get_force_stiffness(self,strain):

        self.inner_force = strain*self.D






