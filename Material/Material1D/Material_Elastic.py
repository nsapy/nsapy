# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis
**用于结构非线性分析的Python模块

'''

import numpy as np
import scipy.sparse as sps

from nsapy.Nsabase import Material1D

class Material_Elastic(Material1D):
    """Material_Elastic"""
    def __init__(self, *arg):
        super(Material_Elastic,self).__init__(*arg)
        self.E0 = float(arg[1])
        # self.v = float(arg[2]) if len(arg)>=3 else 0.3
        # self.G0 = self.E0/2./(1.-self.v)
        self.strain = 0.0
        self.stress = 0.0
        self.tangent = self.E0

    def set_strain(self,strain):
        self.strain = strain

    def get_stress_tangent(self):
        self.stress = self.strain*self.E0
        self.tangent = self.E0

