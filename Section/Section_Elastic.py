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