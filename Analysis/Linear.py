# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis
**用于结构非线性分析的Python模块

'''

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl

from nsapy.Nsabase import Analysis

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
