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

class Analysis_Eigen(Analysis):
    """Analysis_Linear"""
    def __init__(self, domain):
        # super(Analysis_Eigen, self).__init__(domain)
        self.domain = domain

    def execute(self, k):
        Omega_2,self.domain.mode = spsl.eigsh(self.domain.K.tocsc(),k,self.domain.M.tocsc())
        self.domain.Omega = np.sqrt(Omega_2)
        self.domain.T = 2.*np.pi/self.domain.Omega

class Analysis_Linear_Static(Analysis_Linear):
    """Analysis_Linear_Static"""
    def __init__(self, domain):
        super(Analysis_Linear_Static, self).__init__(domain)

    def execute(self,nsteps):
        for i in range(nsteps):
            self.cstep = i
            self.ctime = self.domain.load.values()[0].time[i]
            self.domain.apply_load(i)
            self.domain.apply_cons(i)
            self.domain.U = spsl.spsolve(self.domain.K,self.domain.F)
            self.domain.get_result()

            print 'Current Time is %s.'%self.ctime
            self.write_data_file()

        self.save_data_file()
