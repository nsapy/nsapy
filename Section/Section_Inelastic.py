# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis
**用于结构非线性分析的Python模块

'''

import numpy as np
import scipy.sparse as sps

from nsapy.Nsabase import Section

class Section_Inelastic(Section):
    """Section_Inelastic"""
    def __init__(self, *arg):
        super(Section_Inelastic, self).__init__(*arg)
        self.E0 = self.mat.E0
        self.A = arg[2] if len(arg)>=3 else 0.0
        self.Da = self.E0*self.A

        self.inner_force = np.zeros(5)

    def update(self):
        self.mat.update_history_para()

class Truss(Section_Inelastic):
    """Truss"""
    def __init__(self, *arg):
        super(Truss, self).__init__(*arg)

    def get_force_stiffness(self,strain):

        # 将截面应变传递给材料
        self.mat.set_strain(strain[0])
        # 计算材料应力及切线刚度
        self.mat.get_stress_tangent()
        # 由材料状态计算截面内力及切线刚度
        self.inner_force = np.array([[self.mat.stress*self.A]])
        self.Da = self.mat.tangent*self.A

class Fiber(Truss):
    """Fiber"""
    def __init__(self, *arg):
        super(Fiber, self).__init__(*arg)
        self.coord = arg[3] if len(arg)>=4 else np.zeros(2)

class FiberSection(Section):
    """docstring for FiberSection"""
    def __init__(self, *arg):
        self.tag = arg[0]
        self.fiber = {}

    def add_fiber(self,fiber):
        self.fiber[fiber.tag] = fiber

    def get_force_stiffness(self,strain):

        # 将截面应变传递给材料
        self.mat.set_strain(strain[0])
        # 计算材料应力及切线刚度
        self.mat.get_stress_tangent()
        # 由材料状态计算截面内力及切线刚度
        self.inner_force = np.array([[self.mat.stress*self.A]])
        self.Da = self.mat.tangent*self.A
            