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

    def set_strain(self,strain):
        self.strain = strain

    def update_history_para(self):
        self.mat.update_history_para()

class Truss(Section_Inelastic):
    """Truss"""
    def __init__(self, *arg):
        super(Truss, self).__init__(*arg)
    
    def get_force_stiffness(self):

        # 将截面应变传递给材料
        self.mat.set_strain(self.strain[0])
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

    def get_force_stiffness(self):

        # 将截面应变传递给材料
        self.mat.set_strain(self.strain)
        # 计算材料应力及切线刚度
        self.mat.get_stress_tangent()
        # 由材料状态计算截面内力及切线刚度
        self.inner_force = self.mat.stress*self.A
        self.Da = self.mat.tangent*self.A

class FiberTruss(Section):
    """FiberTruss"""
    def __init__(self, *arg):
        self.tag = arg[0]
        self.fiber = {}
        self.Da = 0.0

    def add_fiber(self,*arg):
        self.fiber[arg[0]] = Fiber(*arg)

    def get_section_geometry(self):
        '''计算截面几何参数：面积及形心位置'''
        self.A = 0.0
        for fi in self.fiber.values():
            self.A += fi.A

    def set_strain(self,strain):
        self.strain = strain

    def get_force_stiffness(self):

        '''计算截面刚度及内力'''
        # 初始置零
        self.Da = 0.0
        self.Db = 0.0
        self.inner_force = np.zeros(2)
        # 遍历所有纤维
        for fi in self.fiber.values():
            # 计算纤维中心点的相对坐标
            y = fi.coord[0]-self.y0
            # 计算纤维应变
            fstrain = self.strain[0]-y*self.strain[1]
            # 将应变传递给纤维
            fi.set_strain(fstrain)
            # 计算纤维内力及切线刚度
            fi.get_force_stiffness()
            # 累加计算截面内力(轴力、弯矩)及切线刚度
            self.inner_force[0] += fi.inner_force
            self.Da += fi.Da
            self.inner_force[1] += fi.inner_force*(-y)
            self.Db += fi.Da**y**2

    def update_history_para(self):
        for fi in self.fiber.values():
            fi.mat.update_history_para()

class FiberSection2D(Section):
    """FiberSection2D"""
    def __init__(self, *arg):
        self.tag = arg[0]
        # 截面纤维
        self.fiber = {}
        # 截面刚度矩阵
        self.D = np.zeros((2,2))
        # 截面应变向量
        self.strain = np.zeros(2)
        
    def add_fiber(self,*arg):
        self.fiber[arg[0]] = Fiber(*arg)
        self.get_section_geometry()
        self.get_initial_stiffness()

    def add_rect(self,b,h,nb,nh,mat):
        '''添加矩形截面纤维'''
        db = b/nb
        dh = h/nh
        a = db*dh
        for i in xrange(nb):
            for j in xrange(nh):
                tag = 1000*i+j
                coord = [dh/2.+dh*j,db/2.+db*i]
                self.fiber[tag] = Fiber(tag,mat,a,coord)
        self.get_section_geometry()
        self.get_initial_stiffness()

    def get_section_geometry(self):
        '''计算截面几何参数：面积及形心位置'''
        self.A = 0.0
        yA = 0.0
        for fi in self.fiber.values():
            self.A += fi.A
            yA += fi.A*fi.coord[0]
        self.y0 = yA/self.A

    def set_strain(self,strain):
        '''设置截面应变'''
        self.strain = strain

    def get_initial_stiffness(self):
        '''计算截面初始刚度'''
        # 初始置零
        self.D = np.zeros((2,2))
        # 遍历所有纤维
        for fi in self.fiber.values():
            # 计算纤维中心点的相对坐标
            y = fi.coord[0]-self.y0
            # 计算纤维应变
            fstrain = self.strain[0]-y*self.strain[1]
            # 将应变传递给纤维
            fi.set_strain(fstrain)
            # 计算纤维内力及切线刚度
            fi.get_force_stiffness()
            # 累加纤维刚度得到截面刚度矩阵
            self.D[0,0] += fi.Da
            self.D[0,1] += -fi.Da*y
            self.D[1,1] += fi.Da*y**2
        self.D[1,0] = self.D[0,1]

    def get_force_stiffness(self):
        '''计算截面刚度及内力'''
        # 初始置零
        self.D = np.zeros((2,2))
        self.inner_force = np.zeros(2)
        # 遍历所有纤维
        for fi in self.fiber.values():
            # 计算纤维中心点的相对坐标
            y = fi.coord[0]-self.y0
            # 计算纤维应变
            fstrain = self.strain[0]-y*self.strain[1]
            # 将应变传递给纤维
            fi.set_strain(fstrain)
            # 计算纤维内力及切线刚度
            fi.get_force_stiffness()
            # 累加计算截面内力(轴力、弯矩)及切线刚度
            self.inner_force[0] += fi.inner_force
            self.inner_force[1] += fi.inner_force*(-y)
            self.D[0,0] += fi.Da
            self.D[0,1] += -fi.Da*y
            self.D[1,1] += fi.Da*y**2
        self.D[1,0] = self.D[0,1]

    def update_history_para(self):
        for fi in self.fiber.values():
            fi.update_history_para()