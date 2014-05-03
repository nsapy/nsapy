# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis
**用于结构非线性分析的Python模块

'''

import numpy as np
import scipy.sparse as sps
from copy import deepcopy

from nsapy.Nsabase import Material1D

class Material_Inelastic(Material1D):
    """非线性单轴材料"""
    def __init__(self, *arg):
        super(Material_Inelastic,self).__init__(*arg)
        self.E0 = float(arg[1]) # 弹性模量
        
        # 初始化状态变量
        self.init_para()
        
    def init_para(self):

        self.para_keys = ['strain','stress','tangent','dstrain']
        self.para = {para:0.0 for para in self.para_keys} # 状态变量记录在此字典内
        self.para['tangent'] = self.E0

        self.his_para = deepcopy(self.para) # 历史状态变量记录在此字典内

    def update_history_para(self):
        '''更新历史变量'''
        self.his_para = deepcopy(self.para)

    def set_strain(self,strain):
        '''获取应变'''
        self.para['strain'] = strain
        self.para['dstrain'] = self.para['strain']-self.his_para['strain']

    def get_stress_tangent(self):
        '''根据给定应变值计算应力及切线刚度'''
        pass

    def get_key_para(self):
        '''导出用于计算截面,单元状态的关键变量,便于更高层次对象的访问'''
        self.strain = self.para['strain']
        self.stress = self.para['stress']
        self.tangent = self.para['tangent']

class Bilinear(Material_Inelastic):
    """双线性弹塑性材料"""
    def __init__(self, *arg):
        super(Bilinear, self).__init__(*arg)
        self.fy = float(arg[2]) # 屈服应力
        self.alpha = float(arg[3]) # 屈服后刚度比
        self.Ep = self.E0*self.alpha # 屈服后刚度
        self.sy = self.fy/self.E0 # 屈服应变

    def get_stress_tangent(self):
        '''根据给定应变值计算应力及切线刚度'''

        # 确定屈服面(线)
        bound_pos = self.fy+self.Ep*(self.para['strain']-self.sy)
        bound_neg = -self.fy+self.Ep*(self.para['strain']+self.sy)
        # 根据前一步的刚度估算应力值
        try_stress = self.his_para['stress']+self.para['dstrain']*self.his_para['tangent']

        # 加载
        if self.para['dstrain'] > 0.0:
            if try_stress < bound_pos: # 弹性状态
                self.para['tangent'] = self.E0
                # 根据更新的刚度值重新计算应力值
                if self.para['tangent'] != self.his_para['tangent']:
                    self.para['stress'] = min(self.his_para['stress']+self.para['dstrain']*self.para['tangent'],bound_pos)
                else:
                    self.para['stress'] = try_stress
            else:
                self.para['stress'] = bound_pos
                self.para['tangent'] = self.Ep
        # 卸载
        elif self.para['dstrain'] < 0.0:
            if try_stress > bound_neg:
                self.para['tangent'] = self.E0
                if self.para['tangent'] != self.his_para['tangent']:
                    self.para['stress'] = max(self.his_para['stress']+self.para['dstrain']*self.para['tangent'],bound_neg)
                else:
                    self.para['stress'] = try_stress
            else:
                self.para['stress'] = bound_neg
                self.para['tangent'] = self.Ep

        else:
            pass

        self.get_key_para()