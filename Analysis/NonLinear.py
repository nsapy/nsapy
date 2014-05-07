# -*- coding:utf-8 -*-
'''
**Python Module for Nonlinear Structural Analysis
**用于结构非线性分析的Python模块

'''

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
from scipy.optimize import root,newton_krylov

from nsapy.Nsabase import Analysis

def newton(func,U0,jac,tol=1.e-6,maxiter=10,ifprint=False):

    '''
    牛顿法求解非线性方程组
    func ---- 函数,平衡方程,变量为U
    U0   ---- 解的初始估计值
    jac  ---- 切线刚度 \partial{func} / \partial{U}
    tol  ---- 收敛容差(相对)
    maxiter  ---- 最多迭代次数
    '''
    U = U0
    conv = 0
    for niter in xrange(maxiter):
        UF = func(U)
        K = jac(U)
        dU = -spsl.spsolve(K,UF)
        nerr = np.linalg.norm(dU)/np.linalg.norm(U)
        if nerr>tol:
            U = U0+dU
            U0 = U
            print K.todense()[-1,-1],dU[-1],U[-1]
        else:
            break

    if ifprint:
        if niter<maxiter-1:
            print 'Converged after %d iterations, norm error = %.2g %%.'%(niter+1,nerr)
            conv = 1
        else:
            print 'Fail to Converge after %d iterations, norm error = %.2g %%.'%(maxiter,nerr)
    return U,conv

def modified_newton(func,U0,jac,tol=1.e-6,maxiter=10,ifprint=False):

    '''
    修正牛顿法求解非线性方程组(荷载步内不更新刚度矩阵)
    func ---- 函数,平衡方程,变量为U
    U0   ---- 解的初始估计值
    jac  ---- 切线刚度 \partial{func} / \partial{U}
    tol  ---- 收敛容差(相对)
    maxiter  ---- 最多迭代次数
    '''
    U = U0
    conv = 0
    K = jac(U)
    solve = spsl.factorized(K) # Makes LU decomposition.

    for niter in xrange(maxiter):
        UF = func(U)
        dU = -solve(UF) # Uses the LU factors.
        nerr = np.linalg.norm(dU)/np.linalg.norm(U)
        if nerr>tol:
            U = U0+dU
        else:
            break

    if ifprint:
        if niter<maxiter-1:
            print 'Converged after %d iterations, norm error = %.4g %%.'%(niter+1,nerr)
            conv = 1
        else:
            print 'Fail to Converge after %d iterations, norm error = %.4g %%.'%(maxiter,nerr)
    return U,conv

class Analysis_NonLinear(Analysis):
    """Analysis_NonLinear"""
    def __init__(self, domain):
        super(Analysis_NonLinear, self).__init__(domain)
        self.results = []
        self.tol = 1e-6
        self.maxiter = 10

class Analysis_NonLinear_Static(Analysis_NonLinear):
    """非线性静力求解"""

    def __init__(self, domain):
        super(Analysis_NonLinear_Static, self).__init__(domain)
        self.cdof_indx = 0        
        self.cdof_value = 0        

    def balance_func(self,U):

        self.domain.U = U
        # 根据变形计算内力(为得到反力向量Ft)
        self.domain.get_result()
        # 重新组装刚度矩阵并施加荷载、边界条件
        self.domain.reassemble_stiffness_matrix()
        self.domain.apply_load(self.cstep)
        self.domain.apply_cons(self.cstep)
        # 计算残差
        res = self.domain.Ft-self.domain.F

        return res

    def balance_func_no_update(self,U):

        self.domain.U = U
        # 根据变形计算内力(为得到反力向量Ft)
        self.domain.get_result()
        # 计算残差
        res = self.domain.Ft-self.domain.F
        
        return res

    def jac_func(self,U):
        # 以刚度矩阵为雅可比矩阵
        return self.domain.K

    def execute(self,nsteps):
        '''执行多步分析'''
        disnorm = np.linalg.norm(self.domain.load.values()[0].dof_value[0])
        for i in range(nsteps):
            self.cstep = i
            self.ctime = self.domain.load.values()[0].time[i]
            self.domain.apply_load(i)
            self.domain.apply_cons(i)

            # 采用上一步的刚度矩阵试算变形作为迭代初始值
            U0 = self.domain.U+spsl.spsolve(self.domain.K,self.domain.dF)
            # 迭代求解

            U,conv = newton(self.balance_func,U0,self.jac_func,self.tol,self.maxiter,ifprint=True)
            # U,conv = modified_newton(self.balance_func_no_update,U0,self.jac_func,self.tol,self.maxiter,ifprint=True)
            
            self.domain.U = U

            # 更新内力,应力应变状态,刚度矩阵
            self.domain.update()
            # 写入结果文件
            self.write_data_file()

            if len(self.results)>0:
                for ri in self.results:
                    ri.save_result()
            
            print 'Current Time is %s.'%self.ctime

        self.save_data_file()

    def disp_func(self,U):

        self.domain.U = U
        # 根据变形计算内力(为得到反力向量Ft)
        self.domain.get_result()
        # 重新组装刚度矩阵并施加边界条件
        self.domain.reassemble_stiffness_matrix()
        self.domain.apply_cons(self.cstep)
        # 计算残差
        res = self.domain.Ft-self.domain.F
        
        return res

    def disp_control(self,step):

        refU = spsl.spsolve(self.domain.K,self.domain.dF)
        factor = cdof_value/refU[cdof_indx]

        for loadi in self.domain.load.values():
            loadi.dof_value_der*=factor

        self.domain.apply_load(i)
        self.domain.apply_cons(i)

    def execute_disp(self,nsteps):
        '''执行多步分析'''
        for i in range(nsteps):
            self.cstep = i
            self.ctime = self.domain.load.values()[0].time[i]
            self.domain.apply_load(i)
            self.domain.apply_cons(i)
            
            conv = 0
            niter = 0
            while conv == 0 and niter<self.maxiter-1:
                            
                self.disp_control(i)
                # 采用上一步的刚度矩阵试算变形作为迭代初始值
                U0 = self.domain.U+spsl.spsolve(self.domain.K,self.domain.dF)
                # 迭代求解
                fnorm = np.linalg.norm(self.domain.F)
                U,conv = newton(self.balance_func,U0,self.jac_func,fnorm,self.tol,self.maxiter,ifprint=True)
                
                niter += 1
                self.domain.U = U

            # 更新内力,应力应变状态,刚度矩阵
            self.domain.update()
            # 写入结果文件
            self.write_data_file()

            if len(self.results)>0:
                for ri in self.results:
                    ri.save_result()
            
            print 'Current Time is %s.'%self.ctime

        self.save_data_file()