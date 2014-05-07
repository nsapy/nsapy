# distutils: language = c++
# distutils: sources = Concrete01.cpp Concrete02.cpp Steel01.cpp Steel02.cpp
cdef extern from "Concrete01.h":
    cdef cppclass Concrete01:
        Concrete01 (int, double, double, double, double) except +
        int setTrialStrain(double, double)
        int commitState()
        int revertToLastCommit()
        double getStrain()      
        double getStress()
        double getTangent()
        double getInitialTangent()

cdef class PyConcrete01(object):
    cdef Concrete01 *thisptr
    cdef int tag
    cdef double stress,strain,tangent
    property tag:
        def __get__(self): return self.tag
        def __set__(self, tag): self.tag = tag
    property stress:
        def __get__(self): return self.stress
        def __set__(self, stress): self.stress = stress
    property strain:
        def __get__(self): return self.strain
        def __set__(self, strain): self.strain = strain
    property tangent:
        def __get__(self): return self.tangent
        def __set__(self, tangent): self.tangent = tangent
    def __cinit__(self, int mtag, double fpc, double eco, double fpcu, double ecu):
        self.thisptr = new Concrete01(mtag, fpc, eco, fpcu, ecu)
        self.tag = mtag
    def __dealloc__(self):
        del self.thisptr
    def update_history_para(self):
        '''更新历史变量'''
        self.thisptr.commitState()
    def set_strain(self,strain):
        '''获取应变'''
        self.strain = strain
    def get_stress_tangent(self):
        '''根据给定应变值计算应力及切线刚度'''
        self.thisptr.setTrialStrain(self.strain, 0.0)
    def get_key_para(self):
        '''导出用于计算截面,单元状态的关键变量,便于更高层次对象的访问'''
        self.strain = self.thisptr.getStrain()
        self.stress = self.thisptr.getStress()
        self.tangent = self.thisptr.getTangent()

cdef extern from "Concrete02.h":
    cdef cppclass Concrete02:
        Concrete02 (int, double, double, double, double, 
            double, double, double) except +
        int setTrialStrain(double, double)
        int commitState()
        int revertToLastCommit()
        double getStrain()      
        double getStress()
        double getTangent()
        double getInitialTangent()

cdef class PyConcrete02(object): 
    cdef Concrete02 *thisptr
    cdef int tag
    cdef double stress,strain,tangent
    property tag:
        def __get__(self): return self.tag
        def __set__(self, tag): self.tag = tag
    property stress:
        def __get__(self): return self.stress
        def __set__(self, stress): self.stress = stress
    property strain:
        def __get__(self): return self.strain
        def __set__(self, strain): self.strain = strain
    property tangent:
        def __get__(self): return self.tangent
        def __set__(self, tangent): self.tangent = tangent

    def __cinit__(self, int mtag, double _fc, double _epsc0, double _fcu, double _epscu, double _rat, double _ft, double _Ets):
        self.thisptr = new Concrete02(mtag, _fc, _epsc0, _fcu, _epscu, _rat, _ft, _Ets)
        self.tag = mtag
    def __dealloc__(self):
        del self.thisptr
    def update_history_para(self):
        '''更新历史变量'''
        self.thisptr.commitState()
    def set_strain(self,strain):
        '''获取应变'''
        self.strain = strain
    def get_stress_tangent(self):
        '''根据给定应变值计算应力及切线刚度'''
        self.thisptr.setTrialStrain(self.strain, 0.0)
    def get_key_para(self):
        '''导出用于计算截面,单元状态的关键变量,便于更高层次对象的访问'''
        self.strain = self.thisptr.getStrain()
        self.stress = self.thisptr.getStress()
        self.tangent = self.thisptr.getTangent()

cdef extern from "Steel01.h":
    cdef cppclass Steel01:
        Steel01 (int, double, double, double, 
            double, double, double, double) except +
        int setTrialStrain(double, double)
        int commitState()
        int revertToLastCommit()
        double getStrain()      
        double getStress()
        double getTangent()
        double getInitialTangent()

cdef class PySteel01(object): 
    cdef Steel01 *thisptr
    cdef int tag
    cdef double stress,strain,tangent
    property tag:
        def __get__(self): return self.tag
        def __set__(self, tag): self.tag = tag
    property stress:
        def __get__(self): return self.stress
        def __set__(self, stress): self.stress = stress
    property strain:
        def __get__(self): return self.strain
        def __set__(self, strain): self.strain = strain
    property tangent:
        def __get__(self): return self.tangent
        def __set__(self, tangent): self.tangent = tangent

    def __cinit__(self, int mtag, double FY, double E, double B,
                    double A1, double A2, double A3, double A4):
        self.thisptr = new Steel01(mtag, FY, E, B, A1, A2, A3, A4)
        self.tag = mtag
    def __dealloc__(self):
        del self.thisptr
    def update_history_para(self):
        '''更新历史变量'''
        self.thisptr.commitState()
    def set_strain(self,strain):
        '''获取应变'''
        self.strain = strain
    def get_stress_tangent(self):
        '''根据给定应变值计算应力及切线刚度'''
        self.thisptr.setTrialStrain(self.strain, 0.0)
    def get_key_para(self):
        '''导出用于计算截面,单元状态的关键变量,便于更高层次对象的访问'''
        self.strain = self.thisptr.getStrain()
        self.stress = self.thisptr.getStress()
        self.tangent = self.thisptr.getTangent()

cdef extern from "Steel02.h":
    cdef cppclass Steel02:
        Steel02 (int,
        double, double, double,
        double, double, double,
        double, double, double, double, double) except +
        int setTrialStrain(double, double)
        int commitState()
        int revertToLastCommit()
        double getStrain()      
        double getStress()
        double getTangent()
        double getInitialTangent()

cdef class PySteel02(object): 
    cdef Steel02 *thisptr
    cdef int tag
    cdef double stress,strain,tangent
    property tag:
        def __get__(self): return self.tag
        def __set__(self, tag): self.tag = tag
    property stress:
        def __get__(self): return self.stress
        def __set__(self, stress): self.stress = stress
    property strain:
        def __get__(self): return self.strain
        def __set__(self, strain): self.strain = strain
    property tangent:
        def __get__(self): return self.tangent
        def __set__(self, tangent): self.tangent = tangent

    def __cinit__(self, int mtag, double fy, double E0, double b,
        double R0, double cR1, double cR2,
        double a1, double a2, double a3, double a4, double sigInit):
        self.thisptr = new Steel02(mtag,fy,E0,b,R0,cR1,cR2,a1,a2,a3,a4,sigInit)
        self.tag = mtag
    def __dealloc__(self):
        del self.thisptr
    def update_history_para(self):
        '''更新历史变量'''
        self.thisptr.commitState()
    def set_strain(self,strain):
        '''获取应变'''
        self.strain = strain
    def get_stress_tangent(self):
        '''根据给定应变值计算应力及切线刚度'''
        self.thisptr.setTrialStrain(self.strain, 0.0)
    def get_key_para(self):
        '''导出用于计算截面,单元状态的关键变量,便于更高层次对象的访问'''
        self.strain = self.thisptr.getStrain()
        self.stress = self.thisptr.getStress()
        self.tangent = self.thisptr.getTangent()