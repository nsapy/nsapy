from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "osmat",
    ext_modules = cythonize('osmat.pyx'),
)
