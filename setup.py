from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy

Ext_Modules=[
    Extension("polyArea", ["polyarea.pyx"],
    language='c++')
]

setup(
    name = 'polyarea app',
    cmdclass = {'build_ext': build_ext},
    ext_modules=Ext_Modules,
    include_dirs=[numpy.get_include()]
)    
