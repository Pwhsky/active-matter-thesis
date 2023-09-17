from setuptools import setup,Extension
from Cython.Build import cythonize
import numpy as np
extensions = [Extension("cython_functions", sources=["cython_functions.pyx"], extra_compile_args=["-Ofast", "-march=native"])]


setup(
     ext_modules=cythonize(extensions),
    include_dirs=[np.get_include()],
)
