
import setuptools
from setuptools import find_packages
from distutils.core import setup
from distutils.core import Extension
from distutils import log
import re, os

packages = find_packages(exclude=('tests', 'doc'))
provides = ['taurex_fastchem', ]

requires = []

install_requires = ['taurex', 'cython',]

entry_points = {'taurex.plugins': 'fastchem = taurex_fastchem'}

from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import glob
import numpy as np
sources = glob.glob('fastchem_src/*.cpp')

clib = Extension("taurex_fastchem.external.fastchem",  # indicate where it should be available !
                      sources=["taurex_fastchem/external/fastchem.pyx",
                               *sources,
                               ],
                      include_dirs=['fastchem_src/'],
                      extra_compile_args=["-O3", "--std=c++11", "-ffast-math", "-Wall"],
                      language="c++"
                      )
ext = cythonize([clib],language_level="3")
setup(name='taurex_fastchem',
      author="Ahmed Faris Al-Refaie",
      author_email="ahmed.al-refaie.12@ucl.ac.uk",
      license="BSD",
      description='Plugin to use FastChem equilibrium chemistry ',
      packages=packages,
      ext_modules=ext,
      entry_points=entry_points,
      provides=provides,
      requires=requires,
 #     cmdclass = {"build_ext": build_ext},
      install_requires=install_requires)