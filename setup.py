
import setuptools
from setuptools import find_packages
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils import log
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
sources = ['fastchem_src/solver_linsol_quadsol.cpp',
 'fastchem_src/solver_nelder_mead_single.cpp',
 'fastchem_src/init_solver.cpp',
 'fastchem_src/solve_fastchem.cpp',
 'fastchem_src/init_read_files.cpp',
 'fastchem_src/mass_action_constant.cpp',
 'fastchem_src/initial_h_density.cpp',
 'fastchem_src/calc_total_h_density.cpp',
 'fastchem_src/calc_densities_private.cpp',
 'fastchem_src/calc_densities.cpp',
 'fastchem_src/fastchem.cpp',
 'fastchem_src/calc_mean_mol_weight.cpp',
 'fastchem_src/solver_newtsol.cpp',
 'taurex_fastchem/glue/init_add_species.cpp',
 'fastchem_src/calc_densities_ph.cpp',
 'fastchem_src/solver_scaling_factor.cpp',
 'fastchem_src/init.cpp',
 'fastchem_src/check.cpp',
 'fastchem_src/bisection.cpp',
 'fastchem_src/get.cpp',
 'fastchem_src/solver.cpp']

clib = Extension("taurex_fastchem.external.fastchem",  # indicate where it should be available !
                      sources=["taurex_fastchem/external/fastchem.pyx",
                               *sources,
                               ],
                      include_dirs=['fastchem_src/'],
                      extra_compile_args=["-O3", "--std=c++11", "-ffast-math", "-Wall","-pedantic","-ffast-math" ,"-march=native"],
                      language="c++"
                      )

input_data = glob.glob('input/*.dat')

data_files = ( 'taurex_fastchem/external/data' , [*input_data, os.path.join('fastchem_src','chem_input','chemical_elements.dat') ] )

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
      data_files=[data_files],
 #     cmdclass = {"build_ext": build_ext},
      install_requires=install_requires)