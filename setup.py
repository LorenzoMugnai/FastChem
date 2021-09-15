
import setuptools
from setuptools import find_packages
from setuptools import setup, Extension
import re, os

packages = find_packages(exclude=('tests', 'doc'))
provides = ['taurex_fastchem', ]

requires = []

install_requires = ['taurex', 'cython',]

entry_points = {'taurex.plugins': 'fastchem = taurex_fastchem'}

from Cython.Build import cythonize
from Cython.Distutils import build_ext
import glob
import numpy as np
import shutil
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
                      extra_compile_args=["-O3", "--std=c++11",],
                      
                      language="c++"
                      )

files = glob.glob(os.path.join('taurex_fastchem','data','input','*.dat'))

if os.path.exists('input'):
     try:
          if os.path.exists(os.path.join('taurex_fastchem','data','input')):
               os.unlink(os.path.join('taurex_fastchem','data','input'))
          shutil.copytree('input', os.path.join('taurex_fastchem','data','input') )
          shutil.copy(os.path.join('fastchem_src','chem_input','chemical_elements.dat'), os.path.join('taurex_fastchem','data','input') )
     except Exception:
          pass


files =  ([f[len('taurex_fastchem')+1:] for f in files])

data_files = [*files, os.path.join('fastchem_src','chem_input','chemical_elements.dat') ]

version = "1.0.0-dev0"

with open("README.md", "r") as fh:
    long_description = fh.read()

pos = long_description.find('# TauREx')

long_description = long_description[pos:]


ext = cythonize([clib],language_level="3")
setup(name='taurex_fastchem',
      author="Ahmed Faris Al-Refaie",
      author_email="ahmed.al-refaie.12@ucl.ac.uk",
      license="BSD",
      description='Python Wrapper for FastChem chemical scheme',
      packages=packages,
      ext_modules=ext,
      entry_points=entry_points,
      provides=provides,
      requires=requires,
      version=version,
      keywords=['exoplanet',
                'chemistry'
                'taurex',
                'plugin',
                'taurex3',
                'atmosphere',
                'atmospheric'],
      long_description=long_description,
      long_description_content_type='text/markdown',
      package_data={'taurex_fastchem': data_files},
 #     cmdclass = {"build_ext": build_ext},
      install_requires=install_requires)