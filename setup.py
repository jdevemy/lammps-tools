#!/usr/bin/env python

from distutils.core import setup

setup(name = 'lammps-tools',
      version = '0.1',
      description = 'LAMMPS tools from TIM',
      author = 'Julien Devemy',
      author_email = 'julien.devemy@univ-bpclermont.fr',
      py_modules = ['dlpoly', 'ff', 'lammps', 'zmat', 'molecule', 'utils3d', 'utilsscript', 'xyz', 'element', 'conf'],
      package_dir = {'': 'lib'},
      data_files=[('data', ['data/elements.csv'])],
      scripts=['create_conf', 'zmat2xyz', 'fusion_mols', \
               'fep', 'fdti', 'nti', 'bar']
     )
