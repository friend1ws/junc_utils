#!/usr/bin/env python

from distutils.core import setup

setup(name='junc_utils',
      version='0.2.0',
      description='Python programs for analyzing splicing junctions',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/sv_utils',
      package_dir = {'': 'lib'},
      packages=['junc_utils'],
      scripts=['junc_utils'],
      license='GPL-3'
     )

