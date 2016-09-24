#!/usr/bin/env python
from __future__ import absolute_import
import os
import sys
import re

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

# Handle encoding
major, minor1, minor2, release, serial = sys.version_info
if major >= 3:
    def rd(filename):
        f = open(filename, encoding="utf-8")
        r = f.read()
        f.close()
        return r
else:
    def rd(filename):
        f = open(filename)
        r = f.read()
        f.close()
        return r

setup(    
    name='ExoSOFTmodel', 
    packages =['ExoSOFTmodel'],
    version="1.0.4", 
    author='Kyle Mede',
    author_email = 'kylemede@gmail.com',
    url = 'https://github.com/kylemede/ExoSOFTmodel',
    license = ['GNU GPLv3'],
    description ='An astronomical model for calculating the predicted astrometry and radial velocity due to a companion.',
    long_description = rd("README-pypi.rst"),
    package_data={"": ["LICENSE", "AUTHORS.rst"]},
    include_package_data=True,
    keywords=['model'],
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python'
        ],
    #include_dirs = [np.get_include()],
    install_requires = ['six','KMlogger','numpy','cython','emcee','corner'],
    #ext_modules=[]
)
