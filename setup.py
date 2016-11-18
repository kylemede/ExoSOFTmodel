#!/usr/bin/env python
from __future__ import absolute_import
import os
import sys
#import re


try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup
    
from distutils.extension import Extension

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

setup(    
    name='ExoSOFTmodel', 
    packages =['ExoSOFTmodel'],
    version="1.1.2", 
    author='Kyle Mede',
    author_email = 'kylemede@gmail.com',
    url = 'https://github.com/kylemede/ExoSOFTmodel',
    license = ['GNU GPLv3'],
    description ='An astronomical model for calculating the predicted astrometry and radial velocity due to a companion.',
    long_description = "For examples of how to use ExoSOFTmodel or"+
    " read its documentation, please visit:\nhttps://github.com/kylemede/ExoSOFTmodel",
    package_data={"": ["LICENSE", "AUTHORS.rst","ExoSOFTmodel/cytools.pyx"]},
    include_package_data=True,
    keywords=['model'],
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python'
        ],
    #include_dirs=['.','./ExoSOFTmodel'],
    install_requires = ['six','KMlogger','numpy','emcee','corner'],
    ext_modules=[ Extension("cytools",["ExoSOFTmodel/cytools.c"]) ]
)
