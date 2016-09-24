ExoSOFTmodel
============

<!--[![Build Status](https://travis-ci.org/kylemede/KMlogger.svg?branch=master)](https://travis-ci.org/kylemede/KMlogger)-->
<!--[![PyPI version](https://badge.fury.io/py/KMlogger.svg)](https://badge.fury.io/py/KMlogger)-->
[![License](https://img.shields.io/badge/license-GPL-blue.svg)](https://github.com/kylemede/ExoSOFTmodel/blob/master/LICENSE)
<!--[![Coverage Status](https://coveralls.io/repos/github/kylemede/KMlogger/badge.svg?branch=master)](https://coveralls.io/github/kylemede/KMlogger?branch=master)-->

An astronomical model for calculating the predicted astrometry and radial velocity due to a companion.

MORE INFORMATION TO COME.

Examples made both showing how to use ExoSOFTmodel with the popular emcee, and 
ExoSOFT.


Installation Notes
==================
Only uses standardized Python libraries and KMlogger, which is also hosted on PyPi.  

NOT ON PIP YET!!
----------------

Thus, the easiest way to install it is:
 
 $pip install ExoSOFTmodel
 
Re-compiling Cython tools
-------------------------

Should you need to modify any of the tools within [cytools.pyx](https://github.com/kylemede/ExoSOFTmodel/blob/master/ExoSOFTmodel/cytools.pyx)
 you need to have cython installed.  Then compile the code within its directory using:
 
 $python setup.py build_ext --inplace

License
-------

Copyright 2016 Kyle Mede and contributors.

ExoSOFTmodel is free software made available under the GNU GPLv3 license. 
For details see the license.txt file.