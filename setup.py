#!/usr/bin/env python

"""
Author: Luca Masera.
Edited from the PC++ project https://bitbucket.org/francesco-asnicar/gene_network_expansion.
Copyright (C) 2017, all rights reserved
This file (setup.py) is part of the PC++ project.
PyPCalg is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""

from distutils.core import setup
from distutils.extension import Extension

setup(
	name="pypcalg",
    ext_modules=[
        Extension(
        	"pypcalg", 
        	["pypcalg.cpp"],
        	libraries = ["boost_python"]
        )
    ]
)
