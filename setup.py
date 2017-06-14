#!/usr/bin/env python

"""
Author: Luca Masera.
Edited from the PC++ project https://bitbucket.org/francesco-asnicar/gene_network_expansion.
Copyright (C) 2017, all rights reserved
This file (setup.py) is part of the PC++ project.

Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
"""

from distutils.core import setup
from distutils.extension import Extension

setup(
    name="pypcalg",
    version='0.1.1',
    description='Python wrap for the skeleton function of the PC-algorithm implemented in the PC++ project',
    url='https://github.com/lucamasera/pyPCalg',
    author='Luca Masera',
    author_email='luca.masera@unitn.it',
    license='BOOST v1.0',
    ext_modules=[
        Extension(
            "pypcalg", 
            ["pypcalg.cpp"],
            libraries = ["boost_python"]
        )
    ]
)
