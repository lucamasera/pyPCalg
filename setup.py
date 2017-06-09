#!/usr/bin/env python

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
