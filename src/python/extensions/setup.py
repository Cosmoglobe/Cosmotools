#/usr/bin/env python2.6
import os
import sys

if not os.path.exists('alpha.o') and "clean" not in sys.argv:
    print "Hand-compiling alpha.o; I don't know how to automate this"
    os.system("gfortran -fPIC -O3 -o quick_qpoint.o -c quick_qpoint.f90")
    os.system("gfortran -fPIC -O3 -o alpha.o -c alpha.f90")

if os.path.exists('alpha.o') and "clean" in sys.argv:
	os.remove("alpha.o")


try:
    import setuptools
except ImportError:
    pass
from numpy.distutils.core import setup, Extension

include = ["-I/usit/titan/u1/joez/usr/lib64/python2.4/site-packages/numpy/core/include/numpy"]

clib = Extension(	name='_quietPyExt',
					sources=['_quietPyExt.c'],
					extra_objects=['alpha.o','quick_qpoint.o'],
					libraries=['gfortran'],
					extra_compile_args=include)


version = "0.1"

distrib = setup(
        version=version,
        author="Joe Zuntz",
        author_email="jaz@astro.ox.ac.uk",
        description="Helper for slice updating",
        license="Academic Free License",
        name="quietPyExt",
        packages=["quietPyExt"],
        ext_modules = [clib],
)
