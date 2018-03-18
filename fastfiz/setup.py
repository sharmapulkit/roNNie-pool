#!/usr/bin/env python

"""
setup.py file for fastfiz
"""

from distutils.core import setup, Extension
import glob

fastfiz_module = Extension('_fastfiz',
                           sources=['src/fastfiz.i',
                                    'src/FastFiz.cpp',
                                    'src/Noise.cpp',
                                    'src/Rules.cpp',
                                    'src/Stopwatch.cpp'],
                           include_dirs=['include'],
                           swig_opts=['-c++', '-Wall', '-Iinclude'],
                           extra_compile_args=['-std=c++11'],
                           extra_link_args=['-lgsl', '-lgslcblas'],
                           )

setup (name = 'fastfiz',
       description = """fastfiz Billiards engine wrapper""",
       ext_modules = [fastfiz_module],
       py_modules = ["fastfiz"],
       )
