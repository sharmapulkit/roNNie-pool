#!/usr/bin/env python

"""
setup.py file for fastfiz
"""

from distutils.core import setup, Extension
import glob

fastfiz_module = Extension('fastfiz._fastfiz',
                           sources=['swig/fastfiz.i',
                                    'src/FastFiz.cpp',
                                    'src/Noise.cpp',
                                    'src/Stopwatch.cpp',
                                    'src/Rules.cpp',
                                    ],
                           include_dirs=['include'],
                           swig_opts=['-c++', '-Wall', '-I./include', '-outdir', '.'],
                           extra_compile_args=['-std=c++11'],
                           extra_link_args=['-lgsl', '-lgslcblas'],
                           )
# rules_module = Extension('fz._rules',
#                           sources=['rules.i',
#                                   'src/Rules.cpp',
#                                   'src/FastFiz.cpp'],
#                           include_dirs=['.', 'include'],
#                           swig_opts=['-c++', '-Wall', '-I./include', '-outdir', '.'],
#                           extra_compile_args=['-std=c++11'],
#                           extra_link_args=['-lgsl', '-lgslcblas'],
#                           )

setup (name = 'fastfiz',
       description = """fastfiz Billiards engine wrapper""",
       ext_modules = [fastfiz_module],
#       py_modules = ["fz"],
       )
