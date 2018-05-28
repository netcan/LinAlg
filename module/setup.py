#!/usr/bin/env python
# encoding: utf-8

from distutils.core import setup, Extension

linalgmodule = Extension('linalg',
        extra_compile_args = ['-std=c++14'],
        sources = ['linalgmodule.cpp'],
        libraries = ['linalg'],
        )

setup (name = 'linalg',
        version = '0.1',
        author = 'netcan',
        author_email = 'netcan1996@gmail.com',
        url = 'http://www.netcan666.com/2018/05/27/%E4%BB%8E%E5%A4%B4%E5%BC%80%E5%A7%8B%E5%AE%9E%E7%8E%B0%E4%B8%80%E4%B8%AA%E7%BA%BF%E6%80%A7%E4%BB%A3%E6%95%B0%E5%BA%93/',
        description = 'A Linear Algebra library for studying by netcan',
        ext_modules = [linalgmodule]
        )
