#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# setup script for the ztf_tiling package.
#
# Author: D. Kaplan (kaplan@uwm.edu)

import os,glob
from setuptools import setup

setup(
    name='ztf_tiling',
    version='0.1',
    description=' ',
    author='David Kaplan',
    author_email='kaplan@uwm.edu',
    packages=['ztf_tiling','ztf_tiling.data'],
    #scripts=glob.glob('scripts/*.py'),
    url = 'https://github.com/dlakaplan/ztf_tiling',
    install_requires=['numpy','astropy'],
    data_files=[('ztf_tiling/data/',glob.glob('data/*.dat'))],
    )
