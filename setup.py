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
    #url = 'https://github.com/MatteoGiomi/dataslicer',
    install_requires=['numpy','healpy','astropy'],
    data_files=[('ztf_tiling/data/',glob.glob('data/*.dat'))],
    )
