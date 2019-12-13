from __future__ import print_function
import argparse
import os
import sys
import numpy
from setuptools import find_packages, setup, Extension
from setuptools.command.test import test as test_command


setup(
    name="mtbench",
    license='BSD2',
    author="Ryan Modrak",
    url="https://github.com/uafgeotools/mtbench",
    packages=find_packages(),
    zip_safe=False,
    keywords=[
        "seismology"
    ],
    python_requires='~=3.5',
    install_requires=[
        "numpy", "scipy", "obspy", 
        "h5py", "retry", "flake8>=3.0", "pytest", "nose",
    ],
)

