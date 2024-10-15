#!/usr/bin/env python
from soi.__version__ import version

from setuptools import setup, find_packages
from distutils.extension import Extension
#from Cython.Build import cythonize

with open('README.md') as f:
    long_description = f.read()


setup(
    name='orthoindex',
    version=version,
    description='OrthoIndex: distinguishing synteny from orthology to out-paralogy',
    url='https://github.com/zhangrengang/orthoindex/',
    author='Zhang, Ren-Gang and Wang, Zhao-Xuan',
    license='GPL-3.0',

    python_requires='>=3.7',
    packages=find_packages(),
    include_package_data=True,
    scripts=[],
    entry_points={
        'console_scripts': ['soi = soi.options:main',
							'soi-syn = soi.mcscan:main',
							'soi-orth = soi.OrthoFinder:main',
        ],
    },
)
