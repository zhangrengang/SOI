#!/usr/bin/env python
from soi.__version__ import version

from setuptools import setup, find_packages
from distutils.extension import Extension

with open('README.md') as f:
    long_description = f.read()


setup(
    name='soi',
    version=version,
    description='SOI: identifying orthologous synteny',
    url='https://github.com/zhangrengang/soi/',
    author='Zhang, Ren-Gang',
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
