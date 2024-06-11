#!/usr/bin/env python

from pathlib import Path

import setuptools

setuptools.setup(name ="spatial_compare",
    packages=setuptools.find_packages(),
    install_requires=[],
    include_package_data=True,
    version=0.1,
    python_requires='>=3.9'
)