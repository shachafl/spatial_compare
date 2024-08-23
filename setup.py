from pathlib import Path
import setuptools

#!/usr/bin/env python



setuptools.setup(
    name="spatial_compare",
    packages=setuptools.find_packages(),
    install_requires=[
        "anndata",
        "jupyter",
        "seaborn",
        "scanpy",
        "scipy",
        "igraph",
        "leidenalg"
    ],
    include_package_data=True,
    version="0.1",
    python_requires=">=3.9"
)