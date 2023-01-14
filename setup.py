# -*- coding: utf-8 -*-

# stdlib modules
import setuptools
# from sys import path

# # internal modules (add to path, so that __init__.py is not executed)
# path.append("pyasim")
# from version import version

# with open("README.md", "r") as fh:
#     long_description = fh.read()

from src.simu._version import VERSION

long_description = "TODO: write this"

setuptools.setup(
    name="simu",
    version=VERSION,
    author="Volker Siepmann",
    author_email="volker.siepmann@gmail.com",
    description="Equation oriented process modelling library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/VolkerSiep/SiMu",
    # package_data={'simu': ['units.xml']},
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.9',
    install_requires=["numpy", "scipy", "casadi", "pyyaml"],
    extras_require={
        "dev": [
            "Sphinx>=2.2", "sphinxcontrib-bibtex>=2.4", "pytest>=5.3",
            "matplotlib"
        ]
    })
