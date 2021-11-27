# -*- coding: utf-8 -*-

# stdlib modules
import setuptools
# from sys import path

# # internal modules (add to path, so that __init__.py is not executed)
# path.append("pyasim")
# from version import version

# with open("README.md", "r") as fh:
#     long_description = fh.read()

version = "0.1a1"
long_description = "TODO: write this"

packages = setuptools.find_packages(include=["mushell.*", "mushell"])

setuptools.setup(
    name="mushell",
    version=version,
    author="Volker Siepmann",
    author_email="volker.siepmann@yara.com",
    description="Equation oriented process modelling library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://yara-alm.visualstudio.com/Production.YTC/_git/Pyasim",  # TODO: github repo
    packages=packages,
    # package_data={'mushell': ['units.xml']},
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.9',
    install_requires=[
        "casadi",
        "pyyaml"
        ],
    extras_require={"doc": ["Sphinx>=2.2",
                            "sphinxcontrib-bibtex<2.0"],
                    # above: bug in V2, causing error (TODO: test now and then to relax this again)
                    "test": ["pytest>=5.3"]}
)
