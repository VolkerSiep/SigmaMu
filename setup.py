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

packages = setuptools.find_packages(include=["simu.*", "simu"])

with open("requirements.txt", encoding="utf-8") as file:
    requirements = file.readlines()

setuptools.setup(
    name="simu",
    version=version,
    author="Volker Siepmann",
    author_email="volker.siepmann@gmail.com",
    description="Equation oriented process modelling library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/VolkerSiep/SiMu",
    packages=packages,
    # package_data={'simu': ['units.xml']},
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.9',
    install_requires=requirements,
    extras_require={
        "dev": [
            "Sphinx>=2.2",
            "sphinxcontrib-bibtex>=2.4",
            "pytest>=5.3"]
    })
