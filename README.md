ligandsplitter
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/ligandsplitter/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/ligandsplitter/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ligandsplitter/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ligandsplitter/branch/main)


A Python package for splitting, creating, and validating ligand files

### Usage

This repository is currently not meant to be used as a standalone program and was made to enhance a series of notebooks that explore molecular docking using Jupyter notebooks. This repository, called basil_dock, is currently being developed and will be released in the near future.

### Installation

When cloning the basil_dock, make sure to initialize the submodule. This can be done in one of two ways:
1. Clone basil_dock before initializing submodule
```
git clone <url to basil_dock>
```
```
cd ligandsplitter
```
```
git submodule update --init --recursive
```
2. Clone basil_dock and initialize submodule at the same time
```
git clone --recurse-submodules <url to basil_dock>
```
### Copyright

Copyright (c) 2024, Lee Schoneman


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.10.
