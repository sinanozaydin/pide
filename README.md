<img src="./docs/figures/pide_logo.png">

[![PyPI version](https://img.shields.io/pypi/v/pide.svg)](https://pypi.org/project/pide/) [![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
 

# About

`pide` is a Python3 library for calculating geophysical parameters (e.g., electrical conductivity, seismic velocity), employing the results from experimental petrology, mineral/rock physics, and thermomechanical modelling studies. `pide` can calculate the theoretical electrical conductivity of any earth material that exists in the literature. `pide` can also calculate seismic velocity utilising the external 'sister' library `santex`. Using these theoretical calculations, users can utilise inversion modules to decode geophysical anomalies compositionally or convert thermomechanical models into geophysical observables. With a given spatial mapping of earth materials, which can preferentially be loaded from a thermomechanical model, `pide`  can be used to build synthetic electrical conductivity and seismic velocity models and generate gravity and magnetic anomalies. Moreover, it is built as a modular tool, so users can easily build their functions.

# Installation

## Linux and MacOS

To install `pide`, the user can simply go their terminal and type the following command:

```bash
pip install pide
```
or alternatively, they can clone the repository, then go to the directory of the source with `cd` and perform:

```bash
pip install .
```

If you want to help develop and change code as it is being used:

```bash
pip install -e .
```
Building from the source is encouraged at this point, since `pide` will be in development stage in the following years (2024).

## Windows

On Windows computers, installation of `pide` may require a few more steps. If the user fails to install with `pip install pide` command, it will likely be caused by the `pycifrw` library, which is required by the dependency `santex` through `orix`. Then, we encourage the users to install the `pycifrw` first if `pip install pide` fails. This, for instance, can be achieved through utilising conda package manager with:

```bash
conda install pycifrw -c conda-forge
```

Alternatively, the user can install `pide` by manually installing the dependencies. This will involve the execution of the following commands:

```bash
pip install numpy matplotlib scipy h5py harmonica pyproj netCDF4 psutil
pip install --no-deps santex
```

then:

```bash
pip install --no-deps pide
```

or alternatively, they can clone the repository, then go to the directory of the source with `cd` and perform:

```bash
pip install --no-deps .
```

# Workflow and how to use

How to use `pide` can be learned through Jupyter notebooks provided in examples/notebooks directory. The general workflow can be tracked through the chart below:

<img src="./docs/paper/figures/pide_workflow.png">

Information on all methods (input/output, examples) can be accesed through the pide object method as follows:
<pre>

pide_object = pide.pide()

pide_object.list_methods()
pide_object.get_method_manual(method_name)

</pre>

# Getting Involved

Since `pide` is an open-source library, users are encouraged to be contribute and become developers of the project. For further information about how to contribute, please refer to the [Contributing Guide](https://github.com/sinanozaydin/pide/blob/post_joss/CONTRIBUTING.md).

# Running Tests

To run the tests for the package, simply go to the test files and run them with python command:

```bash
python3 test_electrical_cond.py
```

# Contacts

| **Sinan Özaydın** | sinan.ozaydin@protonmail.com | sinan.ozaydin@sydney.edu.au
