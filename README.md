# RGE++

A C++ library for the numerical integration of renormalisation group equations (RGEs).

## Features

**RGE++** allows for a fast numerical evaluation of RGEs for various models of high energy physics including

* the standard model (`sm`) up to four-loop accuracy
* the standard model with right-handed neutrinos  (`nusm`) up to two-loop accuracy
* the minimal supersymmetric standard model  (`mssm`) up to three-loop accuracy
* the minimal supersymmetric standard model with right-handed neutrinos  (`numssm`) up to two-loop accuracy
* the Z2 symmetric two Higgs doublet models of type I, II, X, Y  (`thdmi`, `thdmii`, `thdmx`, `thdmy`)  up to three-loop accuracy

In oder to connect the Yukawa matrices with physical observables, **RGE++** also features auxiliary routines for the extraction of CKM and PMNS elements.

**RGE++** internally utilises the Eigen3 library for all algebraic operations of matrix types and the ODEint library for the integration of the RGEs. These functionalities are implemented in base classes. The specific model classes inherit those and serve as containers for the explicit expressions of the RGEs. Thus, additional models can be implemented rather easily.

## Installation

On a UNIX system you first need to install boost, which is best done via the systems package manager.
Then clone the git repository to a local directory on your system by typing

```
    git clone https://github.com/Herren/RGEpp.git
```

In addition to that, you also need the Eigen3 and up-to-date ODEint libraries to be installed on your system. You can obtain them via

```
    git clone https://github.com/boostorg/odeint.git
    git clone https://gitlab.com/libeigen/eigen.git
```

inside **your** local folder named 'rgepp'. This ensures that the include paths in the makefile are correct.

In case you already have installed Eigen3 and/or ODEint on your system or you want to have them in a different folder, adjust the `EIGENPATH` and `ODEINTPATH` variables in `makefile` correspondingly.

## Examples

After installation an easy example can be compiled and run by typing

```
    make sm_example
    ./examples/sm_example 3000 3
```

The command line output should now show the fermion masses and mixing angles of the SM evolved to 3000 GeV at three-loop accuracy.

`numssm_example.cpp` illustrates a simplified SUSY GUT scenario where all Yukawa and gauge couplings are set at the GUT scale. Right-handed neutrinos are consecutively integrated out at their mass scale and the output shows the resulting fermion masses and mixing at the EW scale. The example can be compiled and run by typing

```
    make numssm_example
    ./examples/numssm_example
```

Both code examples are well commented and explained to great detail in the manual.

The makefile calls `clang++` to compile the examples. Should `clang` not be available on your system `g++` also works. Depending on the exact version, it might be necessary
to pass `-ftemplate-depth=1500` as the default value in more recent versions is too low.

## Documentation

**RGE++** comes alongside the paper arXiv:2001.xxxxx which also serves as a manual and is included in the `doc` folder. Just run `pdflatex` (and `bibtex` if you like) for compilation.

The documentation does not only include the examples above but also a detailed description of all model classes and their member functions.

## Scientific credit

If you use **RGE++** as part of your scientific work, please do not only cite arXiv:2001.xxxxx but also give fair credit to the work **RGE++** is built upon. These are

* Eigen3 and ODEint libraries for internal usage
* REAP/MPT package for the routines for extracting fermion observables
* the theorists who calculated the actual RGEs

Table 1 of the manual gives a complete list of publications for each model and auxiliary function implemented in **RGE++**. This list can easily be copied from the `.tex` file and works with inspire. You couldn't have it more convenient :)

## Authors

This project was an essential part of Thomas Deppisch's PhD thesis and first used in JHEP 1901 (2019) 005 Thanks to Stefan Schacht and Martin Spinrath for their extensive testing during this project. For general use **RGE++** has been expanded afterwards and now also features Florian Herren's work on beta functions in the SM and beyond, e.g. Phys.Rev. D97 (2018) no.1, 015016
