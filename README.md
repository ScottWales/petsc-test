Testing PETsc functionality
=========================================

[![Build Status](https://travis-ci.org/ScottWales/petsc-test.png?branch=master)](https://travis-ci.org/ScottWales/petsc-test)

To build:
=========

The build system supports incremental and parallel builds. To compile all
PROGRAMs and tests run:

    $ make

Tests
=====

Tests use [pFunit](http://sourceforge.net/projects/pfunit) & require the PFUNIT environment variable to be set.
Tests are run automatically when building the program, or run

    $ make check

Note that you should use the `pfunit_2.1.0` branch, and it will require
gfortran 4.8 or later (see
http://sourceforge.net/p/pfunit/mailman/message/31811969/)

Automated Travis testing is done using gfortran-4.8, to install on Ubuntu:

    $ sudo apt-add-repository ppa:ubuntu-toolchain-r/test
    $ sudo apt-get update
    $ sudo apt-get install gfortran-4.8
    $ export FC=gfortran-4.8
    $ export OMPI_FC=$FC

