Testing PETsc functionality
=========================================

[![Build Status](https://travis-ci.org/ScottWales/petsc-test.png?branch=master)](https://travis-ci.org/ScottWales/petsc-test)

Vagrant Environment
===================

A [Vagrant](http://www.vagrantup.com) environment is provided to help with development & testing, running `vagrant up` will create a Centos-based VM with gfortran-4.8 and all program pre-requisites running on your local machine. You can then switch to it with `vagrant ssh`.

To build:
=========

The build system supports incremental and parallel builds. To compile all PROGRAMs and tests run:

    $ make

If you are not using the Vagrant environment you can specify the paths to PFunit and PETSc using the `PFUNIT` and `PETSC_DIR` environment variables.

Tests
=====

Tests use [pFunit](http://sourceforge.net/projects/pfunit) & require the PFUNIT environment variable to be set. Tests are run automatically when building the program, or run

    $ make check

Note that you should use the `pfunit_2.1.0` branch, and it will require gfortran 4.8 or later (see http://sourceforge.net/p/pfunit/mailman/message/31811969/)

Automated Travis testing is done using gfortran-4.8, note that all compiler warnings are considered errors in the Travis environment.

