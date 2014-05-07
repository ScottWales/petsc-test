#!/bin/bash
## \file    bootstrap.sh
#  \author  Scott Wales <scott.wales@unimelb.edu.au>
#
#  Copyright 2014 ARC Centre of Excellence for Climate Systems Science
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

# Script for VM setup

yum install --assumeyes gcc gfortran openmpi-devel lapack-devel
module load openmpi-x86_64

# Build PETsc
petsc_ver=3.4.4
cd /tmp
wget -nv http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${petsc_ver}.tar.gz
tar xf petsc-${petsc_ver}.tar.gz
cd /tmp/petsc-${petsc_ver}
./configure
make
make install
