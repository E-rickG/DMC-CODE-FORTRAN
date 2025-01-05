#!/bin/bash

set -e

gfortran -c nr_module.f90
gfortran -c nrutil_module.f90
gfortran -c dmc_module.f90
gfortran -c ran3_module.f90
gfortran -c gasdev_module.f90
gfortran -c potentials_module.f90
gfortran -c montecarlo.f90
echo "Compilation completed!"

