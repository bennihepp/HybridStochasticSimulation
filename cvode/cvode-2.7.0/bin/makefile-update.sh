#! /bin/bash
# -----------------------------------------------------------------------------------
# $Revision: 1.2 $
# $Date: 2007/12/19 20:33:59 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# This script updates example Makefiles before export.
# It is called by the configure script, after an initial export
# Makefile_ex has been created by config.status.
# -----------------------------------------------------------------

infile="${1}"
solver="${2}"
examples="${3}"
examples_bl="${4}"
solver_lib="${5}"
solver_flib="${6}"

sed "s/@SOLVER@/${solver}/" ${infile}       | \
sed "s/@EXAMPLES@/${examples}/"             | \
sed "s/@EXAMPLES_BL@/${examples_bl}/"       | \
sed "s/@SOLVER_LIB@/${solver_lib}/"         | \
sed "s/@SOLVER_FLIB@/${solver_flib}/"       > foo_makefile

mv foo_makefile ${infile}

