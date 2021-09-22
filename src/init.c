/*
######################################################################
#
# init.c
#
# Written by Carter T. Butts <buttsc@uci.edu>, with portions adapted
# from the Writing R Extensions manual
# Last Modified 09/13/21
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later
#
# Part of the R/relevent package
#
# This file contains routines for registering C entry points.
#
######################################################################
*/

#include <stdlib.h> 
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <R.h>
#include "relevent.h"

#define PA(...) __VA_ARGS__
#define RCALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#define RCDEF(name, n, types) {#name, (DL_FUNC) &name, n, (R_NativePrimitiveArgType []) types}

/*
Here we define every R-callable C function that can be reached via the
.Call interface.  The format here is RCALLDEF(funcname,numargs), where the
function name is unquoted and numargs is the number of arguments.  Note
that if these functions ever change their calling signature, this will need
to be updated as well.  If you need to add more functions, just add a new
line to the table, and the rest should work automagically.
*/
static const R_CallMethodDef R_CallDef[] = {
  RCALLDEF(accum_interact_R, 2),
  RCALLDEF(accum_rrl_R, 2),
  RCALLDEF(accum_ps_R, 1),
  RCALLDEF(acl_ps_R, 3),
  RCALLDEF(acl_tri_R, 2),
  RCALLDEF(drem_n2llik_R, 14),
  RCALLDEF(drem_n2llik_samp_R, 16),
  RCALLDEF(drem_gof_R, 14),
  RCALLDEF(lambda_R, 13),
  {NULL, NULL, 0}
};


/*
Here we define every R-callable C function that can be reached via the
.C interface.  The format here is RCDEF(funcname, numargs, types), where
the arguments are respectively the function name (no quotes), the number of
arguments, and a constant array (i.e., {a, b, c}) containing the R type
codes for the respective arguments.  If you change the calling signature of
any of these functions, you'll need to update this too.  If you need to add
more functions, just add another line to the table, and all will be fixed.
Allegedly.
*/
static const R_CMethodDef R_CDeff[] = {
  RCDEF(rem_int_dev_R, 11, PA({ REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP})),
  RCDEF(rem_int_ev_dev_R, 11, PA({ REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP})),
  RCDEF(rem_ord_dev_R, 11, PA({ REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP})),
  {NULL, NULL, 0, NULL}
};


/*Register the symbols when the package is loaded*/
void R_init_relevent(DllInfo *dll)
{
    R_registerRoutines(dll, R_CDeff, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);  /*This prevents others from being called*/
}

