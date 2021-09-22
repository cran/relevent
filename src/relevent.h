/*
######################################################################
#
# relevent.h
#
# Written by Carter T. Butts <buttsc@uci.edu>
# Last Modified 08/29/21
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later
#
# Part of the R/relevent package
#
# This file contains headers for relevent.c.
#
######################################################################
*/
#ifndef RELEVENT_H
  #define RELEVENT_H


/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#define NIDEGSEND  0   /*Norm indegree -> future sending rate*/
#define NIDEGREC   1   /*Norm indegree -> future receiving rate*/
#define NODEGSEND  2   /*Norm outdegree -> future sending rate*/
#define NODEGREC   3   /*Norm outdegree -> future receiving rate*/
#define NTDEGSEND  4   /*Norm total degree -> future sending rate*/
#define NTDEGREC   5   /*Norm total degree -> future receiving rate*/
#define FPSENDSEND 6   /*Fraction past sending -> future sending rate*/
#define FPRECSEND  7   /*Fraction past receipt -> future sending rate*/
#define RRRECSEND  8   /*Recency of receipt -> future sending rate*/
#define RRSENDSEND 9   /*Recency of sending -> future sending rate*/
#define COVSEND    10  /*Covariate effect for sending*/
#define COVREC     11  /*Covariate effect for receiving*/
#define COVSENDREC 12  /*Covariate effect for sending and receiving*/
#define COVEVENT   13  /*Generic event-wise covariate effect*/
#define OTPSEND    14  /*Outbound two-paths -> future sending rate*/
#define ITPSEND    15  /*Incoming two-paths -> future sending rate*/
#define OSPSEND    16  /*Outbound shared partners -> future sending rate*/
#define ISPSEND    17  /*Inbound shared partners -> future sending rate*/
#define FESEND     18  /*Fixed effects for sending*/
#define FEREC      19  /*Fixed effects for receiving*/
#define FESENDREC  20  /*Fixed effects for sending and receiving*/
#define PSABBA     21  /*P-Shift (turn receiving): AB->BA (dyadic)*/
#define PSABB0     22  /*P-Shift (turn receiving): AB->B0 (non-dyadic)*/
#define PSABBY     23  /*P-Shift (turn receiving): AB->BY (dyadic)*/
#define PSA0X0     24  /*P-Shift (turn claiming): A0->X0 (non-dyadic)*/
#define PSA0XA     25  /*P-Shift (turn claiming): A0->XA (non-dyadic)*/
#define PSA0XY     26  /*P-Shift (turn claiming): A0->XY (non-dyadic)*/
#define PSABX0     27  /*P-Shift (turn usurping): AB->X0 (non-dyadic)*/
#define PSABXA     28  /*P-Shift (turn usurping): AB->XA (dyadic)*/
#define PSABXB     29  /*P-Shift (turn usurping): AB->XB (dyadic)*/
#define PSABXY     30  /*P-Shift (turn usurping): AB->XY (dyadic)*/
#define PSA0AY     31  /*P-Shift (turn continuing): A0->AY (non-dyadic)*/
#define PSABA0     32  /*P-Shift (turn continuing): AB->A0 (non-dyadic)*/
#define PSABAY     33  /*P-Shift (turn continuing): AB->AY (dyadic)*/


/*INTERNAL ROUTINES---------------------------------------------------------*/

double acl_adj(SEXP acl, int src, int dest);

double rrl_rank(SEXP rrl, int src, int dest, int mode);

void logrm_irr(SEXP lrm, int n, SEXP rrl, double coef, int mode);

void logrm_irr_samp(SEXP lrv, int ns, int *tail, int *head, int n, SEXP rrl, double coef, int mode);

void logrm_rceff(SEXP lrm, int m, int nvar, int n, int it, int v, double *inparm, double *outparm, double coef, int mode);

void logrm_normint(SEXP lrm, int n, SEXP acl, double *deg, double coef, int mode);

void logrm_rceff_samp(SEXP lrv, int ns, int *tail, int *head, int m, int nvar, int n, int it, int v, double *inparm, double *outparm, double coef, int mode);

void logrm_ladj(SEXP lrm, int n, SEXP adj, double coef, int mode);

void logrm_ladj_samp(SEXP lrv, int ns, int *tail, int *head, SEXP adj, double coef, int mode);

void lambda(SEXP pv, int it, SEXP effects, int n, int m, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrm);

void lambda_samp(SEXP pv, int it, SEXP effects, int nv, int m, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrv, int ns, int *tail, int *head);

int pshiftclassify(int osrc, int odest, int nsrc, int ndest);


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

SEXP accum_interact_R(SEXP elist, SEXP oldacl);

SEXP accum_rrl_R(SEXP elist, SEXP oldrrl);

SEXP accum_ps_R(SEXP elist);

SEXP acl_ps_R(SEXP elist, SEXP n, SEXP oldps);

SEXP acl_tri_R(SEXP acl, SEXP oldtri);

SEXP drem_n2llik_R(SEXP pv, SEXP effects, SEXP edgelist, SEXP n, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrm, SEXP ordinal, SEXP condnum);

SEXP drem_n2llik_samp_R(SEXP pv, SEXP effects, SEXP edgelist, SEXP n, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrv, SEXP tail, SEXP head, SEXP ordinal, SEXP condnum);

SEXP drem_gof_R(SEXP pv, SEXP effects, SEXP edgelist, SEXP n, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrm, SEXP ordinal, SEXP condnum);

void rem_int_dev_R(double *par, int *pnpar, double *evm, int *pm, double *statsa, int *pnet, int *suppm, int *calcderiv, double *val, double *grad, double *hess);

void rem_int_ev_dev_R(double *par, int *pnpar, double *ev, double *statsm, int *pnet, int *suppv, int *calcderiv, double *val, double *grad, double *hess, int *initvals);

void rem_ord_dev_R(double *par, int *pnpar, int *evm, int *pm, double *statsa, int *pnet, int *suppm, int *calcderiv, double *val, double *grad, double *hess);

SEXP lambda_R(SEXP pv, SEXP iter, SEXP effects, SEXP n, SEXP nev, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrm);


#endif
