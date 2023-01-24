/*
######################################################################
#
# relevent.c
#
# Written by Carter T. Butts <buttsc@uci.edu>
# Last Modified 01/24/23
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later
#
# Part of the R/relevent package
#
# This file contains routines for preparing/fitting relational event
# models.
#
######################################################################
*/

#include <stdio.h>
#include <stdlib.h>
#include <Rinternals.h>
#include <R.h>
#include <Rmath.h>
#include "relevent.h"
#include "utils.h"


double acl_adj(SEXP acl, int src, int dest)
/*Given an accumulated communication list (for one iteration only), return the number of calls from src to dest.*/
{
  SEXP srclist,destlist,destval;
  char buf[20];
  double adj;
  
  /*Rprintf("Querying acl_adj for %d->%d\n",src+1,dest+1);*/
  /*Check to see if the source is present*/
  snprintf(buf,sizeof(buf),"%d",src+1);
  if((srclist=getListElement(acl,buf))==R_NilValue)
    return 0.0;

  /*If still here, check to see if the destination is present*/
  snprintf(buf,sizeof(buf),"%d",dest+1);
  if((destlist=getListElement(srclist,buf))==R_NilValue)
    return 0.0;

  /*Finally, if we're _still_ here, return the appropriate value*/
  PROTECT(destval=coerceVector(destlist,REALSXP));
  adj=REAL(destval)[0];
  /*Rprintf("\tFound value of %.0f\n",adj);*/

  UNPROTECT(1);
  return adj;
}


double rrl_rank(SEXP rrl, int src, int dest, int mode)
/*Given a recency ranked communication list (for one iteration only), return the rank of the call from src to dest in terms of src's outgoing calls (mode=0), or dest's incoming calls (mode=1).  Presumably rrl has already been subsetted in the appropriate way for this to make sense.*/
{
  SEXP srclist,destlist,names;
  char buf[20];
  double count;
  
  if(mode==0){ /*First look for source, then see how recent dest is*/
    /*Check to see if the source is present*/
    snprintf(buf,sizeof(buf),"%d",src+1);
    if((srclist=getListElement(rrl,buf))==R_NilValue)
      return DBL_MAX;

    /*If still here, check to see if the destination is present*/
    PROTECT(names=coerceVector(srclist,INTSXP));  /*Get name vector*/
    for(count=0.0;(count<(double)length(names))&& (INTEGER(names)[(int)count]!=dest+1); count++);
    if(INTEGER(names)[(int)count]==dest+1){
      UNPROTECT(1);
      return count+1.0;        /*If we found it, return the rank*/
    }else{
      UNPROTECT(1);
      return DBL_MAX;          /*Otherwise, return something huge*/
    }
  }else{       /*First look for dest, then see how recent source is*/
    /*Check to see if the destination is present*/
    snprintf(buf,sizeof(buf),"%d",dest+1);
    if((destlist=getListElement(rrl,buf))==R_NilValue)
      return DBL_MAX;

    /*If still here, check to see if the source is present*/
    PROTECT(names=coerceVector(destlist,INTSXP));  /*Get name vector*/
    for(count=0.0;(count<(double)length(names))&& (INTEGER(names)[(int)count]!=src+1); count++);
    if(INTEGER(names)[(int)count]==src+1){
      UNPROTECT(1);
      return count+1.0;        /*If we found it, return the rank*/
    }else{
      UNPROTECT(1);
      return DBL_MAX;          /*Otherwise, return something huge*/
    }
  }
}


SEXP accum_interact_R(SEXP elist, SEXP oldacl)
/*Given a list of events (in (time,src,dest) form), return a list of vectors with the accumulated number of interactions for each actor at each event.  If the second argument, oldacl, is non-null, then oldacl is expanded to add any unprocessed events in elist; if elist is shorter than oldacl, an error is thrown.  While a new list object is returned, its contents are mostly ported over from the old list, so the old list should not then be accessed directly unless you know what you are doing (since it will point to the same memory locations).*/
{
  int pc=0,m,i,om;
  SEXP acl,list,srcl,destl,elem;

  /*Check to be sure the inputs are valid*/
  if(oldacl!=R_NilValue)
    om=0;
  else
    om=length(oldacl);
  if(nrows(elist)<=om)
    error("Passed an edgelist to accum_interact_R that is shorter than the old acl it was intended to update!  Don't do that.\n");

  /*Allocate memory for the (new acl)*/
  m=nrows(elist);
  PROTECT(elist=coerceVector(elist,STRSXP)); pc++;
  PROTECT(acl=allocVector(VECSXP,m)); pc++;  

  /*Create the acl by iteration*/
  PROTECT(elem=allocVector(VECSXP,0)); pc++;
  SET_VECTOR_ELT(acl,0,elem);
  for(i=1;i<m;i++){
    if(i<=om){              /*Transfer over the old acl*/
      SET_VECTOR_ELT(acl,i-1,VECTOR_ELT(oldacl,i-1));
    }else{                  /*Now, expand*/
      PROTECT(list=duplicate(VECTOR_ELT(acl,i-1))); pc++;  /*Copy the old list*/
      srcl=getListElement(list,CHAR(STRING_ELT(elist,i-1+m))); /*Try to get src*/
      if(srcl==R_NilValue){               /*If no src, create a new list*/
        PROTECT(srcl=allocVector(VECSXP,0)); pc++;
        PROTECT(elem=allocVector(INTSXP,1)); pc++;
        INTEGER(elem)[0]=1;
        PROTECT(srcl=setListElement(srcl,CHAR(STRING_ELT(elist,i-1+2*m)),elem)); pc++;
        PROTECT(list=setListElement(list,CHAR(STRING_ELT(elist,i-1+m)),srcl)); pc++;
      }else{                             /*Otherwise, set dest*/
        destl=getListElement(srcl,CHAR(STRING_ELT(elist,i-1+2*m)));  /*Get dest*/
        if(destl==R_NilValue){
          PROTECT(elem=allocVector(INTSXP,1)); pc++;
          INTEGER(elem)[0]=1;
          PROTECT(srcl=setListElement(srcl,CHAR(STRING_ELT(elist,i-1+2*m)),elem)); pc++;
          list=setListElement(list,CHAR(STRING_ELT(elist,i-1+m)),srcl);
        }else{
          PROTECT(elem=coerceVector(destl,INTSXP)); pc++;
          INTEGER(elem)[0]++;
          srcl=setListElement(srcl,CHAR(STRING_ELT(elist,i-1+2*m)),elem);
        }
      }
      SET_VECTOR_ELT(acl,i,list);   /*Install the updated list*/
      if(pc>1000){             /*If the stack is getting out of hand, fix it*/
        UNPROTECT(pc-3); pc=3;
      }
    }
  }

  /*Unprotect and return*/
  UNPROTECT(pc);
  return acl;
}


SEXP accum_ps_R(SEXP elist)
/*Build matrix of accumulated participation shift (P-shift) counts.  Output is 
in the form of a matrix whose rows contain counts for each of the 13 P-shift 
types, such that the ith row contains the counts immediately prior to the
resolution of event i.  Thus, there are m+1 rows in all, where m is the number
of events.  The P-shifts are given in the order used in Gibson's 2003 Social 
Forces paper, namely:

  (Turn Receiving)  [0] AB->BA, [1] AB->B0, [2] AB->BY,
  (Turn Claiming)   [3] A0->X0, [4] A0->XA, [5] A0->XY,
  (Turn Usurping)   [6] AB->X0, [7] AB->XA, [8] AB->XB, [9] AB->XY,
  (Turn Continuing) [10] A0->AY, [11] AB->A0, [12] AB->AY  

(This uses Gibson's notation, in which A is the initial source, B is the
initial target, X is a new (shifted) speaker, Y is a new (shifted) target,
and 0 is used where no well-defined speaker or target is present.  Here, this
occurs when NA is given for source or destination.)

It is worth noting that not all adjacent event pairs induce P-shifts, and hence
the shift counts will not increment with every event.  In particular, the first
event does not induce a shift (since there is no prior event), and neither does
a repetition of a previous event (e.g., AB->AB or A0->A0).  The full set is
thus affinely independent in general, although they will have a near (or 
even full) dimension of affine dependence on most data sets.  You probably
don't want to use all the stats at once, although we compute them all (since
the cost of doing so is trivial).
*/
{
  int i,j,m,pc=0,src,dest,osrc,odest;
  SEXP psmat,elem;
  double *psm;
  
  /*Set things up*/
  m=nrows(elist);
  PROTECT(elist=coerceVector(elist,STRSXP)); pc++;
  PROTECT(psmat=allocVector(REALSXP,(m+1)*13)); pc++;
  psm=REAL(psmat);
  for(i=0;i<2;i++)
    for(j=0;j<13;j++)
      psm[i+j*(m+1)]=0.0;

  /*Get first source/dest pair*/
  PROTECT(elem=allocVector(STRSXP,1));
  SET_STRING_ELT(elem,0,STRING_ELT(elist,m));
  PROTECT(elem=coerceVector(elem,INTSXP));
  src=asInteger(elem);   
  PROTECT(elem=allocVector(STRSXP,1));
  SET_STRING_ELT(elem,0,STRING_ELT(elist,2*m));
  PROTECT(elem=coerceVector(elem,INTSXP));  
  dest=asInteger(elem);
  UNPROTECT(4);

  /*Move sequentially through the event list, accumulating shifts*/
  /*Rprintf("Processing event list.\n");*/
  for(i=2;i<m+1;i++){  /*Note: psm[0:1,]==0.0, and we go out to m*/
    /*Rprintf("\tIteration %d\n",i);*/
    /*Begin with the old counts*/
    for(j=0;j<13;j++)
      psm[i+j*(m+1)]=psm[(i-1)+j*(m+1)];
    /*Save the old source/dest pair*/
    osrc=src;
    odest=dest;
    /*Get the current source/dest pair*/
    /*Rprintf("\t\tGetting new src/dest...\n");*/
    PROTECT(elem=allocVector(STRSXP,1));                /*Possibly the worst*/
    SET_STRING_ELT(elem,0,STRING_ELT(elist,i-1+m));     /*String->Int*/
    PROTECT(elem=coerceVector(elem,INTSXP));            /*conversion ever*/
    src=asInteger(elem);                                /*written.  I blame*/
    PROTECT(elem=allocVector(STRSXP,1));                /*society.*/
    SET_STRING_ELT(elem,0,STRING_ELT(elist,i-1+2*m));
    PROTECT(elem=coerceVector(elem,INTSXP));  
    dest=asInteger(elem);
    UNPROTECT(4);
    /*Rprintf("\t\t\tsrc=%d, dest=%d\n",src,dest);*/
    /*Now, classify the shift*/
    if((src!=NA_INTEGER)&&(osrc!=NA_INTEGER)&&((src!=osrc)||(dest!=odest))){  /*First, verify that there was a shift*/
      if(odest==NA_INTEGER){                /*A0->??*/
        if(dest==NA_INTEGER){                   /*A0->X0*/
          psm[i+3*(m+1)]++;
        }else if(dest==osrc){                   /*A0->XA*/
          psm[i+4*(m+1)]++;
        }else if(src!=osrc){                    /*A0->XY*/
          psm[i+5*(m+1)]++;
        }else{                                  /*A0->AY*/
          psm[i+10*(m+1)]++;
        }
      }else{                                /*AB->??*/
        if(src==osrc){                        /*AB->A?*/
          if(dest==NA_INTEGER){                 /*AB->A0*/
            psm[i+11*(m+1)]++;
          }else{                                /*AB->AY*/
            psm[i+12*(m+1)]++;
          }
        }else if(src==odest){                 /*AB->B?*/
          if(dest==osrc){                       /*AB->BA*/
            psm[i+0*(m+1)]++;
          }else if(dest==NA_INTEGER){           /*AB->B0*/
            psm[i+1*(m+1)]++;
          }else{                                /*AB->BY*/
            psm[i+2*(m+1)]++;
          }
        }else{                                /*AB->X?*/
          if(dest==NA_INTEGER){                 /*AB->X0*/
            psm[i+6*(m+1)]++;
          }else if(dest==osrc){                 /*AB->XA*/
            psm[i+7*(m+1)]++;
          }else if(dest==odest){                /*AB->XB*/
            psm[i+8*(m+1)]++;
          }else{                                /*AB->XY*/
            psm[i+9*(m+1)]++;
          }
        }
      }
    }
  }

  /*Unprotect and return*/
  /*Rprintf("Done!\n");*/
  UNPROTECT(pc);
  return psmat;
}


int pshiftclassify(int osrc, int odest, int nsrc, int ndest)
/*Given interactions (osrc,odest)->(nsrc,ndest), identify the corresponding
P-shift using Gibson's 2003 classification.  Shifts are returned as integers
matching the order of the shift from 0 to 12, as per the scheme used elsewhere
in this code; this is also the order he uses in his paper.  If a shift is not
present (or not well-defined), a negative value is returned.*/
{
  if((osrc==NA_INTEGER)||(nsrc==NA_INTEGER))   /*Undefined case*/
    return -1;
  if((osrc==nsrc)&&(odest==ndest))             /*No change->no shift*/
    return -1;
  if(odest==NA_INTEGER){                /*A0->??*/
    if(ndest==NA_INTEGER){                   /*A0->X0*/
      return 3;
    }else if(ndest==osrc){                   /*A0->XA*/
      return 4;
    }else if(nsrc!=osrc){                    /*A0->XY*/
      return 5;
    }else{                                   /*A0->AY*/
      return 10;
    }
  }else{                                /*AB->??*/
    if(nsrc==osrc){                        /*AB->A?*/
      if(ndest==NA_INTEGER){                 /*AB->A0*/
        return 11;
      }else{                                 /*AB->AY*/
        return 12;
      }
    }else if(nsrc==odest){                 /*AB->B?*/
      if(ndest==osrc){                        /*AB->BA*/
        return 0;
      }else if(ndest==NA_INTEGER){            /*AB->B0*/
        return 1;
      }else{                                  /*AB->BY*/
        return 2;
      }
    }else{                                /*AB->X?*/
      if(ndest==NA_INTEGER){                  /*AB->X0*/
        return 6;
      }else if(ndest==osrc){                  /*AB->XA*/
        return 7;
      }else if(ndest==odest){                 /*AB->XB*/
        return 8;
      }else{                                  /*AB->XY*/
        return 9;
      }
    }
  }
  /*We should never get here!!!*/
  return -2;
}


SEXP acl_ps_R(SEXP elist, SEXP n, SEXP oldps)
/*Build a list of P-shift changescores over time.  The output is in the form
of an acl, namely a list with the structure
  shifttype
    $iter
      $ego
        $alter
          $count (always 1)
Shift types are in uppercase, four-byte form, following Gibson's notation (see
also above).  Where alter does not appear in ego's list for a given 
shift type/iteration combination, this means that an event from ego to alter
in that iteration would not result in an additional P-shift of the indicated
type.  (There are, hence, no ties for the first event, since it cannot create
P-shifts.)  Otherwise, this routine has the same basic assumptions as 
accum_ps_R.

If oldps is given, it is assumed to contain a structure of the identical type
to the output of this function, for a subset of elist; in this case, oldps
is then extended to include the new events.  (Obviously, elist needs to include
everything in oldps.)
*/
{
  int i,j,k,m,om=0,pc=0,src,dest,in,pst,pctemp,pctemp2;
  SEXP psl,elem,elem2,lptrs[13];
  char ss[50],sd[50],*snames[13] = {"ABBA","ABB0","ABBY","A0X0","A0XA","A0XY", "ABX0","ABXA","ABXB","ABXY","A0AY","ABA0","ABAY"};
  
  /*Set things up*/
  m=nrows(elist);
  PROTECT(n=coerceVector(n,INTSXP)); pc++;
  in=INTEGER(n)[0];
  PROTECT(elist=coerceVector(elist,STRSXP)); pc++;
  PROTECT(psl=allocVector(VECSXP,0)); pc++;
  for(i=0;i<13;i++){        /*Build the PS list*/
    PROTECT(elem=allocVector(VECSXP,m)); pc++;
    PROTECT(psl=setListElement(psl,snames[i],elem)); pc++;
  }

  /*Initialize/move over existing stuff*/
  if(oldps!=R_NilValue){
    om=length(VECTOR_ELT(oldps,0));
    if(m<om)
     error("New elist length (%d) was shorter than old one (%d) in acl_ps_R.  Don't do that.\n",m,om);
    for(i=0;i<13;i++)
      for(j=0;j<om;j++){
        SET_VECTOR_ELT(VECTOR_ELT(psl,i),j,VECTOR_ELT(VECTOR_ELT(oldps,i),j));
      }
  }else{
    /*Write empty lists in place for the first iteration (no shifts)*/
    for(i=0;i<13;i++){
      PROTECT(elem=allocVector(VECSXP,0)); pc++;
      SET_VECTOR_ELT(VECTOR_ELT(psl,i),0,elem);
    }  
    om=1;
  }

  /*Walk through the event list, accumulating changescores as we go*/
  for(i=om;i<m;i++){
//    Rprintf("Event %d\n",i);
    /*Get source/dest pair from the last event*/
    PROTECT(elem=allocVector(STRSXP,1));
    SET_STRING_ELT(elem,0,STRING_ELT(elist,i-1+m));
    PROTECT(elem=coerceVector(elem,INTSXP));
    src=asInteger(elem);   
    PROTECT(elem=allocVector(STRSXP,1));
    SET_STRING_ELT(elem,0,STRING_ELT(elist,i-1+2*m));
    PROTECT(elem=coerceVector(elem,INTSXP));  
    dest=asInteger(elem);
    UNPROTECT(4);
//    Rprintf("\t%d -> %d\n",src,dest);
    /*Create the lists for this iteration*/
    pctemp=pc;
    for(j=0;j<13;j++){
      PROTECT(lptrs[j]=allocVector(VECSXP,0)); pc++;
    }
    /*Consider all possible moves, and classify accordingly*/
    for(j=1;j<=in;j++){
      for(k=1;k<=in;k++)
        if(j!=k){              /*We assume no self-ties here*/
          pst=pshiftclassify(src,dest,j,k);  /*What shift is this?*/
          if(pst>=0){    /*If a shift would occur, add it to the list*/
            pctemp2=pc;
            snprintf(ss,50,"%d",j);   /*Create source string*/
            snprintf(sd,50,"%d",k);   /*Create dest string*/
            PROTECT(elem=allocVector(REALSXP,1)); pc++;  /*Create new entry*/
            REAL(elem)[0]=1.0;
            elem2=getListElement(lptrs[pst],ss);  /*Get src list*/
            if(elem2==R_NilValue){      /*If src not found, create*/
              PROTECT(elem2=allocVector(VECSXP,0)); pc++;
            }
            PROTECT(elem2=setListElement(elem2,sd,elem)); pc++;  /*set entry*/
            PROTECT(lptrs[pst]=setListElement(lptrs[pst],ss,elem2)); pc++;
            SET_VECTOR_ELT(VECTOR_ELT(psl,pst),i,lptrs[pst]);
            UNPROTECT(pc-pctemp2);
            pc=pctemp2;
          }  
        }
    }
    /*Write the lists into place*/
    for(j=0;j<13;j++)
      SET_VECTOR_ELT(VECTOR_ELT(psl,j),i,lptrs[j]);
    /*Unprotect the temp lists, now that they are assigned*/
    UNPROTECT(pc-pctemp);
    pc=pctemp;
  }

  /*Unprotect and return*/
//  Rprintf("Done!\n");
  UNPROTECT(pc);
  return psl;
}


SEXP accum_rrl_R(SEXP elist, SEXP oldrrl)
/*Build list of accumulated incoming and outgoing recency-ranked communications.
Specifically, the output is a nested list whose first dimension is in vs out,
second dimension is iteration, third dimension is vertex, and fourth
dimension (where applicable) is a recency ordered vector of communication
partners (most recent to least recent).  This structure is built from
the edgelist matrix, rather than the acl, due to the fact that this task
is much easier to perform for the former than the latter.

If oldrrl is given, then entries from this are used before adding entries for
newer time points; elist should include the history indexed by oldrrl, plus
any new time steps.*/
{
  int i,m,om=0,pc=0,src,dest;
  SEXP ircl,orcl,rrl,elem,elem2,il,ol;
  
  /*Perform initial setup*/
  m=nrows(elist);
  PROTECT(elist=coerceVector(elist,STRSXP)); pc++;
  PROTECT(ircl=allocVector(VECSXP,m)); pc++;
  PROTECT(orcl=allocVector(VECSXP,m)); pc++;
  if(oldrrl!=R_NilValue){                    /*Transfer old list elements*/
    elem=getListElement(oldrrl,"in");
    if(length(ircl)<length(elem))
      error("New elist shorter than old one....\n");
    for(i=0;i<length(elem);i++)
      SET_VECTOR_ELT(ircl,i,VECTOR_ELT(elem,i));
    elem=getListElement(oldrrl,"out");
    if(length(orcl)<length(elem))
      error("New elist shorter than old one....\n");
    for(i=0;i<length(elem);i++)
      SET_VECTOR_ELT(orcl,i,VECTOR_ELT(elem,i));
    om=length(elem);
  }else{
    om=0;
  }

  /*Build incoming/outgoing recency-ordered contact lists*/
  if(om==0){
    PROTECT(elem=allocVector(VECSXP,0)); pc++;
    SET_VECTOR_ELT(ircl,0,elem);
    PROTECT(elem=allocVector(VECSXP,0)); pc++;
    SET_VECTOR_ELT(orcl,0,elem);
    om++;
  }
  for(i=om;i<m;i++){
    PROTECT(il=duplicate(VECTOR_ELT(ircl,i-1))); pc++;  /*Copy the old lists*/
    PROTECT(ol=duplicate(VECTOR_ELT(orcl,i-1))); pc++;
    PROTECT(elem=allocVector(STRSXP,1)); pc++;          /*Possibly the worst*/
    SET_STRING_ELT(elem,0,STRING_ELT(elist,i-1+m));     /*String->Int*/
    PROTECT(elem=coerceVector(elem,INTSXP)); pc++;      /*conversion ever*/
    src=asInteger(elem);                                /*written.  I blame*/
    PROTECT(elem=allocVector(STRSXP,1)); pc++;          /*society.*/
    SET_STRING_ELT(elem,0,STRING_ELT(elist,i-1+2*m));
    PROTECT(elem=coerceVector(elem,INTSXP)); pc++;
    dest=asInteger(elem);
    /*Update the outgoing comm structure*/
    elem=getListElement(ol,CHAR(STRING_ELT(elist,i-1+m)));
    if(length(elem)==0){
      PROTECT(elem=allocVector(INTSXP,1)); pc++;
      INTEGER(elem)[0]=dest;
    }else{
      PROTECT(elem2=vecRemove(elem,(double)dest)); pc++;
      PROTECT(elem2=coerceVector(elem2,INTSXP)); pc++;
      PROTECT(elem=allocVector(INTSXP,1)); pc++;
      INTEGER(elem)[0]=dest;
      PROTECT(elem=vecAppend(elem,elem2)); pc++;
    }
    PROTECT(ol=setListElement(ol,CHAR(STRING_ELT(elist,i-1+m)),elem)); pc++;
    /*Update the incoming comm structure*/
    elem=getListElement(il,CHAR(STRING_ELT(elist,i-1+2*m)));
    if(length(elem)==0){
      PROTECT(elem=allocVector(INTSXP,1)); pc++;
      INTEGER(elem)[0]=src;
    }else{
      PROTECT(elem2=vecRemove(elem,(double)src)); pc++;
      PROTECT(elem2=coerceVector(elem2,INTSXP)); pc++;
      PROTECT(elem=allocVector(INTSXP,1)); pc++;
      INTEGER(elem)[0]=src;
      PROTECT(elem=vecAppend(elem,elem2)); pc++;
    }
    PROTECT(il=setListElement(il,CHAR(STRING_ELT(elist,i-1+2*m)),elem)); pc++;
    /*Write updated lists in place*/
    SET_VECTOR_ELT(ircl,i,il);
    SET_VECTOR_ELT(orcl,i,ol);
    if(pc>1000){             /*If the stack is getting out of hand, fix it*/
      UNPROTECT(pc-5); pc=5;
    }
  }

  /*Create the rrl list*/
  PROTECT(rrl=allocVector(VECSXP,0)); pc++;
  PROTECT(rrl=setListElement(rrl,"in",ircl)); pc++;
  PROTECT(rrl=setListElement(rrl,"out",orcl)); pc++;

  /*Unprotect and return*/
  UNPROTECT(pc);
  return rrl;
}


SEXP acl_tri_R(SEXP acl, SEXP oldtri)
/*Return a list of pseudo-adjacencies reflecting triadic (really, 2-path) 
properties:
 $ "top"|"tip"|"sop"|"sip"
                          $iter
                                $ ego
                                      $ alter
                                              $ count

If oldtri is given, then the triangle structure contained is extended for the
length of acl.  This allows a triadic structure to be updated without completely
recomputing it.
*/
{
  int i,j,k,l,g,pc=0,pc2=0,om=0;
  SEXP tri,top,tip,sop,sip,list,alt,lnam,altnam,alt2,alt2nam,elem,elem2,elem3;

  /*Allocate memory for the triad lists*/
  PROTECT(top=allocVector(VECSXP,length(acl))); pc++;
  PROTECT(tip=allocVector(VECSXP,length(acl))); pc++;
  PROTECT(sop=allocVector(VECSXP,length(acl))); pc++;
  PROTECT(sip=allocVector(VECSXP,length(acl))); pc++;
  if(oldtri!=R_NilValue){                                /*If there is an old object, move list elements over*/
    alt=getListElement(oldtri,"top");
    if(length(top)<length(alt))
      error("New acl shorter than old one....\n");
    for(i=0;i<length(alt);i++)
      SET_VECTOR_ELT(top,i,VECTOR_ELT(alt,i));
    alt=getListElement(oldtri,"tip");
    if(length(tip)<length(alt))
      error("New acl shorter than old one....\n");
    for(i=0;i<length(alt);i++)
      SET_VECTOR_ELT(tip,i,VECTOR_ELT(alt,i));
    alt=getListElement(oldtri,"sop");
    if(length(sop)<length(alt))
      error("New acl shorter than old one....\n");
    for(i=0;i<length(alt);i++)
      SET_VECTOR_ELT(sop,i,VECTOR_ELT(alt,i));
    alt=getListElement(oldtri,"sip");
    if(length(sip)<length(alt))
      error("New acl shorter than old one....\n");
    for(i=0;i<length(alt);i++)
      SET_VECTOR_ELT(sip,i,VECTOR_ELT(alt,i));
    om=length(alt);
  }else{                                           /*Otherwise, use new lists*/
    om=0;
  }

  /*Walk through the (new) iterations*/
  for(i=om;i<length(acl);i++){
    list=VECTOR_ELT(acl,i);
    lnam=getAttrib(list,R_NamesSymbol);
    /*Get the inbound/outbound two-paths*/
    if(length(list)>1){             /*Must have two verts w/outedges*/
      for(j=0;j<length(list);j++){   /*For each ego...*/
        alt=VECTOR_ELT(list,j);
        altnam=getAttrib(alt,R_NamesSymbol);
        for(k=0;k<length(alt);k++){  /*...check each alter*/
          alt2=getListElement(list,CHAR(STRING_ELT(altnam,k)));  /*Alt's alts*/
          alt2nam=getAttrib(alt2,R_NamesSymbol);
          for(l=0;l<length(alt2);l++){
            /*Each element of alt2 != list[j] contributes a 2-path*/
            if(strcmp(CHAR(STRING_ELT(lnam,j)),CHAR(STRING_ELT(alt2nam,l)))!=0){
              /*Handle outgoing two-paths*/
              elem=VECTOR_ELT(top,i);
              pc2=0;
              if(length(elem)==0){        /*Create from scratch*/
                PROTECT(elem=allocVector(VECSXP,0)); pc2++;   /*iter*/
                PROTECT(elem2=allocVector(VECSXP,0)); pc2++;  /*src*/
                PROTECT(elem3=allocVector(INTSXP,1)); pc2++;  /*dest*/
                INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt2,l)));
                PROTECT(elem2=setListElement(elem2,CHAR(STRING_ELT(alt2nam,l)), elem3)); pc2++;
                PROTECT(elem=setListElement(elem,CHAR(STRING_ELT(lnam,j)), elem2)); pc2++;
              }else{                     /*Add to existing*/
                elem2=getListElement(elem,CHAR(STRING_ELT(lnam,j)));
                if(length(elem2)==0){     /*Need to create src's list*/
                  PROTECT(elem2=allocVector(VECSXP,0)); pc2++;  /*src*/
                  PROTECT(elem3=allocVector(INTSXP,1)); pc2++;  /*dest*/
                  INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt2,l)));
                  PROTECT(elem2=setListElement(elem2, CHAR(STRING_ELT(alt2nam,l)), elem3)); pc2++;
                  PROTECT(elem=setListElement(elem,CHAR(STRING_ELT(lnam,j)), elem2)); pc2++;
                }else{                    /*Need to set dest w/in src*/
                  elem3=getListElement(elem2,CHAR(STRING_ELT(alt2nam,l)));
                  if(length(elem3)==0){
                    PROTECT(elem3=allocVector(INTSXP,1)); pc2++;
                    INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt2,l)));
                  }else{
                    PROTECT(elem3=coerceVector(elem3,INTSXP)); pc2++;
                    INTEGER(elem3)[0]+=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt2,l)));
                  }
                  PROTECT(elem2=setListElement(elem2, CHAR(STRING_ELT(alt2nam,l)), elem3)); pc2++;
                  PROTECT(elem=setListElement(elem,CHAR(STRING_ELT(lnam,j)), elem2)); pc2++;
                }
              }
              SET_VECTOR_ELT(top,i,elem);   /*Write list into place*/
              UNPROTECT(pc2);
              /*Handle incoming two-paths (same, but dest/src swap)*/
              elem=VECTOR_ELT(tip,i);
              pc2=0;
              if(length(elem)==0){        /*Create from scratch*/
                PROTECT(elem=allocVector(VECSXP,0)); pc2++;   /*iter*/
                PROTECT(elem2=allocVector(VECSXP,0)); pc2++;  /*dest*/
                PROTECT(elem3=allocVector(INTSXP,1)); pc2++;  /*src*/
                INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt2,l)));
                PROTECT(elem2=setListElement(elem2,CHAR(STRING_ELT(lnam,j)), elem3)); pc2++;
                PROTECT(elem=setListElement(elem,CHAR(STRING_ELT(alt2nam,l)), elem2)); pc2++;
              }else{                     /*Add to existing*/
                elem2=getListElement(elem,CHAR(STRING_ELT(alt2nam,l)));
                if(length(elem2)==0){     /*Need to create src's list*/
                  PROTECT(elem2=allocVector(VECSXP,0)); pc2++;  /*dest*/
                  PROTECT(elem3=allocVector(INTSXP,1)); pc2++;  /*src*/
                  INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt2,l)));
                  PROTECT(elem2=setListElement(elem2, CHAR(STRING_ELT(lnam,j)), elem3)); pc2++;
                  PROTECT(elem=setListElement(elem,CHAR(STRING_ELT(alt2nam,l)), elem2)); pc2++;
                }else{                    /*Need to set src w/in dest*/
                  elem3=getListElement(elem2,CHAR(STRING_ELT(lnam,j)));
                  if(length(elem3)==0){
                    PROTECT(elem3=allocVector(INTSXP,1)); pc2++;
                    INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt2,l)));
                  }else{
                    PROTECT(elem3=coerceVector(elem3,INTSXP)); pc2++;
                    INTEGER(elem3)[0]+=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt2,l)));
                  }
                  PROTECT(elem2=setListElement(elem2, CHAR(STRING_ELT(lnam,j)), elem3)); pc2++;
                  PROTECT(elem=setListElement(elem,CHAR(STRING_ELT(alt2nam,l)), elem2)); pc2++;
                }
              }
              SET_VECTOR_ELT(tip,i,elem);   /*Write list into place*/
              UNPROTECT(pc2);
            }
          }
        }
      }
    }
    if(pc>1000){             /*If the stack is getting out of hand, fix it*/
      UNPROTECT(pc-4); pc=4;
    }
    /*Get the outbound shared partners*/
    lnam=getAttrib(list,R_NamesSymbol);
    for(j=0;j<length(list);j++){    /*Walk each pair with out-edges*/
      for(k=j+1;k<length(list);k++){
        alt=VECTOR_ELT(list,j);
        altnam=getAttrib(alt,R_NamesSymbol);
        alt2=VECTOR_ELT(list,k);
        alt2nam=getAttrib(alt2,R_NamesSymbol);
        for(l=0;l<length(alt);l++){    /*Each intersection is a sp*/
          for(g=0;g<length(alt2);g++){
            if(strcmp(CHAR(STRING_ELT(altnam,l)), CHAR(STRING_ELT(alt2nam,g)))==0){
              elem=VECTOR_ELT(sop,i);      /*Get sop[[iter]]*/
              pc2=0;
              /*Set j->k*/
              if(length(elem)==0){
                PROTECT(elem=allocVector(VECSXP,0)); pc2++;
                PROTECT(elem2=allocVector(VECSXP,0)); pc2++;
                PROTECT(elem3=allocVector(INTSXP,1)); pc2++;
                INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,l)), asInteger(VECTOR_ELT(alt2,g)));
                PROTECT(elem2=setListElement(elem2,CHAR(STRING_ELT(lnam,k)), elem3)); pc2++;
                PROTECT(elem=setListElement(elem,CHAR(STRING_ELT(lnam,j)), elem2)); pc2++;
              }else{
                elem2=getListElement(elem,CHAR(STRING_ELT(lnam,j)));
                if(length(elem2)==0){
                  PROTECT(elem2=allocVector(VECSXP,0)); pc2++;
                  PROTECT(elem3=allocVector(INTSXP,1)); pc2++;
                  INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,l)), asInteger(VECTOR_ELT(alt2,g)));
                  PROTECT(elem2=setListElement(elem2, CHAR(STRING_ELT(lnam,k)), elem3)); pc2++;
                }else{
                  elem3=getListElement(elem2,CHAR(STRING_ELT(lnam,k)));
                  if(length(elem3)==0){
                    PROTECT(elem3=allocVector(INTSXP,1)); pc2++;
                    INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,l)), asInteger(VECTOR_ELT(alt2,g)));
                  }else{
                    PROTECT(elem3=coerceVector(elem3,INTSXP)); pc2++;
                    INTEGER(elem3)[0]+=MIN(asInteger(VECTOR_ELT(alt,l)), asInteger(VECTOR_ELT(alt2,g)));
                  }
                  PROTECT(elem2=setListElement(elem2, CHAR(STRING_ELT(lnam,k)), elem3)); pc2++;
                }
                PROTECT(elem=setListElement(elem,CHAR(STRING_ELT(lnam,j)), elem2)); pc2++;
              }
              /*Set k->j (note that elem must exist by now!)*/
              elem2=getListElement(elem,CHAR(STRING_ELT(lnam,j)));
              if(length(elem2)==0){
                PROTECT(elem2=allocVector(VECSXP,0)); pc2++;
                PROTECT(elem3=allocVector(INTSXP,1)); pc2++;
                INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,l)), asInteger(VECTOR_ELT(alt2,g)));
                PROTECT(elem2=setListElement(elem2, CHAR(STRING_ELT(lnam,k)), elem3)); pc2++;
              }else{
                elem3=getListElement(elem2,CHAR(STRING_ELT(lnam,k)));
                if(length(elem3)==0){
                  PROTECT(elem3=allocVector(INTSXP,1)); pc2++;
                  INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,l)), asInteger(VECTOR_ELT(alt2,g)));
                }else{
                  PROTECT(elem3=coerceVector(elem3,INTSXP)); pc2++;
                  INTEGER(elem3)[0]+=MIN(asInteger(VECTOR_ELT(alt,l)), asInteger(VECTOR_ELT(alt2,g)));
                }
                PROTECT(elem2=setListElement(elem2, CHAR(STRING_ELT(lnam,k)), elem3)); pc2++;
              }
              PROTECT(elem=setListElement(elem,CHAR(STRING_ELT(lnam,j)), elem2)); pc2++;
              SET_VECTOR_ELT(sop,i,elem);  /*Write list into place*/
              UNPROTECT(pc2);
            }
          }
        }
      }
    }    
    if(pc>1000){             /*If the stack is getting out of hand, fix it*/
      UNPROTECT(pc-4); pc=4;
    }
    /*Get the inbound shared partners*/
    for(j=0;j<length(list);j++){
      alt=VECTOR_ELT(list,j);
      altnam=getAttrib(alt,R_NamesSymbol);
      if(length(alt)>1){           /*list[j] adds 1 sp to all pairs in alt*/
        for(k=0;k<length(alt);k++){
          for(l=0;l<length(alt);l++){
            if(k!=l){
              elem=VECTOR_ELT(sip,i);      /*Get sip[[iter]]*/
              pc2=0;
              if(length(elem)==0){         /*Have to create all from scratch*/
                PROTECT(elem=allocVector(VECSXP,0)); pc2++;
                PROTECT(elem2=allocVector(VECSXP,0)); pc2++;
                PROTECT(elem3=allocVector(INTSXP,1)); pc2++;
                INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt,l)));
                PROTECT(elem2=setListElement(elem2,CHAR(STRING_ELT(altnam,l)), elem3)); pc2++;
                PROTECT(elem=setListElement(elem,CHAR(STRING_ELT(altnam,k)), elem2)); pc2++;
              }else{                       /*Add to sip[[iter]]*/
                elem2=getListElement(elem,CHAR(STRING_ELT(altnam,k)));
                if(length(elem2)==0){
                  PROTECT(elem2=allocVector(VECSXP,0)); pc2++;
                  PROTECT(elem3=allocVector(INTSXP,1)); pc2++;
                  INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt,l)));
                  PROTECT(elem2=setListElement(elem2,CHAR(STRING_ELT(altnam,l)), elem3)); pc2++;
                }else{
                  elem3=getListElement(elem2,CHAR(STRING_ELT(altnam,l)));
                  if(length(elem3)==0){
                    PROTECT(elem3=allocVector(INTSXP,1)); pc2++;
                    INTEGER(elem3)[0]=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt,l)));
                  }else{
                    PROTECT(elem3=coerceVector(elem3,INTSXP)); pc2++;
                    INTEGER(elem3)[0]+=MIN(asInteger(VECTOR_ELT(alt,k)), asInteger(VECTOR_ELT(alt,l)));
                  }
                  PROTECT(elem2=setListElement(elem2,CHAR(STRING_ELT(altnam,l)), elem3)); pc2++;
                }
                PROTECT(elem=setListElement(elem,CHAR(STRING_ELT(altnam,k)), elem2)); pc2++;
              }
              SET_VECTOR_ELT(sip,i,elem);  /*Write list into place*/
              UNPROTECT(pc2);
            }
          }
        }
      }
    }
    if(pc>1000){             /*If the stack is getting out of hand, fix it*/
      UNPROTECT(pc-4); pc=4;
    }
  }
  
  /*Write the lists into tri*/
  PROTECT(tri=allocVector(VECSXP,0)); pc++;
  PROTECT(tri=setListElement(tri,"top",top)); pc++;
  PROTECT(tri=setListElement(tri,"tip",tip)); pc++;
  PROTECT(tri=setListElement(tri,"sop",sop)); pc++;
  PROTECT(tri=setListElement(tri,"sip",sip)); pc++;
  
  /*Unprotect and return*/
  UNPROTECT(pc);
  return tri;
}


void logrm_irr(SEXP lrm, int n, SEXP rrl, double coef, int mode)
{
  int pc=0,i,j;
  SEXP list,src;

  /*Rprintf("Entered logrm_irr; cf=%f, m=%d\n",coef,mode);*/
  /*Get the source IDs*/
  PROTECT(src = coerceVector(getAttrib(rrl,R_NamesSymbol),INTSXP)); pc++;
  
  /*Compute the multipliers*/
  if(mode==0){                        /*Recency affects i->j*/
    for(i=0;i<length(rrl);i++){
      PROTECT(list=coerceVector(VECTOR_ELT(rrl,i),INTSXP)); pc++;
      for(j=0;j<length(list);j++)
        REAL(lrm)[INTEGER(src)[i]-1+(INTEGER(list)[j]-1)*n]+=coef/(j+1.0);
    }
  }else{                           /*Recency affects j->i*/
    for(i=0;i<length(rrl);i++){
      PROTECT(list=coerceVector(VECTOR_ELT(rrl,i),INTSXP)); pc++;
      for(j=0;j<length(list);j++)
        REAL(lrm)[INTEGER(list)[j]-1+(INTEGER(src)[i]-1)*n]+=coef/(j+1.0);
    }
  }

  /*Unprotect and return*/
  UNPROTECT(pc);
}

void logrm_irr_samp(SEXP lrv, int ns, int *tail, int *head, int n, SEXP rrl, double coef, int mode)
{
  int pc=0,i,j,k,flag;
  SEXP list,src;

  /*Rprintf("Entered logrm_irr; cf=%f, m=%d\n",coef,mode);*/
  /*Get the source IDs*/
  PROTECT(src = coerceVector(getAttrib(rrl,R_NamesSymbol),INTSXP)); pc++;
  
  /*Compute the multipliers*/
  if(mode==0){                        /*Recency affects i->j*/
    for(i=0;i<ns;i++){
      flag=0;
      for(j=0;(!flag)&&(j<length(src));j++){
        if(tail[i]==INTEGER(src)[j]){
          PROTECT(list=coerceVector(VECTOR_ELT(rrl,j),INTSXP)); pc++;
          for(k=0;k<length(list);k++){
            if(head[i]==INTEGER(list)[k]){
              REAL(lrv)[i]+=coef/(k+1.0);
              flag++;
            }
          }
        }
      }
    }
  }else{                           /*Recency affects j->i*/
    for(i=0;i<ns;i++){
      flag=0;
      for(j=0;(!flag)&&(j<length(src));j++){
        if(head[i]==INTEGER(src)[j]){
          PROTECT(list=coerceVector(VECTOR_ELT(rrl,j),INTSXP)); pc++;
          for(k=0;k<length(list);k++){
            if(tail[i]==INTEGER(list)[k]){
              REAL(lrv)[i]+=coef/(k+1.0);
              flag++;
            }
          }
        }
      }
    }
  }

  /*Unprotect and return*/
  UNPROTECT(pc);
}


void logrm_rceff(SEXP lrm, int m, int nvar, int n, int it, int v, double *inparm, double *outparm, double coef, int mode)
{
  int i,j;

  /*Rprintf("Entered logrm_rceff; n=%d, cf=%f, m=%d\n",n,coef,mode);*/
  /*Add the row/column/combined effects to the current log rate matrix*/
  switch(mode){
    case 0:                      /*Row effects mode*/
      for(i=0;i<n;i++)
        for(j=0;j<n;j++){
//          Rprintf("\t%d %d\n",i,j);
//          Rprintf("\t\toutparm[%d]=%f\n",i,outparm[it+v*m+i*nvar*m]);
//          Rprintf("\t\tREAL(lrm)[%d+%d*%d]=%f\n",i,j,n,REAL(lrm)[i+j*n]);
          REAL(lrm)[i+j*n]+=outparm[it+v*m+i*nvar*m]*coef;
        }
      break;
//    Rprintf("\tDone with rceff loop\n");
    case 1:                      /*Column effects mode*/
      for(i=0;i<n;i++)
        for(j=0;j<n;j++)
          REAL(lrm)[i+j*n]+=inparm[it+v*m+j*nvar*m]*coef;
      break;
    case 2:                      /*Row/column effects mode (prod)*/
      for(i=0;i<n;i++)
        for(j=0;j<n;j++)
          REAL(lrm)[i+j*n]+=outparm[it+v*m+i*nvar*m]*inparm[it+v*m+j*nvar*m]* coef;
      break;
    case 3:                      /*Row/column effects mode (sum)*/
      for(i=0;i<n;i++)
        for(j=0;j<n;j++)
          REAL(lrm)[i+j*n]+=(outparm[it+v*m+i*nvar*m]+inparm[it+v*m+j*nvar*m])* coef;
      break;
    case 4:                      /*Event-wise effects mode*/
      for(i=0;i<n;i++)
        for(j=0;j<n;j++)
          REAL(lrm)[i+j*n]+=outparm[it+v*m+i*nvar*m+j*n*nvar*m]*coef;
      break;
  }
//  Rprintf("Leaving rceff\n");
}

void logrm_rceff_samp(SEXP lrv, int ns, int *tail, int *head, int m, int nvar, int n, int it, int v, double *inparm, double *outparm, double coef, int mode)
{
  int i;

//  Rprintf("Entered logrm_rceff_samp; ns=%d, cf=%f, m=%d\n",ns,coef,mode);
  /*Add the row/column/combined effects to the current log rate matrix*/
  switch(mode){
    case 0:                      /*Row effects mode*/
      for(i=0;i<ns;i++){
        REAL(lrv)[i]+=outparm[it+v*m+tail[i]*nvar*m]*coef;
      }
      break;
//    Rprintf("\tDone with rceff loop\n");
    case 1:                      /*Column effects mode*/
      for(i=0;i<ns;i++){
//        Rprintf("t=%d, h=%d, out=%f, in=%f\n",tail[i],head[i],outparm[it+v*m+tail[i]*nvar*m], inparm[it+v*m+head[i]*nvar*m]);
        REAL(lrv)[i]+=inparm[it+v*m+head[i]*nvar*m]*coef;
      }
//      Rprintf("\tDone with rceff loop\n");
      break;
    case 2:                      /*Row/column effects mode (prod)*/
      for(i=0;i<ns;i++){
        REAL(lrv)[i]+=outparm[it+v*m+tail[i]*nvar*m]* inparm[it+v*m+head[i]*nvar*m]*coef;
      }
      break;
    case 3:                      /*Row/column effects mode (sum)*/
      for(i=0;i<ns;i++){
//        Rprintf("t=%d, h=%d, out=%f, in=%f\n",tail[i],head[i],outparm[it+v*m+tail[i]*nvar*m], inparm[it+v*m+head[i]*nvar*m]);
        REAL(lrv)[i]+=(outparm[it+v*m+tail[i]*nvar*m]+ inparm[it+v*m+head[i]*nvar*m])*coef;
      }
//      Rprintf("\tDone with rceff loop\n");
      break;
    case 4:                      /*Eventwise effects mode*/
      for(i=0;i<ns;i++){
//        Rprintf("t=%d, h=%d, cov=%f\n",tail[i],head[i],outparm[it+v*m+tail[i]*nvar*m+head[i]*n*nvar*m]);
        REAL(lrv)[i]+=outparm[it+v*m+tail[i]*nvar*m+head[i]*n*nvar*m]*coef;
      }
//      Rprintf("\tDone with rceff loop\n");
      break;
  }
//  Rprintf("Leaving rceff\n");
}


void logrm_normint(SEXP lrm, int n, SEXP acl, double *deg, double coef, int mode)
/*Note: it is presumed that the acl used here only contains the iteration of interest; deg must be the vector (presumably degree) with respect to which normalization is to be performed.  When mode=0, the i->j interaction is normalized by deg[i], while deg[j] is used when mode=1.  Modes 2 and 3 are as per modes 1 and 2, only the rate for the i->j interaction is modified by the j->i interaction instead (and i,j are reversed generally).*/
{
  int i,j;

  /*Rprintf("Entered logrm_normint; cf=%f, m=%d\n",coef,mode);*/
  /*Compute the rate multipliers*/
  switch(mode){
    case 0:                     /*i->j normalized by deg[i]*/
      for(i=0;i<n;i++)
        for(j=0;j<n;j++)
          if(i!=j){
            if(deg[i]==0.0)
              REAL(lrm)[i+j*n]+=coef/((double)n-1.0);
            else
              REAL(lrm)[i+j*n]+=coef*acl_adj(acl,i,j)/deg[i];
          }
      break;
    case 1:                     /*i->j normalized by deg[j]*/
      for(i=0;i<n;i++)
        for(j=0;j<n;j++)
          if(i!=j){
            if(deg[j]==0.0)
              REAL(lrm)[i+j*n]+=coef/((double)n-1.0);
            else
              REAL(lrm)[i+j*n]+=coef*acl_adj(acl,i,j)/deg[j];
          }
      break;
    case 2:                     /*j->i normalized by deg[j]*/
      for(i=0;i<n;i++)
        for(j=0;j<n;j++)
          if(i!=j){
            if(deg[j]==0.0)
              REAL(lrm)[i+j*n]+=coef/((double)n-1.0);
            else
              REAL(lrm)[i+j*n]+=coef*acl_adj(acl,j,i)/deg[j];
          }
      break;
    case 3:                     /*j->i normalized by deg[i]*/
      for(i=0;i<n;i++)
        for(j=0;j<n;j++)
          if(i!=j){
            if(deg[i]==0.0)
              REAL(lrm)[i+j*n]+=coef/((double)n-1.0);
            else
              REAL(lrm)[i+j*n]+=coef*acl_adj(acl,j,i)/deg[i];
          }
      break;
  }
}

void logrm_normint_samp(SEXP lrv, int ns, int *tail, int *head, int n, SEXP acl, double *deg, double coef, int mode)
/*Note: it is presumed that the acl used here only contains the iteration of interest; deg must be the vector (presumably degree) with respect to which normalization is to be performed.  When mode=0, the i->j interaction is normalized by deg[i], while deg[j] is used when mode=1.  Modes 2 and 3 are as per modes 1 and 2, only the rate for the i->j interaction is modified by the j->i interaction instead (and i,j are reversed generally).*/
{
  int i;

  /*Rprintf("Entered logrm_normint; cf=%f, m=%d\n",coef,mode);*/
  /*Compute the rate multipliers*/
  switch(mode){
    case 0:                     /*i->j normalized by deg[i]*/
      for(i=0;i<ns;i++){
        if(deg[tail[i]]==0.0)
          REAL(lrv)[i]+=coef/((double)n-1.0);
        else
          REAL(lrv)[i]+=coef*acl_adj(acl,tail[i],head[i])/deg[tail[i]];
      }
      break;
    case 1:                     /*i->j normalized by deg[j]*/
      for(i=0;i<ns;i++){
        if(deg[head[i]]==0.0)
          REAL(lrv)[i]+=coef/((double)n-1.0);
        else
          REAL(lrv)[i]+=coef*acl_adj(acl,tail[i],head[i])/deg[head[i]];
      }
      break;
    case 2:                     /*j->i normalized by deg[j]*/
      for(i=0;i<ns;i++){
        if(deg[head[i]]==0.0)
          REAL(lrv)[i]+=coef/((double)n-1.0);
        else
          REAL(lrv)[i]+=coef*acl_adj(acl,head[i],tail[i])/deg[head[i]];
      }
      break;
    case 3:                     /*j->i normalized by deg[i]*/
      for(i=0;i<ns;i++){
        if(deg[tail[i]]==0.0)
          REAL(lrv)[i]+=coef/((double)n-1.0);
        else
          REAL(lrv)[i]+=coef*acl_adj(acl,head[i],tail[i])/deg[tail[i]];
      }
      break;
  }
}


void logrm_ladj(SEXP lrm, int n, SEXP adj, double coef, int mode)
{
  int pc=0,i,j,src;
  SEXP sids,dids,edges;

//  Rprintf("Entered logrm_ladj; cf=%f, m=%d\n",coef,mode);
  /*Get vertex IDs*/
  PROTECT(sids=coerceVector(getAttrib(adj, R_NamesSymbol),INTSXP)); pc++;
  
  /*Compute the multipliers*/
  if(mode==0){  /*Adjust src->dest by strength of src->dest connection*/
    for(i=0;i<length(sids);i++){
//      Rprintf("%d ->\n",INTEGER(sids)[i]);
      src=INTEGER(sids)[i];
      PROTECT(edges=coerceVector(VECTOR_ELT(adj,i),REALSXP)); pc++;
      PROTECT(dids=coerceVector(getAttrib(VECTOR_ELT(adj,i),R_NamesSymbol), INTSXP)); pc++;
      for(j=0;j<length(dids);j++){
//        Rprintf("\t%d (%f)\n",INTEGER(dids)[i],REAL(edges)[j]);
        REAL(lrm)[src-1+(INTEGER(dids)[j]-1)*n]+=coef*REAL(edges)[j];
      }
    }
  }else{  /*Adjust dest->src by strength of src->dest connection*/
    for(i=0;i<length(sids);i++){
      src=INTEGER(sids)[i];
      PROTECT(edges=coerceVector(VECTOR_ELT(adj,i),REALSXP)); pc++;
      PROTECT(dids=coerceVector(getAttrib(VECTOR_ELT(adj,i),R_NamesSymbol), INTSXP)); pc++;
      for(j=0;j<length(dids);j++)
        REAL(lrm)[INTEGER(dids)[j]-1+(src-1)*n]+=coef*REAL(edges)[j];
    }
  }

  UNPROTECT(pc);
}

void logrm_ladj_samp(SEXP lrv, int ns, int *tail, int *head, SEXP adj, double coef, int mode)
/*This version of the routine is based on the assumption that only sampled edges are to be evaluated.*/
{
  int pc=0,i,j,k,flag;
  SEXP sids,dids,edges;

  /*Rprintf("Entered logrm_ladj; cf=%f, m=%d\n",coef,mode);*/
  /*Get vertex IDs*/
  PROTECT(sids=coerceVector(getAttrib(adj, R_NamesSymbol),INTSXP)); pc++;
  
  /*Compute the multipliers*/
  if(mode==0){  /*Adjust src->dest by strength of src->dest connection*/
    for(i=0;i<ns;i++){
      flag=0;
      for(j=0;(j<length(sids))&&(!flag);j++)
        if(tail[i]==INTEGER(sids)[j]-1){
          PROTECT(dids=coerceVector(getAttrib(VECTOR_ELT(adj,j),R_NamesSymbol), INTSXP)); pc++;
          for(k=0;(k<length(dids))&&(!flag);k++)
            if(head[i]==INTEGER(dids)[k]-1){
              PROTECT(edges=coerceVector(VECTOR_ELT(adj,j),REALSXP)); pc++;
              REAL(lrv)[i]+=coef*REAL(edges)[k];
              flag++;
            }
        }
    }
  }else{  /*Adjust dest->src by strength of src->dest connection*/
    for(i=0;i<ns;i++){
      flag=0;
      for(j=0;(j<length(sids))&&(!flag);j++)
        if(head[i]==INTEGER(sids)[j]-1){
          PROTECT(dids=coerceVector(getAttrib(VECTOR_ELT(adj,j),R_NamesSymbol), INTSXP)); pc++;
          for(k=0;(k<length(dids))&&(!flag);k++)
            if(tail[i]==INTEGER(dids)[k]-1){
              PROTECT(edges=coerceVector(VECTOR_ELT(adj,j),REALSXP)); pc++;
              REAL(lrv)[i]+=coef*REAL(edges)[k];
              flag++;
            }
        }
    }
  }

  UNPROTECT(pc);
}


void lambda(SEXP pv, int it, SEXP effects, int nv, int m, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrm)
/*Calculate dyadic event log-rates under the relational event model, for a
specified iteration.  These are stored in the log-rate matrix, lrm.*/
{
  int pc=0,i,j,k,*e,pvc=0,nvar=0,ct;
  SEXP dv=R_NilValue,fe=R_NilValue,cptr=R_NilValue;

//  Rprintf("Entering lambda, it=%d, nv=%d, len(lrm)=%d, typeof(lrm)=%d\n",it,nv, length(lrm),TYPEOF(lrm));
  /*Initialize the log-rate matrix*/
  for(i=0;i<nv;i++)
    for(j=0;j<nv;j++)
      REAL(lrm)[i+j*nv]=0.0;
    
  /*Compute the effects*/
  e=INTEGER(effects);
  for(i=0;i<length(effects);i++){
//    Rprintf("\t%d=%d\n",i,e[i]);
    if(e[i])
      switch(i){
        case NIDEGSEND:    /*Norm indegree -> future sending rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=REAL(cumideg)[it+j*m]/(double)it;
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
//          Rprintf("\tlength(dv)=%d, dv[0]=%f, dv[1]=%f, dv[2]=%f\n",length(dv), REAL(dv)[0],REAL(dv)[1],REAL(dv)[2]);
          logrm_rceff(lrm,1,1,nv,0,0,REAL(dv),REAL(dv),REAL(pv)[pvc],0);
//          Rprintf("Returned from logrm_rceff\n");
          pvc++;     /*Increment the parameter vector count*/
          break;
        case NIDEGREC:     /*Norm indegree -> future receiving rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=REAL(cumideg)[it+j*m]/(double)it;
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
          logrm_rceff(lrm,1,1,nv,0,0,REAL(dv),REAL(dv),REAL(pv)[pvc],1);
          pvc++;     /*Increment the parameter vector count*/
          break;
        case NODEGSEND:    /*Norm outdegree -> future sending rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=REAL(cumodeg)[it+j*m]/(double)it;
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
          logrm_rceff(lrm,1,1,nv,0,0,REAL(dv),REAL(dv),REAL(pv)[pvc],0);
          pvc++;     /*Increment the parameter vector count*/
          break;
        case NODEGREC:     /*Norm outdegree -> future receiving rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=REAL(cumodeg)[it+j*m]/(double)it;
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
          logrm_rceff(lrm,1,1,nv,0,0,REAL(dv),REAL(dv),REAL(pv)[pvc],1);
          pvc++;     /*Increment the parameter vector count*/
          break;
        case NTDEGSEND:    /*Norm total degree -> future sending rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=(REAL(cumideg)[it+j*m]+REAL(cumodeg)[it+j*m])/ (2.0*(double)it);
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
          logrm_rceff(lrm,1,1,nv,0,0,REAL(dv),REAL(dv),REAL(pv)[pvc],0);
          pvc++;     /*Increment the parameter vector count*/
          break;
        case NTDEGREC:     /*Norm total degree -> future receiving rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=(REAL(cumideg)[it+j*m]+REAL(cumodeg)[it+j*m])/ (2.0*(double)it);
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
          logrm_rceff(lrm,1,1,nv,0,0,REAL(dv),REAL(dv),REAL(pv)[pvc],1);
          pvc++;     /*Increment the parameter vector count*/
          break;
        case FPSENDSEND:   /*Fraction past sending -> future sending rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          for(j=0;j<nv;j++)
            REAL(dv)[j]=REAL(cumodeg)[it+j*m];
          logrm_normint(lrm,nv,acl,REAL(dv),REAL(pv)[pvc++],0);
          break;
        case FPRECSEND:    /*Fraction past receipt -> future sending rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          for(j=0;j<nv;j++)
            REAL(dv)[j]=REAL(cumideg)[it+j*m];
          logrm_normint(lrm,nv,acl,REAL(dv),REAL(pv)[pvc++],3);
          break;
        case RRRECSEND:    /*Recency of receipt -> future sending rate*/
          logrm_irr(lrm,nv,VECTOR_ELT(getListElement(rrl,"in"),it), REAL(pv)[pvc++],0);
          break;
        case RRSENDSEND:   /*Recency of sending -> future sending rate*/
          logrm_irr(lrm,nv,VECTOR_ELT(getListElement(rrl,"out"),it), REAL(pv)[pvc++],0);
          break;
        case COVSEND:      /*Covariate effect for sending*/
          /*Determine the time slice to use*/
          ct=LOGICAL(getListElement(covar,"ctimed"))[0];
          j=(ct ? it : 0);
          /*Get the number of variables*/
          nvar=INTEGER(getListElement(covar,"ncov"))[0];
          /*Proceed with the updates*/
          cptr=getListElement(covar,"CovSnd");
          for(k=0;k<nvar;k++)
            logrm_rceff(lrm,ct*m+(1-ct),nvar,nv,j,k,REAL(cptr),REAL(cptr), REAL(pv)[pvc++],0);
          break;
        case COVREC:       /*Covariate effect for receiving*/
          /*Determine the time slice to use*/
          ct=LOGICAL(getListElement(covar,"ctimed"))[1];
          j=(ct ? it : 0);
          /*Get the number of variables*/
          nvar=INTEGER(getListElement(covar,"ncov"))[1];
          /*Proceed with the updates*/
          cptr=getListElement(covar,"CovRec");
          for(k=0;k<nvar;k++)
            logrm_rceff(lrm,ct*m+(1-ct),nvar,nv,j,k,REAL(cptr),REAL(cptr), REAL(pv)[pvc++],1);
          break;
        case COVSENDREC:   /*Covariate effect for sending and receiving*/
          /*Determine the time slice to use*/
          ct=LOGICAL(getListElement(covar,"ctimed"))[2];
          j=(ct ? it : 0);
          /*Get the number of variables*/
          nvar=INTEGER(getListElement(covar,"ncov"))[2];
          /*Proceed with the updates*/
          cptr=getListElement(covar,"CovInt");
          for(k=0;k<nvar;k++)
            logrm_rceff(lrm,ct*m+(1-ct),nvar,nv,j,k,REAL(cptr),REAL(cptr), REAL(pv)[pvc++],3);
          break;
        case COVEVENT:    /*Generic event-wise covariate effect*/
          /*Determine the time slice to use*/
          ct=LOGICAL(getListElement(covar,"ctimed"))[3];
          j=(ct ? it : 0);
          /*Get the number of variables*/
          nvar=INTEGER(getListElement(covar,"ncov"))[3];
          /*Proceed with the updates*/
          cptr=getListElement(covar,"CovEvent");
          for(k=0;k<nvar;k++)
            logrm_rceff(lrm,ct*m+(1-ct),nvar,nv,j,k,REAL(cptr),REAL(cptr), REAL(pv)[pvc++],4);
          break;
        case OTPSEND:      /*Outbound two-paths -> future sending rate*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(tri,"top"),it), REAL(pv)[pvc++],0);
          break;
        case ITPSEND:      /*Incoming two-paths -> future sending rate*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(tri,"tip"),it), REAL(pv)[pvc++],0);
          break;
        case OSPSEND:      /*Outbound shared partners -> future sending rate*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(tri,"sop"),it), REAL(pv)[pvc++],0);
          break;
        case ISPSEND:      /*Inbound shared partners -> future sending rate*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(tri,"sip"),it), REAL(pv)[pvc++],0);
          break;
        case FESEND:       /*Fixed effects for sending*/
          if(length(fe)==0){
            PROTECT(fe=allocVector(REALSXP,nv)); pc++;
          }
          for(j=1;j<nv;j++)
            REAL(fe)[j]=REAL(pv)[pvc+j-1];
          REAL(fe)[0]=0.0;
          logrm_rceff(lrm,1,1,nv,0,0,REAL(fe),REAL(fe),1.0,0);
          pvc+=nv-1;
          break;
        case FEREC:        /*Fixed effects for receiving*/
          if(length(fe)==0){
            PROTECT(fe=allocVector(REALSXP,nv)); pc++;
          }
          for(j=1;j<nv;j++)
            REAL(fe)[j]=REAL(pv)[pvc+j-1];
          REAL(fe)[0]=0.0;
          logrm_rceff(lrm,1,1,nv,0,0,REAL(fe),REAL(fe),1.0,1);
          pvc+=nv-1;
          break;
        case FESENDREC:    /*Fixed effects for sending and receiving*/
          if(length(fe)==0){
            PROTECT(fe=allocVector(REALSXP,nv)); pc++;
          }
          for(j=1;j<nv;j++)
            REAL(fe)[j]=REAL(pv)[pvc+j-1];
          REAL(fe)[0]=0.0;
          logrm_rceff(lrm,1,1,nv,0,0,REAL(fe),REAL(fe),1.0,3);
          pvc+=nv-1;
          break; 
        case PSABBA:       /*P-Shift (turn receiving): AB->BA (dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"ABBA"),it), REAL(pv)[pvc++],0);
          break;
        case PSABB0:       /*P-Shift (turn receiving): AB->B0 (non-dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"ABB0"),it), REAL(pv)[pvc++],0);
          break;
        case PSABBY:       /*P-Shift (turn receiving): AB->BY (dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"ABBY"),it), REAL(pv)[pvc++],0);
          break;
        case PSA0X0:       /*P-Shift (turn claiming): A0->X0 (non-dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"A0X0"),it), REAL(pv)[pvc++],0);
          break;
        case PSA0XA:       /*P-Shift (turn claiming): A0->XA (non-dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"A0XA"),it), REAL(pv)[pvc++],0);
          break;
        case PSA0XY:       /*P-Shift (turn claiming): A0->XY (non-dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"A0XY"),it), REAL(pv)[pvc++],0);
          break;
        case PSABX0:       /*P-Shift (turn usurping): AB->X0 (non-dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"ABX0"),it), REAL(pv)[pvc++],0);
          break;
        case PSABXA:       /*P-Shift (turn usurping): AB->XA (dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"ABXA"),it), REAL(pv)[pvc++],0);
          break;
        case PSABXB:       /*P-Shift (turn usurping): AB->XB (dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"ABXB"),it), REAL(pv)[pvc++],0);
          break;
        case PSABXY:       /*P-Shift (turn usurping): AB->XY (dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"ABXY"),it), REAL(pv)[pvc++],0);
          break;
        case PSA0AY:       /*P-Shift (turn continuing): A0->AY (non-dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"A0AY"),it), REAL(pv)[pvc++],0);
          break;
        case PSABA0:       /*P-Shift (turn continuing): AB->A0 (non-dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"ABA0"),it), REAL(pv)[pvc++],0);
          break;
        case PSABAY:       /*P-Shift (turn continuing): AB->AY (dyadic)*/
          logrm_ladj(lrm,nv,VECTOR_ELT(getListElement(ps,"ABAY"),it), REAL(pv)[pvc++],0);
          break;
      }
  }

  /*Unprotect and return*/
  UNPROTECT(pc);
}

void lambda_samp(SEXP pv, int it, SEXP effects, int nv, int m, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrv, int ns, int *tail, int *head)
{
  int pc=0,i,j,k,*e,pvc=0,nvar=0;
  SEXP dv=R_NilValue,fe=R_NilValue,cptr=R_NilValue;

//  Rprintf("Entering lambda, it=%d, nv=%d, len(lrv)=%d, typeof(lrv)=%d\n",it,nv, length(lrv),TYPEOF(lrv));
  /*Initialize the log-rate matrix*/
  for(i=0;i<ns;i++)
    REAL(lrv)[i]=0.0;
    
  /*Compute the effects*/
  e=INTEGER(effects);
  for(i=0;i<length(effects);i++){
//    Rprintf("\t%d=%d\n",i,e[i]);
    if(e[i])
      switch(i){
        case NIDEGSEND:    /*Norm indegree -> future sending rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=REAL(cumideg)[it+j*m]/(double)it;
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
//          Rprintf("\tlength(dv)=%d, dv[0]=%f, dv[1]=%f, dv[2]=%f\n",length(dv), REAL(dv)[0],REAL(dv)[1],REAL(dv)[2]);
          logrm_rceff_samp(lrv,ns,tail,head,1,1,nv,0,0,REAL(dv),REAL(dv), REAL(pv)[pvc],0);
//          Rprintf("Returned from logrm_rceff\n");
          pvc++;     /*Increment the parameter vector count*/
          break;
        case NIDEGREC:     /*Norm indegree -> future receiving rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=REAL(cumideg)[it+j*m]/(double)it;
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
          logrm_rceff_samp(lrv,ns,tail,head,1,1,nv,0,0,REAL(dv),REAL(dv), REAL(pv)[pvc],1);
          pvc++;     /*Increment the parameter vector count*/
          break;
        case NODEGSEND:    /*Norm outdegree -> future sending rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=REAL(cumodeg)[it+j*m]/(double)it;
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
          logrm_rceff_samp(lrv,ns,tail,head,1,1,nv,0,0,REAL(dv),REAL(dv), REAL(pv)[pvc],0);
          pvc++;     /*Increment the parameter vector count*/
          break;
        case NODEGREC:     /*Norm outdegree -> future receiving rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=REAL(cumodeg)[it+j*m]/(double)it;
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
          logrm_rceff_samp(lrv,ns,tail,head,1,1,nv,0,0,REAL(dv),REAL(dv), REAL(pv)[pvc],1);
          pvc++;     /*Increment the parameter vector count*/
          break;
        case NTDEGSEND:    /*Norm total degree -> future sending rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=(REAL(cumideg)[it+j*m]+REAL(cumodeg)[it+j*m])/ (double)it;
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
          logrm_rceff_samp(lrv,ns,tail,head,1,1,nv,0,0,REAL(dv),REAL(dv), REAL(pv)[pvc],0);
          pvc++;     /*Increment the parameter vector count*/
          break;
        case NTDEGREC:     /*Norm total degree -> future receiving rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          if(it>0){                             /*Normalize the degree*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=(REAL(cumideg)[it+j*m]+REAL(cumodeg)[it+j*m])/ (double)it;
          }else{                                /*Initialize if needed*/
            for(j=0;j<nv;j++)
              REAL(dv)[j]=1.0/(nv-1.0);
          }
//          Rprintf("\tlength(dv)=%d, dv[0]=%f, dv[1]=%f, dv[2]=%f\n",length(dv), REAL(dv)[0],REAL(dv)[1],REAL(dv)[2]);
          logrm_rceff_samp(lrv,ns,tail,head,1,1,nv,0,0,REAL(dv),REAL(dv), REAL(pv)[pvc],1);
          pvc++;     /*Increment the parameter vector count*/
          break;
        case FPSENDSEND:   /*Fraction past sending -> future sending rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          for(j=0;j<nv;j++)
            REAL(dv)[j]=REAL(cumodeg)[it+j*m];
          logrm_normint_samp(lrv,ns,tail,head,nv,acl,REAL(dv),REAL(pv)[pvc++], 0);
          break;
        case FPRECSEND:    /*Fraction past receipt -> future sending rate*/
          if(length(dv)==0){
            PROTECT(dv=allocVector(REALSXP,nv)); pc++;
          }
          for(j=0;j<nv;j++)
            REAL(dv)[j]=REAL(cumideg)[it+j*m];
          logrm_normint_samp(lrv,ns,tail,head,nv,acl,REAL(dv),REAL(pv)[pvc++], 3);
          break;
        case RRRECSEND:    /*Recency of receipt -> future sending rate*/
          logrm_irr_samp(lrv,ns,tail,head,nv, VECTOR_ELT(getListElement(rrl,"in"), it), REAL(pv)[pvc++],0);
          break;
        case RRSENDSEND:   /*Recency of sending -> future sending rate*/
          logrm_irr_samp(lrv,ns,tail,head,nv, VECTOR_ELT(getListElement(rrl,"out"), it), REAL(pv)[pvc++],0);
          break;
        case COVSEND:      /*Covariate effect for sending*/
          /*Determine the time slice to use*/
          j=(LOGICAL(getListElement(covar,"ctimed"))[0] ? it : 0);
          /*Get the number of variables*/
          nvar=INTEGER(getListElement(covar,"ncov"))[0];
          /*Proceed with the updates*/
          cptr=getListElement(covar,"CovSnd");
          for(k=0;k<nvar;k++)
            logrm_rceff_samp(lrv,ns,tail,head,j*m+(1-j),nvar,nv,j,k,REAL(cptr), REAL(cptr),REAL(pv)[pvc++],0);
          break;
        case COVREC:       /*Covariate effect for receiving*/
          /*Determine the time slice to use*/
          j=(LOGICAL(getListElement(covar,"ctimed"))[1] ? it : 0);
          /*Get the number of variables*/
          nvar=INTEGER(getListElement(covar,"ncov"))[1];
          /*Proceed with the updates*/
          cptr=getListElement(covar,"CovRec");
          for(k=0;k<nvar;k++)
            logrm_rceff_samp(lrv,ns,tail,head,j*m+(1-j),nvar,nv,j,k,REAL(cptr), REAL(cptr),REAL(pv)[pvc++],1);
          break;
        case COVSENDREC:   /*Covariate effect for sending and receiving*/
          /*Determine the time slice to use*/
          j=(LOGICAL(getListElement(covar,"ctimed"))[2] ? it : 0);
          /*Get the number of variables*/
          nvar=INTEGER(getListElement(covar,"ncov"))[2];
          /*Proceed with the updates*/
          cptr=getListElement(covar,"CovInt");
          for(k=0;k<nvar;k++)
            logrm_rceff_samp(lrv,ns,tail,head,j*m+(1-j),nvar,nv,j,k,REAL(cptr), REAL(cptr),REAL(pv)[pvc++],3);
          break;
        case COVEVENT:    /*Generic element-wise covariate effect*/
          /*Determine the time slice to use*/
          j=(LOGICAL(getListElement(covar,"ctimed"))[3] ? it : 0);
          /*Get the number of variables*/
          nvar=INTEGER(getListElement(covar,"ncov"))[3];
          /*Proceed with the updates*/
          cptr=getListElement(covar,"CovEvent");
          for(k=0;k<nvar;k++)
            logrm_rceff_samp(lrv,ns,tail,head,j*m+(1-j),nvar,nv,j,k,REAL(cptr), REAL(cptr),REAL(pv)[pvc++],4);
          break;
        case OTPSEND:      /*Outbound two-paths -> future sending rate*/
          logrm_ladj_samp(lrv,ns,tail,head, VECTOR_ELT(getListElement(tri,"top"),it), REAL(pv)[pvc++],0);
          break;
        case ITPSEND:      /*Incoming two-paths -> future sending rate*/
          logrm_ladj_samp(lrv,ns,tail,head, VECTOR_ELT(getListElement(tri,"tip"),it), REAL(pv)[pvc++],0);
          break;
        case OSPSEND:      /*Outbound shared partners -> future sending rate*/
          logrm_ladj_samp(lrv,ns,tail,head, VECTOR_ELT(getListElement(tri,"sop"),it), REAL(pv)[pvc++],0);
          break;
        case ISPSEND:      /*Inbound shared partners -> future sending rate*/
          logrm_ladj_samp(lrv,ns,tail,head, VECTOR_ELT(getListElement(tri,"sip"),it), REAL(pv)[pvc++],0);
          break;
        case FESEND:       /*Fixed effects for sending*/
          if(length(fe)==0){
            PROTECT(fe=allocVector(REALSXP,nv)); pc++;
          }
          for(j=1;j<nv;j++)
            REAL(fe)[j]=REAL(pv)[pvc+j-1];
          REAL(fe)[0]=0.0;
          logrm_rceff_samp(lrv,ns,tail,head,1,1,nv,0,0,REAL(fe),REAL(fe),1.0,0);
          pvc+=nv-1;
          break;
        case FEREC:        /*Fixed effects for receiving*/
          if(length(fe)==0){
            PROTECT(fe=allocVector(REALSXP,nv)); pc++;
          }
          for(j=1;j<nv;j++)
            REAL(fe)[j]=REAL(pv)[pvc+j-1];
          REAL(fe)[0]=0.0;
          logrm_rceff_samp(lrv,ns,tail,head,1,1,nv,0,0,REAL(fe),REAL(fe),1.0,1);
          pvc+=nv-1;
          break;
        case FESENDREC:    /*Fixed effects for sending and receiving*/
          if(length(fe)==0){
            PROTECT(fe=allocVector(REALSXP,nv)); pc++;
          }
          for(j=1;j<nv;j++)
            REAL(fe)[j]=REAL(pv)[pvc+j-1];
          REAL(fe)[0]=0.0;
          logrm_rceff_samp(lrv,ns,tail,head,1,1,nv,0,0,REAL(fe),REAL(fe),1.0,3);
          pvc+=nv-1;
          break; 
        case PSABBA:       /*P-Shift (turn receiving): AB->BA (dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head, VECTOR_ELT(getListElement(ps,"ABBA"),it), REAL(pv)[pvc++],0);
          break;
        case PSABB0:       /*P-Shift (turn receiving): AB->B0 (non-dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"ABB0"),it), REAL(pv)[pvc++],0);
          break;
        case PSABBY:       /*P-Shift (turn receiving): AB->BY (dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"ABBY"),it), REAL(pv)[pvc++],0);
          break;
        case PSA0X0:       /*P-Shift (turn claiming): A0->X0 (non-dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"A0X0"),it), REAL(pv)[pvc++],0);
          break;
        case PSA0XA:       /*P-Shift (turn claiming): A0->XA (non-dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"A0XA"),it), REAL(pv)[pvc++],0);
          break;
        case PSA0XY:       /*P-Shift (turn claiming): A0->XY (non-dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"A0XY"),it), REAL(pv)[pvc++],0);
          break;
        case PSABX0:       /*P-Shift (turn usurping): AB->X0 (non-dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"ABX0"),it), REAL(pv)[pvc++],0);
          break;
        case PSABXA:       /*P-Shift (turn usurping): AB->XA (dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"ABXA"),it), REAL(pv)[pvc++],0);
          break;
        case PSABXB:       /*P-Shift (turn usurping): AB->XB (dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"ABXB"),it), REAL(pv)[pvc++],0);
          break;
        case PSABXY:       /*P-Shift (turn usurping): AB->XY (dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"ABXY"),it), REAL(pv)[pvc++],0);
          break;
        case PSA0AY:       /*P-Shift (turn continuing): A0->AY (non-dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"A0AY"),it), REAL(pv)[pvc++],0);
          break;
        case PSABA0:       /*P-Shift (turn continuing): AB->A0 (non-dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"ABA0"),it), REAL(pv)[pvc++],0);
          break;
        case PSABAY:       /*P-Shift (turn continuing): AB->AY (dyadic)*/
          logrm_ladj_samp(lrv,ns,tail,head,          VECTOR_ELT(getListElement(ps,"ABAY"),it), REAL(pv)[pvc++],0);
          break;
      }
  }

  /*Unprotect and return*/
  UNPROTECT(pc);
}


SEXP drem_n2llik_R(SEXP pv, SEXP effects, SEXP edgelist, SEXP n, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrm, SEXP ordinal, SEXP condnum)
/*Calculate deviance for the relational event model.*/
{
  int pc=0,i,m,nv,j,k,ncond;
  double lrsum,*el,tdelta;
  SEXP ll,aclit;

  /*Set things up*/
  PROTECT(ll=allocVector(REALSXP,1)); pc++;
  REAL(ll)[0]=0.0;
  PROTECT(lrm=coerceVector(lrm,REALSXP)); pc++;
  PROTECT(pv=coerceVector(pv,REALSXP)); pc++;
  PROTECT(effects=coerceVector(effects,LGLSXP)); pc++;
  m=nrows(edgelist);
  PROTECT(edgelist=coerceVector(edgelist,REALSXP)); pc++;
  PROTECT(n=coerceVector(n,INTSXP)); pc++;
  nv=INTEGER(n)[0];
  PROTECT(cumideg=coerceVector(cumideg,REALSXP)); pc++;
  PROTECT(cumodeg=coerceVector(cumodeg,REALSXP)); pc++;
  PROTECT(ordinal=coerceVector(ordinal,LGLSXP)); pc++;
  el=REAL(edgelist);
  PROTECT(condnum=coerceVector(condnum,INTSXP)); pc++;
  ncond=INTEGER(condnum)[0];

  /*Find the log-likelihood*/
  for(i=ncond;i<m;i++){
    if(length(acl)>0)
      aclit=VECTOR_ELT(acl,i);
    else
      aclit=R_NilValue;
//      Rprintf("Calling lambda(%d)\n",i);
    lambda(pv,i,effects,nv,m,aclit,cumideg,cumodeg,rrl,covar, ps,tri,lrm);
//    Rprintf("Finished lambda(%d); calculating normalizing factor\n",i);
    lrsum=-DBL_MAX;
    for(j=0;j<nv;j++)
      for(k=0;k<nv;k++)
        if(j!=k)
          lrsum=logsum(REAL(lrm)[j+k*nv],lrsum);
//    if(INTEGER(ordinal)[0]||(i<m-1))
//      Rprintf("Calculating log-likelihood; (%d,%d) log-rate is %f, norm is %f\n", (int)el[i+m]-1,(int)el[i+2*m]-1, REAL(lrm)[(int)el[i+m]-1+((int)el[i+2*m]-1)*nv],lrsum);
//    else
//      Rprintf("Calculating log-likelihood contribution for right-censoring; norm is %f\n",lrsum);
    if(INTEGER(ordinal)[0]){                   /*log( lambda/sum(lambda') )*/
      REAL(ll)[0]+=REAL(lrm)[(int)el[i+m]-1+((int)el[i+2*m]-1)*nv]-lrsum;
    }else if(i<(m-1)){            /*log( lambda exp(-tdelta*sum(lambda')) )*/
      if(i>0)
        tdelta=el[i]-el[i-1];
      else
        tdelta=el[0];
      REAL(ll)[0]+=REAL(lrm)[(int)el[i+m]-1+((int)el[i+2*m]-1)*nv] - tdelta*exp(lrsum);
    }else{                               /*log( exp(-tdelta*sum(lambda')) )*/
      REAL(ll)[0]-=(el[i]-el[i-1])*exp(lrsum);
    }
  }

  REAL(ll)[0]*=-2.0;
  UNPROTECT(pc);
  return ll;
}


SEXP drem_n2llik_samp_R(SEXP pv, SEXP effects, SEXP edgelist, SEXP n, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrv, SEXP tail, SEXP head, SEXP ordinal, SEXP condnum)
/*Calculate deviance for the relational event model, using dyad sampling.  This
degrades the estimates slightly, but can greatly conserve computational
resources.*/
{
  int pc=0,i,m,nv,j,ncond;
  double lrsum,*el,tdelta;
  SEXP ll,aclit;

  /*Set things up*/
  PROTECT(ll=allocVector(REALSXP,1)); pc++;
  REAL(ll)[0]=0.0;
  PROTECT(lrv=coerceVector(lrv,REALSXP)); pc++;
  PROTECT(tail=coerceVector(tail,INTSXP)); pc++;
  PROTECT(head=coerceVector(head,INTSXP)); pc++;
  PROTECT(pv=coerceVector(pv,REALSXP)); pc++;
  PROTECT(effects=coerceVector(effects,LGLSXP)); pc++;
  m=nrows(edgelist);
  PROTECT(edgelist=coerceVector(edgelist,REALSXP)); pc++;
  PROTECT(n=coerceVector(n,INTSXP)); pc++;
  nv=INTEGER(n)[0];
  PROTECT(cumideg=coerceVector(cumideg,REALSXP)); pc++;
  PROTECT(cumodeg=coerceVector(cumodeg,REALSXP)); pc++;
  PROTECT(ordinal=coerceVector(ordinal,LGLSXP)); pc++;
  el=REAL(edgelist);
  PROTECT(condnum=coerceVector(condnum,INTSXP)); pc++;
  ncond=INTEGER(condnum)[0];

  /*Find the log-likelihood*/
  for(i=ncond;i<m;i++){
    if((i<m-1)||INTEGER(ordinal)[0]){
      INTEGER(tail)[0]=(int)REAL(edgelist)[i+m]-1;
      INTEGER(head)[0]=(int)REAL(edgelist)[i+2*m]-1;
    }else{
      INTEGER(tail)[0]=(int)REAL(edgelist)[i-1+m]-1;
      INTEGER(head)[0]=(int)REAL(edgelist)[i-1+2*m]-1;
    }
    if(length(acl)>0)
      aclit=VECTOR_ELT(acl,i);
    else
      aclit=R_NilValue;
    lambda_samp(pv,i,effects,nv,m,aclit,cumideg,cumodeg,rrl,covar, ps,tri,lrv,length(lrv),INTEGER(tail),INTEGER(head));
    lrsum=-DBL_MAX;
    for(j=0;j<length(lrv);j++)
      lrsum=logsum(REAL(lrv)[j],lrsum);
    lrsum+=-log((double)length(lrv))+2.0*log((double)nv);
    if(INTEGER(ordinal)[0]){
      REAL(ll)[0]+=REAL(lrv)[0]-lrsum;
    }else if(i<m-1){
      if(i>0)
        tdelta=el[i]-el[i-1];
      else
        tdelta=el[0];
      REAL(ll)[0]+=REAL(lrv)[0]-tdelta*exp(lrsum);
    }else{
      REAL(ll)[0]-=(el[i]-el[i-1])*exp(lrsum);
    }
  }

  REAL(ll)[0]*=-2.0;
  UNPROTECT(pc);
  return ll;
}


SEXP drem_gof_R(SEXP pv, SEXP effects, SEXP edgelist, SEXP n, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrm, SEXP ordinal, SEXP condnum)
/*Compute goodness-of-fit stats for a fitted relational event model, including
the event-wise deviance residuals and predicted events.*/
{
  int pc=0,i,m,nv,j,k,ncond;
  double lrsum,maxlrm,*el,obslr,lrjk,ldt,*dc;
  SEXP ll,pred,outlist,aclit,erank,devcen;
  
  /*Set things up*/
  m=nrows(edgelist);
  PROTECT(n=coerceVector(n,INTSXP)); pc++;
  nv=INTEGER(n)[0];
  PROTECT(condnum=coerceVector(condnum,INTSXP)); pc++;
  ncond=INTEGER(condnum)[0];
  PROTECT(ordinal=coerceVector(ordinal,LGLSXP)); pc++;
  PROTECT(ll=allocVector(REALSXP,m-1+INTEGER(ordinal)[0]-ncond)); pc++;
  PROTECT(devcen=allocVector(REALSXP,1)); pc++;
  dc=REAL(devcen);
  PROTECT(pred=allocVector(INTSXP,2*(m-1+INTEGER(ordinal)[0]-ncond))); pc++;    
  PROTECT(erank=allocVector(INTSXP,m-1+INTEGER(ordinal)[0]-ncond)); pc++;    
  PROTECT(lrm=coerceVector(lrm,REALSXP)); pc++;
  PROTECT(pv=coerceVector(pv,REALSXP)); pc++;
  PROTECT(effects=coerceVector(effects,LGLSXP)); pc++;
  PROTECT(edgelist=coerceVector(edgelist,REALSXP)); pc++;
  PROTECT(cumideg=coerceVector(cumideg,REALSXP)); pc++;
  PROTECT(cumodeg=coerceVector(cumodeg,REALSXP)); pc++;
  PROTECT(outlist=allocVector(VECSXP,0)); pc++;
  el=REAL(edgelist);

  /*Calculate deviance residuals and predictions*/
  for(i=ncond;i<m-1+INTEGER(ordinal)[0];i++){
    if(length(acl)>0)
      aclit=VECTOR_ELT(acl,i);
    else
      aclit=R_NilValue;
    lambda(pv,i,effects,nv,m,aclit,cumideg,cumodeg,rrl,covar,ps,tri,lrm);
    if(INTEGER(ordinal)[0])
      lrsum=-DBL_MAX;
    else
      lrsum=0.0;
    maxlrm=-DBL_MAX;
    obslr=REAL(lrm)[(int)el[i+m]-1+((int)el[i+2*m]-1)*nv]; /*Obs hazard*/
    INTEGER(erank)[i-ncond]=1;              /*Initialize the edge rank*/
    if(i>0)                                /*Record time since prior event*/
      ldt=log(el[i]-el[i-1]);
    else
      ldt=log(el[i]);
    for(j=0;j<nv;j++)
      for(k=0;k<nv;k++)
        if(j!=k){
          lrjk=REAL(lrm)[j+k*nv];           /*Pull the j->k log rate*/
          if(INTEGER(ordinal)[0])
            lrsum=logsum(lrjk,lrsum);       /*Increment the log rate sum*/
          else
            lrsum+=exp(lrjk+ldt);           /*Increment the survival sum*/
          if(lrjk>maxlrm){                  /*Record the most probable edge*/
            INTEGER(pred)[i-ncond]=j+1;
            INTEGER(pred)[i+m-1+INTEGER(ordinal)[0]-ncond]=k+1;
            maxlrm=lrjk;                    /*Update the maximum*/
          }
          if(lrjk>obslr)     /*If j->k more likely, increment obs event rank*/
            INTEGER(erank)[i-ncond]++;
        }
    REAL(ll)[i-ncond]=-2.0*(obslr-lrsum);
  }
  /*Deal with censoring factor, if not ordinal; note, i carries from above*/
  if(!(INTEGER(ordinal)[0])){
    if(length(acl)>0)
      aclit=VECTOR_ELT(acl,i);
    else
      aclit=R_NilValue;
    lambda(pv,i,effects,nv,m,aclit,cumideg,cumodeg,rrl,covar,ps,tri,lrm);
    lrsum=0.0;
    if(i>0)                                /*Record time since prior event*/
      ldt=log(el[i]-el[i-1]);
    else
      ldt=log(el[i]);
    for(j=0;j<nv;j++)
      for(k=0;k<nv;k++)
        if(j!=k){
          lrjk=REAL(lrm)[j+k*nv];           /*Pull the j->k log rate*/
          lrsum+=exp(lrjk+ldt);             /*Increment the survival sum*/
        }
    *dc=2.0*lrsum;
  }else{
    *dc=0.0;
  }

  /*Add deviance residuals and prediction information to the output list*/
  PROTECT(outlist=setListElement(outlist,"residuals",ll)); pc++;
  PROTECT(outlist=setListElement(outlist,"predicted",pred)); pc++;
  PROTECT(outlist=setListElement(outlist,"obs.rank",erank)); pc++;
  PROTECT(outlist=setListElement(outlist,"dev.censor",devcen)); pc++;

  /*Unprotect and return*/
  UNPROTECT(pc);
  return outlist;
}


void rem_int_dev_R(double *par, int *pnpar, double *evm, int *pm, double *statsa, int *pnet, int *suppm, int *calcderiv, double *val, double *grad, double *hess)
/*Deviance calculation for a single general-form event sequence (interval timing
likelihood).  Arguments are as follows:
  par - vector of parameters (pre-mapped, if that was required)
  pnpar - length of par
  evm - event matrix (m x 2, 1st col is type code (0=exog), 2nd col is time)
  pm - number of events in the sequence (m)
  statsa - stats array (m x net x npar)
  pnet - number of event types (not including the 0th type)
  suppm - support matrix (m x net, 1 if possible, 0 if not possible)
  calcderiv - 1 if derivatives should be calculated, 0 otherwise
  val - pre-allocated slot for the deviance value (length 1)
  grad - pre-allocated slot for the gradient (length npar)
  hess - pre-allocated slot for the hessian (npar x npar)

  Note that although this function is intended to be used to obtain the
  deviance, it itself actually calculates the log likelihood; the conversion
  to deviance is done on the R side.  
*/
{
  int i,j,k,l,m,net,npar;
  double dt,lp,dtelp;
 
  /*Initialize stuff*/
  m=*pm;
  net=*pnet;
  npar=*pnpar;
  *val=0.0;
  if(*calcderiv){
    for(i=0;i<npar;i++){
      grad[i]=0.0;
      for(j=0;j<npar;j++)
        hess[i+j*npar]=0.0;
    }
  }
  
  /*Perform the deviance calculations*/
  for(i=0;i<m;i++){
    if(i>0)                                     /*Get the time since last event*/
      dt=evm[i+m]-evm[i-1+m];
    else
      dt=evm[m];
    for(j=0;j<net;j++)                        /*Walk through each event type...*/
      if(suppm[i+j*m]){                       /*...ignoring impossible events. */
        /*Calculate the linear predictor for this (potential) event*/
        for(k=0,lp=0.0;k<npar;k++)
          lp+=par[k]*statsa[i+j*m+k*m*net];
        dtelp=dt*exp(lp);               /*Save a few operations by precomputing*/
        /*If this is the observed event, increment value and gradient*/
        if(((int)evm[i])==j+1){
          *val+=lp;                                /*Add log hazard to deviance*/
          if(*calcderiv)                            /*Add raw stats to gradient*/
            for(k=0;k<npar;k++)
              grad[k]+=statsa[i+j*m+k*m*net];
        }
        /*Observed or not, add increments from survival function*/
        *val-=dtelp;
        if(*calcderiv){
          for(k=0;k<npar;k++){
            grad[k]-=statsa[i+j*m+k*m*net]*dtelp;
            for(l=k;l<npar;l++){
              hess[k+l*npar]-=statsa[i+j*m+k*m*net]*statsa[i+j*m+l*m*net]*dtelp;
              hess[l+k*npar]=hess[k+l*npar];
            }
          }
        }
      }
  }
}


void rem_int_ev_dev_R(double *par, int *pnpar, double *ev, double *statsm, int *pnet, int *suppv, int *calcderiv, double *val, double *grad, double *hess, int *initvals)
/*Deviance calculation for a single event in a general-form event sequence (interval timing
likelihood).  Arguments are as follows:
  par - vector of parameters (pre-mapped, if that was required)
  pnpar - length of par
  ev - event vector (length 2, 1st element is type code (0=exog), 2nd is time since last event)
  statsm - stats matrix (net x npar)
  pnet - number of event types (not including the 0th type)
  suppv - support vector (length net, 1 if possible, 0 if not possible)
  calcderiv - 1 if derivatives should be calculated, 0 otherwise
  val - pre-allocated slot for the deviance value (length 1)
  grad - pre-allocated slot for the gradient (length npar)
  hess - pre-allocated slot for the hessian (npar x npar)
  initvals - 1 if we should initialize val/grad/etc., 0 otherwise

  Note that although this function is intended to be used to obtain the
  deviance, it itself actually calculates the log likelihood; the conversion
  to deviance is done on the R side.  
*/
{
  int i,j,k,l,net,npar;
  double dt,lp,dtelp;
 
  /*Initialize stuff*/
  net=*pnet;
  npar=*pnpar;
  if(*initvals){  /*If needed, initialize value/gradient/hessian*/
    *val=0.0;
    if(*calcderiv){
      for(i=0;i<npar;i++){
        grad[i]=0.0;
        for(j=0;j<npar;j++)
          hess[i+j*npar]=0.0;
      }
    }
  }
  
  /*Perform the deviance calculations*/
  dt=ev[1];                                   /*Get the time since last event*/
  for(j=0;j<net;j++)                          /*Walk through each event type...*/
    if(suppv[j]){                             /*...ignoring impossible events. */
      /*Calculate the linear predictor for this (potential) event*/
      for(k=0,lp=0.0;k<npar;k++)
        lp+=par[k]*statsm[j+k*net];
      dtelp=dt*exp(lp);               /*Save a few operations by precomputing*/
      /*If this is the observed event, increment value and gradient*/
      if(((int)ev[0])==j+1){
        *val+=lp;                                /*Add log hazard to deviance*/
        if(*calcderiv)                           /*Add raw stats to gradient*/
          for(k=0;k<npar;k++)
            grad[k]+=statsm[j+k*net];
      }
      /*Observed or not, add increments from survival function*/
      *val-=dtelp;
      if(*calcderiv){
        for(k=0;k<npar;k++){
          grad[k]-=statsm[j+k*net]*dtelp;
          for(l=k;l<npar;l++){
            hess[k+l*npar]-=statsm[j+k*net]*statsm[j+l*net]*dtelp;
            hess[l+k*npar]=hess[k+l*npar];
          }
        }
      }
    }
}


void rem_ord_dev_R(double *par, int *pnpar, int *evm, int *pm, double *statsa, int *pnet, int *suppm, int *calcderiv, double *val, double *grad, double *hess)
/*Deviance calculation for a single general-form event sequence (ordinal timing
likelihood).  Arguments are as follows:
  par - vector of parameters (pre-mapped, if that was required)
  pnpar - length of par
  evm - event vector (length m)
  pm - number of events in the sequence (m)
  statsa - stats array (m x net x npar)
  pnet - number of event types (not including the 0th type)
  suppm - support matrix (m x net, 1 if possible, 0 if not possible)
  calcderiv - 1 if derivatives should be calculated, 0 otherwise
  val - pre-allocated slot for the deviance value (length 1)
  grad - pre-allocated slot for the gradient (length npar)
  hess - pre-allocated slot for the hessian (npar x npar)
  
  Note that although this function is intended to be used to obtain the
  deviance, it itself actually calculates the log likelihood; the conversion
  to deviance is done on the R side.  
*/
{
  int i,j,k,l,m,net,npar;
  double lp,elp,selp,*tselp=NULL,*ttselp=NULL,tselpdselp;
 
  /*Initialize stuff*/
  m=*pm;
  net=*pnet;
  npar=*pnpar;
  *val=0.0;
  if(*calcderiv){
    tselp=(double *)R_alloc(npar,sizeof(double));
    ttselp=(double *)R_alloc(npar*npar,sizeof(double));
    for(i=0;i<npar;i++){
      grad[i]=0.0;
      for(j=0;j<npar;j++)
        hess[i+j*npar]=0.0;
    }
  }
  
  /*Perform the deviance calculations*/
  for(i=0;i<m;i++)
    if(evm[i]>0){                                /*Exogenous events are ignored*/
      selp=0.0;                                    /*Initialize summation terms*/
      if(*calcderiv){
        for(j=0;j<npar;j++){
          tselp[j]=0.0;
          for(k=0;k<npar;k++)
            ttselp[j+k*npar]=0.0;
        }
      }
      for(j=0;j<net;j++)                      /*Walk through each event type...*/
        if(suppm[i+j*m]){                     /*...ignoring impossible events. */
          /*Calculate the linear predictor for this (potential) event*/
          for(k=0,lp=0.0;k<npar;k++)
            lp+=par[k]*statsa[i+j*m+k*m*net];
          elp=exp(lp);                        /*Store to save a few exp() calls*/
          /*If this is the observed event, increment value and gradient*/
          if(evm[i]==j+1){
            *val+=lp;                              /*Add log hazard to deviance*/
            if(*calcderiv)                          /*Add raw stats to gradient*/
              for(k=0;k<npar;k++)
                grad[k]+=statsa[i+j*m+k*m*net];
          }
          /*If needed, accumulate more elements for the normalizing factor*/
          selp+=elp;                                       /*Accum sum(exp(lp))*/
          if(*calcderiv){
            for(k=0;k<npar;k++){                      /*Accum sum(t(a)*exp(lp))*/
              tselp[k]+=statsa[i+j*m+k*m*net]*elp;
              for(l=k;l<npar;l++){               /*Accum sum(t(a)*t(b)*exp(lp))*/
                ttselp[k+l*npar]+=statsa[i+j*m+k*m*net]*statsa[i+j*m+l*m*net]* elp;
              }
            }
          }
        }
      /*Add normalizing factor elements*/
      *val-=log(selp);
      if(*calcderiv){
        for(j=0;j<npar;j++){
          tselpdselp=tselp[j]/selp;
          grad[j]-=tselpdselp;
          for(k=j;k<npar;k++)
            hess[j+k*npar]-=(ttselp[j+k*npar]-tselpdselp*tselp[k])/selp;
        }
      }
    }
  if(*calcderiv){
    for(j=0;j<npar;j++)                     /*Fill in lower triangle of hessian*/
      for(k=j+1;k<npar;k++)
        hess[k+j*npar]=hess[j+k*npar];
  }
}


SEXP lambda_R(SEXP pv, SEXP iter, SEXP effects, SEXP n, SEXP nev, SEXP acl, SEXP cumideg, SEXP cumodeg, SEXP rrl, SEXP covar, SEXP ps, SEXP tri, SEXP lrm)
/*Calculate dyadic event log-rates under the relational event model, for a
specified iteration.  These are stored in the log-rate matrix, lrm, which
is also returned.  The arguments are as follows:

pv - parameter vector
iter - the iteration to use (employed when consulting past history)
effects - logical vector of included effects
n - number of vertices
nev - number of events in the current history
acl - accumulated local communication list
cumideg - accumulated indegree
cumodeg - accumulated outdegree
rrl - recency lists
covar - covariates
ps - pshift indicators
tri - accumulated triadic effects
lrm - n x n hazard matrix

Note that this is just an interface to lambda.  Using it dynamically is
not efficient, but it's a thing that can be done.
*/
{
  int pc=0,it,nv,m;

  /*Validate the inputs*/
  PROTECT(nev=coerceVector(nev,INTSXP)); pc++;
  m=INTEGER(nev)[0];
  PROTECT(n=coerceVector(n,INTSXP)); pc++;
  nv=INTEGER(n)[0];
  PROTECT(iter=coerceVector(iter,INTSXP)); pc++;
  it=INTEGER(iter)[0];
  if(it<1)
    error("Can't compute on iteration number <1.\n");
  if(it>m)
    error("Can't compute on iteration number >nev.\n");
  PROTECT(lrm=coerceVector(lrm,REALSXP)); pc++;
  PROTECT(pv=coerceVector(pv,REALSXP)); pc++;
  PROTECT(effects=coerceVector(effects,LGLSXP)); pc++;
  PROTECT(cumideg=coerceVector(cumideg,REALSXP)); pc++;
  PROTECT(cumodeg=coerceVector(cumodeg,REALSXP)); pc++;

  /*Call lambda*/
  lambda(pv, it-1, effects, nv, m, acl, cumideg, cumodeg, rrl, covar, ps, tri, lrm);

  /*Unprotect and return*/
  UNPROTECT(pc);
  return lrm;
}
