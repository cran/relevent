######################################################################
#
# data.proc.R
#
# Written by Carter T. Butts <buttsc@uci.edu>.
#
# Last Modified 11/23/22
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/relevent package
#
# Contents:
#
# as.sociomatrix.eventlist
#
######################################################################


#Simple function to convert an event list to a sociomatrix, with the i,j
#cell value being the number of (i,j) events in the list
as.sociomatrix.eventlist <-function (eventlist, n=NULL) 
{
    #If the number of nodes is not given, attempt to obtain from an attr
    if(is.null(n))
      n<-attr(eventlist,"n")
    if(is.null(n))
      stop("n not provided, and could not be obtained from eventlist.\n")
    #Validate the event list
    if(NCOL(eventlist)!=3)
      stop("eventlist must be a three-column (event, sender, receiver) matrix.\n")
    #Drop any missing data
    eventlist <- eventlist[!apply(is.na(eventlist[,2:3]),1,any),,drop=FALSE]
    #Carry on
    g <- matrix(0, n, n)             #Create the adjacency matrix
    if (NROW(eventlist) > 0) {       #Tabulate event counts
        tabmat <- table(as.data.frame(eventlist[, -1, drop = FALSE]))
        sndnum<-as.numeric(dimnames(tabmat)[[1]])     #Senders used
        recnum<-as.numeric(dimnames(tabmat)[[2]])     #Receivers used
        if(any(is.na(sndnum))||any(is.na(recnum)))    #Should be numeric
          stop("eventlist contains non-numeric vertex IDs.\n")
        if(!all(sndnum%in%(1:n)))                     #Should be in range
          stop("Illegal sender ID - should be in [1,",n,"]")
        if(!all(recnum%in%(1:n)))                     #Should be in range
          stop("Illegal receiver ID - should be in [1,",n,"]")
        g[sndnum, recnum] <- tabmat
    }
    g
}

