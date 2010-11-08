######################################################################
#
# zzz.R
#
# Written by Carter T. Butts <buttsc@uci.edu>.
#
# Last Modified 09/07/10
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/relevent package
#
# .onLoad is run when the package is loaded with library(relevent)
#
######################################################################

.onLoad <- function(lib, pkg){
   #library.dynam("relevent", pkg, lib)
    if(R.version$major=="1"){
     ehelp <- help(package="relevent")$info[[2]][[2]]
     cat(paste("'",ehelp[4],"'\n",
               "Version ",ehelp[2],
               " created on ",ehelp[3],".\n", sep=""))
    }else{
     ehelp <- help(package="relevent")$info[[1]]
     cat(paste(substring(ehelp[4],first=16),"\n",
               "Version ",substring(ehelp[2],first=16),
               " created on ",
                substring(ehelp[3],first=16),".\n", sep=""))
    }
    cat(paste("copyright (c) 2007, Carter T. Butts, University of California-Irvine\n",sep=""))
    cat('Type help(package="relevent") to get started.\n')
}
