######################################################################
#
# zzz.R
#
# Written by Carter T. Butts <buttsc@uci.edu>.
#
# Last Modified 02/29/12
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/relevent package
#
# .onLoad is run when the package is loaded with library(relevent)
#
######################################################################

.onLoad <- function(lib, pkg){
   temp<-packageDescription("relevent")
   msg<-(paste(temp$Package,": ",temp$Title,"\n",
               "Version ",temp$Version,
               " created on ",
                temp$Date,".\n", sep=""))
   msg<-paste(msg,"copyright (c) 2007, Carter T. Butts, University of California-Irvine\n",sep="")
   msg<-paste(msg,'Type help(package="relevent") to get started.\n',sep="")
   packageStartupMessage(msg)
}
