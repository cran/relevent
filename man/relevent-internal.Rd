\name{relevent-internal}
\alias{accum.interact}
\alias{acl.adjmat}
\alias{acl.adj}
\alias{acl.deg}
\alias{accum.ps}
\alias{acl.ps}
\alias{acl.tri}
\alias{accum.rrl}
\alias{covarPrep}
\title{Internal Relevent Package Functions}


\description{
  Internal relevent package functions.
}
\usage{
accum.interact(elist, old.acl=NULL)
acl.adj(acl, iter, src, dest)
acl.adjmat(acl, n, iter)
acl.deg(acl, n, cmode=c("in", "out", "total"), old.deg=NULL)
accum.ps(elist)
acl.ps(elist, n, old.ps=NULL)
acl.tri(acl, old.tri=NULL)
accum.rrl(elist, old.rrl=NULL)
covarPrep(covar, n, m, effects=NULL)
}
\details{
  Most of these are not to be called by the user; they can be employed by the cogniscenti, but they may change without notice (so use at own risk).
}
\arguments{
 \item{elist}{a three-column (time, source, destination) dyadic event list, sorted in ascending temporal order.}
 \item{acl}{a nested list structure, of the form iteration by ego by alter, containing accumulated dyadic event counts at each event onset.}
 \item{iter}{iteration (i.e., event) number.}
 \item{src}{integer denoting the sender of a dyadic event.}
 \item{dest}{integer denoting the receiver of a dyadic event.}
 \item{n}{the number of actors eligible to send/receive events.}
 \item{cmode}{the type of degree to be calculated.}
 \item{old.acl}{previously computed acl structure to which new events should be added.}
 \item{old.deg}{previously computed cumulative degree structure to which events should be added.}
 \item{old.ps}{previously computed P-shift structure to which events should be added.}
 \item{old.tri}{previously computed triad structure to which events should be added.}
 \item{old.rrl}{previously computed recency structure to which events should be added.}
 \item{covar}{a covariate list, of the form passed normally to \code{\link{rem.dyad}}.}
 \item{m}{the intended number of events for which data should be checked.}
 \item{effects}{a logical effect inclusion vector of the type used internally by \code{\link{rem.dyad}}.}
%\item{\dots}{further arguments passed to or used by methods.}
}
\seealso{\code{\link{rem}}}
\keyword{internal}
