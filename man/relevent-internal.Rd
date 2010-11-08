\name{network-internal}
\alias{accum.interact}
\alias{acl.adjmat}
\alias{acl.adj}
\alias{acl.deg}
\alias{accum.ps}
\alias{acl.ps}
\alias{acl.tri}
\alias{accum.rrl}
\title{Internal Relevent Package Functions}


\description{
  Internal relevent package functions.
}
\usage{
accum.interact(elist)
acl.adj(acl, iter, src, dest)
acl.adjmat(acl, n, iter)
acl.deg(acl, n, cmode=c("in", "out", "total"))
accum.ps(elist)
acl.ps(elist)
acl.tri(acl)
accum.rrl(elist)
}
\details{
  Most of these are not to be called by the user.
}
\arguments{
 \item{elist}{a three-column (time, source, destination) dyadic event list, sorted in ascending temporal order.}
 \item{acl}{a nested list structure, of the form iteration by ego by alter, containing accumulated dyadic event counts at each event onset.}
 \item{iter}{iteration (i.e., event) number.}
 \item{src}{integer denoting the sender of a dyadic event.}
 \item{dest}{integer denoting the receiver of a dyadic event.}
 \item{n}{the number of actors eligible to send/receive events.}
 \item{cmode}{the type of degree to be calculated.}
%\item{\dots}{further arguments passed to or used by methods.}
}
\seealso{\code{\link{rem}}}
\keyword{internal}
