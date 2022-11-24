\name{as.sociomatrix.eventlist}
\alias{as.sociomatrix.eventlist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert an Event List Into a Sociomatrix
}
\description{
Convert a dyadic event list into an adjacency matrix, such that the \eqn{i,j} cell value is the number of \eqn{(i,j)} events in the list.
}
\usage{
as.sociomatrix.eventlist(eventlist, n = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eventlist}{
a three-column numeric matrix (or equivalent), containing the event list to be converted.
}
  \item{n}{
the number of vertices.  If omitted, this is assumed to be contained in an attribute called \code{"n"} attached to \code{eventlist}.
}
}
\details{
An event list must be a three-column matrix (or something that can be treated as one), whose second and third columns must contain vertex IDs; these can be given as characters, but must be coercable with \code{as.numeric} to numeric form.  Vertex IDs must be integers from \code{1:n}, where \code{n} is either supplied as an argument, or attached as an attribute of the \code{eventlist} object.  The first column of an eventlist matrix conventionally contains the event time, and is ignored; the second and third should contain the IDs of the senders and receivers of events (respectively).  Rows with missing values for one or both vertex IDs are removed during processing (but \code{NA}s in the first column have no effect, since the event timing information is not used).

The resulting output is an \code{n} by \code{n} adjacency matrix, whose \code{i,j} cell is the total number of events in \code{eventlist} from vertex \code{i} to vertex \code{j}.  This can be useful for visualizing or otherwise analyzing the time-marginalized structure of a dyadic interaction network.
}
\value{
A sociomatrix containing the time-aggregated event counts.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rem.dyad}}
}
\examples{
#Create a simple event list
el <- cbind(1:6, c(4,4,4,1,3,4), c(1,2,3,3,1,1))

#Convert to matrix form
as.sociomatrix.eventlist(el, 4)

#Can also store n as an attribute
attr(el, "n") <- 4
as.sociomatrix.eventlist(el)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{array}
\keyword{graphs}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
