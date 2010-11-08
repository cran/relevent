\name{rem}
\alias{rem}
\alias{rem.ord.dev}
\alias{rem.int.dev}
\alias{rem.ord.nlp}
\alias{rem.int.nlp}
\alias{print.rem}
\alias{summary.rem}
\alias{print.summary.rem}
\title{
Fit a Relational Event Model to Single or Multiple Sequence Data
}
\description{
  Fits a relational event model to general event sequence data, using either the ordinal or interval time likelihoods.  Maximum likelihood and posterior mode methods are supported, as are local (per sequence) parameters and sequences with exogenous events.
}
\usage{
rem(eventlist, statslist, supplist = NULL, timing = c("ordinal", 
    "interval"), estimator = c("BPM", "MLE", "BMCMC", "BSIR"), 
    prior.param = list(mu = 0, sigma = 1000, nu = 4), mcmc.draws = 1500, 
    mcmc.thin = 25, mcmc.burn = 2000, mcmc.chains = 3, mcmc.sd = 0.05, 
    mcmc.ind.int = 50, mcmc.ind.sd = 10, sir.draws = 1000, sir.expand = 10, 
    sir.nu = 4, verbose = FALSE)
\method{print}{rem}(x, ...)
\method{summary}{rem}(object, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eventlist}{
a two-column matrix (or list thereof) containing the observed event sequence and timing information.
}
  \item{statslist}{
an event number by event type by statistic array (or list thereof) containing the sufficient statistics for the model to be estimated.
}
  \item{supplist}{
an event number by event type logical array (or list thereof) indicating which events were potentially observable at each point in the event history.
}
  \item{timing}{
the type of timing information to be used during estimation; \code{"ordinal"} indicates that only event order should be employed, while \code{"interval"} uses the exact inter-event times.
}
  \item{estimator}{
the type of estimator to be used; \code{"MLE"} selects maximum likelihood estimation, \code{"BPM"} selects Bayesian posterior mode estimation, \code{"BMCMC"} selects Bayesian posterior mean estimation via MCMC, and \code{"BSIR"} selects Bayesian posterior mean estimation via simulated importance resampling. 
}
  \item{prior.param}{
for the Bayesian methods, the prior parameters to be employed; currently, these are the location, scale, and degrees of freedom parameters for independent t priors, and may be given as vectors (to set different priors for each parameter).  (By default, a diffuse, heavy-tailed t distribution is used.)
}
  \item{mcmc.draws}{
total number of posterior draws to take when using the BMCMC method.
}
  \item{mcmc.thin}{
thinning interval for MCMC draws (BMCMC method).
}
  \item{mcmc.burn}{
number of burn-in iterations to use for each MCMC chain (BMCMC method).
}
  \item{mcmc.chains}{
number of MCMC chains to use (BMCMC method).
}
  \item{mcmc.sd}{
standard deviation for the random walk Metropolis sampler (BMCMC method).
}
  \item{mcmc.ind.int}{
interval at which to take draws from the independence sampler (versus the random walk Metropolis sampler).  (BMCMC method).
}
  \item{mcmc.ind.sd}{
standard deviation for the MCMC independence sampler (BMCMC method).
}
  \item{sir.draws}{
number of SIR draws to take (BSIR method).
}
  \item{sir.expand}{
expansion factor for the SIR sample; intitial sample size is \code{sir.draws} multiplied by \code{sir.expand}.
}
  \item{sir.nu}{
degrees of freedom parameter for the SIR sampling distribution.
}
  \item{verbose}{
logical; should verbose progress information be displayed?
}
  \item{x}{ an object of class \code{rem}.}
  \item{object}{ an object of class \code{rem}.}
  \item{...}{additional arguments.}
}
\details{
\code{rem} fits a general relational event model to one or more event sequences (or \dQuote{histories}), using either full interval or ordinal timing information.  Although particularly applicable to \dQuote{egocentric} relational event data, \code{rem} can be used to fit nearly any standard relational event model; the function depends heavily on user-supplied statistics, however, and thus lacks the built-in functionality of a routine like \code{\link{rem.dyad}}.  Four estimation methods are currently supported: maximum likelihood estimation, Bayesian posterior mode estimation, Bayesian posterior mean estimation via MCMC, and Bayesian posterior mean estimation via sampling importance resampling (SIR).  For the Bayesian methods, adjustable independent t priors are employed.  For both mode-based methods, estimates of uncertainty (standard errors or posterior standard deviations) are approximated using the appropriate inverse hessian matrix; for the two simulation-based methods, posterior standard deviations are estimated from the resulting sample.  

Irrespective of whether Bayesian or frequentist methods are used, the relevant likelihood is either based entirely on the order of events (\code{timing="ordinal"}) or on the realized event times (\code{timing="interval"}).  In the latter case, all event times are understood to be relative to the onset of observation (i.e., observation starts at time 0), and the last event time given is taken to be the end of the observation period.  (This should generally be marked as exogenous -- see below.)

Event source/target/content are handled generically by \code{rem} via \emph{event types}.  Each event must be of a given type, and any number of types may be employed (up to limits of time and memory).  Effects within the relational event model are associated with user-supplied statistics, of which any number may again be supplied (model identification notwithstanding).  At each point in the event history, it is possible that only particular types of events may be realized; this constraint can be specified by means of an optional user-supplied support structure.  Finally, it is also possible that an event sequence may be punctuated by \emph{exogenous events,} which are unmodeled but which may affect the endogenous event dynamics.  These are supported by means of a tacit \dQuote{exogenous} event type, which is handled by the estimation routine as appropriate for the specified likelihood.

Observed event data is supplied to \code{rem} via the \code{eventlist} argument.  For each event history, the observed events are indicated by a two-column matrix, whose \eqn{i}th row contains respectively the event type (as an integer ranging from 1 to the number of event types, inclusive) and the event time for the \eqn{i}th event in the history.  (The second column may be omitted in the ordinal case, and will in any event be ignored.)  Events must be given in ascending temporal order; if multiple histories are being modeled simultaneously (e.g., as with egocentric relational event samples), then \code{eventlist} should be a list with one matrix per event history.  Exogenous events, if present, are indicated by specifying an event type of 0.  (Note that the \dQuote{type} of an exogenous event is irrelevant, since any such properties of exogenous events are handled via the model statistics.)  If exact timing information is used, the hazard for the first event implicitly begins at time 0, and observation implicitly ends with the time of the last event (which should properly be coded as exogenous, unless the sampling design was based on observation of an endogenous terminal event).  Where applicable, censoring due to the sampling interval is accounted for in the data likelihood (assuming that the user has set the model statistics appropriately).

Statistics for the relational event model are specified in a manner somewhat analogous to that of \code{eventlist}.  Like the latter, \code{statslist} is generally a list with one element per event history, or a single element where only a single history is to be examined.  Each element of \code{statslist} should be a list containing either one or two three-dimensional arrays, with the first dimension indexing event order (from first event to last, including exogenous events where applicable), the second indexing event type (in order corresponding to the integer values of \code{eventlist}), and the third indexing the model statistics.  The \eqn{ijk}th cell of a \code{statslist} array is thus the value of the \eqn{k}th statistic prospectively impacting the hazard of observing an event of type \eqn{j} as the \eqn{i}th event in the history (given the previous \eqn{i-1} realized events).  Models estimated by \emph{rem} are regular in the sense that one parameter is estimated per statistic; intuitively, a large value of a \eqn{ijk}th \code{statslist} cell associatd with a large (positive) parameter represents an increased hazard of observing a type \eqn{j} event at the \eqn{i}th point in the respective history, while the same statistic associated with a highly negative parameter represents a correspondingly diminished hazard of observing said event.  (The total hazard of a given event type is equal to \eqn{\exp(\theta^T s_{ij})}{exp(theta^T s_ij)}, where \eqn{\theta}{theta} is the vector of model parameters and \eqn{s_{ij}}{s_ij} is the corresponding vector of sufficient statistics for a type \eqn{j} event given the \eqn{i-1} previously realized events; see the reference below for details.)  It is up to the user to supply these statistics, and moreover to ensure that they are well-behaved (e.g., not linearly dependent).  An array within a \code{statslist} element may be designated as \emph{global} or \emph{local} by assigning it to the appropriately named list element.  Statistics belonging to a global array are assumed to correspond to parameters that are homogeneous across event histories, and are estimated in a pooled fashion; if global arrays are supplied, they must be given for every element of \code{statslist} (and must carry the same statistics and event types, although these statistics will not typically take the same values).  Statistics belonging to a local array, on the other hand, are taken as idiosyncratic to the event history in question, and their corresponding parameters are estimated locally.  Both local and global statistics may be employed simultaneously if desired, but at least one must be specified in any case.  \code{rem} will return an error if passed a \code{statslist} with obvious inconsistencies.

If desired, support constraints for the event histories can be specified using \code{supplist}.  \code{supplist} should be a list with one element per history, each of which should be an event order by event type logical matrix.  The \eqn{ij}th cell of this matrix should be \code{TRUE} if an event of type \code{j} was a possible next event given the preceding code \eqn{i-1} events, and \code{FALSE} otherwise.  (By default, all events are assumed to be possible at all times.)  As with the model statistics, the elements of the support list must be user supplied, and will often be history-dependent.  (E.g., in a model for spell-based data, event types will come in onset/termination pairs, with terminal events necessarily being preceded by corresponding onset events.)

Given the above structure, \code{rem} will attempt to find a maximum likelihood or posterior estimate for the model parameters, as appropriate given \code{estimator}.  In the latter case, the prior parameters for each parameter may be set using \code{prior.param}.  Each parameter is taken to be \emph{a priori} t distributed, with the indicated location, scale, and degree of freedom parameters; by default, a fairly diffuse and heavy-tailed prior is used.  By specifying the elements of \code{prior.param} as vectors, it is possible to employ different priors for each model parameter.  In this case, the vector elements are used in the order of the statistics (first global, then each local in order by event history).  Standard errors or posterior standard deviation estimates are returned as appropriate, along with various goodness-of-fit indices.  (Bear in mind that the \dQuote{p-values} shown in the summary method for the posterior mode case are based on posterior quantiles (under an assumption of asymptotic normality), and should be interpreted in this fashion.)

For the MCMC sampling method, a combined independence and random walk Metropolis scheme is employed.  Proposals are multivariate Gaussian, with standard deviations as set via the appropriate arguments.  (These may be given as vectors, with one entry per parameter, if desired.)  Gewke and Gelman-Rubin MCMC diagnostics (produced by the \code{coda} package) are computed, and are stored as elements \code{geweke} and \code{gelman.rubin} within the model fit object.  The posterior draws themselves are stored as an element called \code{draws} within the model fit object, with corresponding log-posterior values \code{lp}.

The SIR method initially seeks the posterior mode (identically to the BPM method), and obtains approximate scale information using the Hessian of the log-posterior surface.  This is used to generate a set of approximate posterior draws via a multivariate t distribution centered on the posterior mode, with degrees of freedom given by \code{sir.nu}.  This crude sample is then refined by importance resampling, the final result of which is stored as element \code{draws} (with log-posterior vector \code{lp}) in the model fit object.  As with the BMCMC procedure, posterior mean and standard deviations are estimated from the final sample, although the mode information is retained in elements \code{coef.mode} and \code{cov.hess}.

As a general matter, the MLE and BPM methods are most dependent upon asymptotic assumptions, but are also (usually) the least computationally complex.  BMCMC requires no such assumptions, but can be extremely slow (and, like all MCMC methods, depends upon the quality of the MCMC sample).  The BSIR method is something of a compromise between BPM and BMCMC, starting with a mode approximation but refining it in the direction of the true posterior surface; as one might expect, its cost is also intermediate between these extremes.  For well-behaved models on large data sets, all methods are likely to produce nearly identical results.  The simulation-based methods (particularly BMCMC) may be safer in less salutary circumstances.  (Tests conducted by the author have so far obtained the best overall results from the BPM, particularly vis a vis estimates of uncertainty -- this advice may or may not generalize, however.)
}
\value{
An object of class \code{rem}, for which \code{\link{print}} and \code{\link{summary}} methods currently exist.
}
\references{
Butts, C.T.  (2008).  \dQuote{A Relational Event Framework for Social Action.}  \emph{Sociological Methodology}, 38(1). 
}
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
%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
