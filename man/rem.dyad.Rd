\name{rem.dyad}
\alias{rem.dyad.lambda}
\alias{rem.dyad.lprior}
\alias{rem.dyad.nlpost}
\alias{rem.dyad.n2llik}
\alias{rem.dyad.n2llik.samp}
\alias{rem.dyad.nlpost.samp}
\alias{rem.dyad.gof}
\alias{rem.dyad}
\alias{rem.dyad}
\alias{print.rem.dyad}
\alias{print.summary.rem.dyad}
\alias{summary.rem.dyad}
\alias{simulate.rem.dyad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit a Relational Event Model to Dyadic Data }
\description{
  Fits a relational event model to dyadic edgelist data, using either the ordinal or temporal likelihood.  Maximum likelihood, posterior mode, and posterior importance resampling methods are supported.
}
\usage{
rem.dyad(edgelist, n, effects = NULL, ordinal = TRUE, acl = NULL,
    cumideg = NULL, cumodeg = NULL, rrl = NULL, covar = NULL, ps = NULL,
    tri = NULL, optim.method = "BFGS", optim.control = list(), 
    coef.seed = NULL, hessian = FALSE, sample.size = Inf, verbose = TRUE, 
    fit.method = c("BPM", "MLE", "BSIR"), conditioned.obs = 0, 
    prior.mean = 0, prior.scale = 100, prior.nu = 4, sir.draws = 500, 
    sir.expand = 10, sir.nu = 4, gof = TRUE)
\method{print}{rem.dyad}(x, ...)
\method{summary}{rem.dyad}(object, ...) 
\method{simulate}{rem.dyad}(object, nsim = object$m, seed = NULL, 
    coef = NULL, covar = NULL, edgelist = NULL, redraw.timing = FALSE,
    redraw.events = FALSE, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edgelist}{ a three-column edgelist matrix, with each row containing (in order) the time/order, sender, and receiver for the event in question, or \code{NULL} to create a model skeleton (useful for simulation). }
  \item{n}{ number of senders/receivers. }
  \item{effects}{ a character vector indicating which effects to use; see below for specification. }
  \item{ordinal}{ logical; should the ordinal likelihood be used?  (If \code{FALSE}, the temporal likelihood is used instead.) }
  \item{acl}{ optionally, a pre-computed acl structure. }
  \item{cumideg}{ optionally, a pre-computed cumulative indegree structure. }
  \item{cumodeg}{ optionally, a pre-computed cumulative outdegree stucture. }
  \item{rrl}{ optionally, a pre-computed recency-ranked communications list. }
  \item{covar}{ an optional list of sender/receiver/event covariates. }
  \item{ps}{ optionally, a pre-computed p-shift matrix. }
  \item{tri}{ optionally, a pre-computed triad statistic structure. }
  \item{optim.method}{ the method to be used by \code{\link{optim}}. }
  \item{optim.control}{ additional control parameters to \code{\link{optim}}. }
  \item{coef.seed}{an optional vector of coefficients to use as the starting point for the optimization process; if \code{edgelist==NULL}, this is the vector of embedded coefficients for the model skeleton.}
   \item{hessian}{ logical; compute the hessian of the log-likelihood/posterior surface? }
  \item{sample.size}{ sample size to use when estimating the sum of event rates. }
  \item{verbose}{ logical; deliver progress reports? }
  \item{fit.method}{ method to use when fitting the model. }
  \item{conditioned.obs}{ the number of initial observations on which to condition when fitting the model (defaults to 0).}
  \item{prior.mean}{ for Bayesian estimation, location vector for prior distribution (multivariate-t).  (Can be a single value.)}
  \item{prior.scale}{ for Bayesian estimation, scale vector for prior distribution.  (Can be a single value.) }
  \item{prior.nu}{ for Bayesian estimation, degrees of freedom for prior distribution. (Setting this to \code{Inf} results in a Gaussian prior.) }
  \item{sir.draws}{ for sampling importance resampling method, the number of posterior draws to take (post-resampling). }
  \item{sir.expand}{ for sampling importance resampling method, the expansion factor to use in the initial (pre-resampling) sample; sample size is \code{sir.expand*sir.draws}.}
  \item{sir.nu}{ for sampling importance resampling method, the degrees of freedom for the t distribution used to obtain initial (pre-resampling) sample. }
  \item{gof}{logical; calculate goodness-of-fit information?}
  \item{x}{ an object of class \code{rem.dyad}.}
  \item{object}{ an object of class \code{rem.dyad}.}
  \item{nsim}{ number of events to simulate (defaults to the observed sequence length in the fitted model). }
  \item{seed}{ random number seed to use for simulation. }
  \item{coef}{ optional vector of coefficients to override those in the fitted model object, for simulation purposes.}
  \item{redraw.timing}{ logical; should any prespecified events in \code{edgelist} have their timings redrawn during simulation? }
  \item{redraw.events}{ logical; should any prespecified events in \code{edgelist} have their senders and receivers redrawn during simulation? }
  \item{...}{additional arguments.}
}
\details{
\code{rem.dyad} fits a (dyadic) relational event model to an event sequence, using either the full temporal or ordinal data likelihoods.  Three estimation methods are currently supported: maximum likelihood estimation, Bayesian posterior mode estimation, and Bayesian sampling importance resampling.  For the Bayesian methods, an adjustable multivariate-t (or, if \code{prior.nu==Inf}, Gaussian) prior is employed.  In the case of Bayesian sampling importance resampling, the posterior mode (and the hessian of the posterior about it) is used as the basis for a multivariate-t sample, which is then resampled via SIR methods to obtain an approximate set of posterior draws.  While this approximation is not guaranteed to work well, it is generally more robust than pure mode approximations (or, in the case of the MLE, estimates of uncertainty derived from the inverse hessian matrix).

Whether Bayesian or frequentist methods are used, the relevant likelihood is either based entirely on the order of events (\code{ordinal=TRUE}) or on the realized event times (\code{ordinal=FALSE}).  In the latter case, all event times are understood to be relative to the onset of observation (i.e., observation starts at time 0), and the last event time given is taken to be the end of the observation period.  (If an event is also specified, this event is ignored.)

Effects to be fit by \code{rem.dyad} are determined by the eponymous \code{effects} argument, a character vector which lists the effects to be used.  These are as follows:
\itemize{
\item \code{NIDSnd}: Normalized indegree of \eqn{v} affects \eqn{v}'s future sending rate
\item \code{NIDRec}: Normalized indegree of \eqn{v} affects \eqn{v}'s future receiving rate
\item \code{NODSnd}: Normalized outdegree of \eqn{v} affects \eqn{v}'s future sending rate
\item \code{NODRec}: Normalized outdegree of \eqn{v} affects \eqn{v}'s future receiving rate
\item \code{NTDegSnd}: Normalized total degree of \eqn{v} affects \eqn{v}'s future sending rate
\item \code{NTDegRec}: Normalized total degree of \eqn{v} affects \eqn{v}'s future receiving rate
\item \code{FrPSndSnd}: Fraction of \eqn{v}'s past actions directed to \eqn{v'} affects \eqn{v}'s future rate of sending to \eqn{v'}
\item \code{FrRecSnd}: Fraction of \eqn{v}'s past receipt of actions from \eqn{v'} affects \eqn{v}'s future rate of sending to \eqn{v'}
\item \code{RRecSnd}: Recency of receipt of actions from \eqn{v'} affects \eqn{v}'s future rate of sending to \eqn{v'}
\item \code{RSndSnd}: Recency of sending to \eqn{v'} affects \eqn{v}'s future rate of sending to \eqn{v'}
\item \code{CovSnd}: Covariate effect for outgoing actions (requires a \code{covar} entry of the same name)
\item \code{CovRec}: Covariate effect for incoming actions (requires a \code{covar} entry of the same name)
\item \code{CovInt}: Covariate effect for both outgoing and incoming actions (requires a \code{covar} entry of the same name)
\item \code{CovEvent}: Covariate effect for each \eqn{(v,v')} action (requires a \code{covar} entry of the same name)
\item \code{OTPSnd}: Number of outbound two-paths from \eqn{v} to \eqn{v'} affects \eqn{v}'s future rate of sending to \eqn{v'}
\item \code{ITPSnd}: Number of incoming two-paths from \eqn{v'} to \eqn{v} affects \eqn{v}'s future rate of sending to \eqn{v'}
\item \code{OSPSnd}: Number of outbound shared partners for \eqn{v} and \eqn{v'} affects \eqn{v}'s future rate of sending to \eqn{v'}
\item \code{ISPSnd}: Number of inbound shared partners for \eqn{v} and \eqn{v'} affects \eqn{v}'s future rate of sending to \eqn{v'}
\item \code{FESnd}: Fixed effects for outgoing actions
\item \code{FERec}: Fixed effects for incoming actions
\item \code{FEInt}: Fixed effects for both outgoing and incoming actions
\item \code{PSAB-BA}: P-Shift effect (turn receiving) -- AB->BA (dyadic)
\item \code{PSAB-B0}: P-Shift effect (turn receiving) -- AB->B0 (non-dyadic)
\item \code{PAAB-BY}: P-Shift effect (turn receiving) -- AB->BY (dyadic)
\item \code{PSA0-X0}: P-Shift effect (turn claiming) -- A0->X0 (non-dyadic)
\item \code{PSA0-XA}: P-Shift effect (turn claiming) -- A0->XA (non-dyadic)
\item \code{PSA0-XY}: P-Shift effect (turn claiming) -- A0->XY (non-dyadic)
\item \code{PSAB-X0}: P-Shift effect (turn usurping) -- AB->X0 (non-dyadic)
\item \code{PSAB-XA}: P-Shift effect (turn usurping) -- AB->XA (dyadic)
\item \code{PSAB-XB}: P-Shift effect (turn usurping) -- AB->XB (dyadic)
\item \code{PSAB-XY}: P-Shift effect (turn usurping) -- AB->XY (dyadic)
\item \code{PSA0-AY}: P-Shift effect (turn continuing) -- A0->AY (non-dyadic)
\item \code{PSAB-A0}: P-Shift effect (turn continuing) -- AB->A0 (non-dyadic)
\item \code{PSAB-AY}: P-Shift effect (turn continuing) -- AB->AY (dyadic)
}
Note that not all effects may lead to identified models in all cases - it is up to the user to ensure that the postulated model makes sense.

Data to be used by \code{rem.dyad} must consist of an edgelist matrix, whose rows contain information on successive events.  This matrix must have three columns, containing (respectively) the event times, sender IDs (as integers from 1 to \code{n}), and receiver IDs (also from 1 to \code{n}).  As already noted, event times should be relative to onset of observation where the temporal likelihood is being used; otherwise, only event order is employed.  In the temporal likelihood case, the last row should contain the time for the termination of the observation period -- any event on this row is ignored.  If \code{conditioned.obs>0}, the relevant number of initial observations is taken as fixed, and the likelihood of the remaining sequence is calculated conditional on these values; this can be useful when analyzing an event history with no clear starting point.

If covariates effects are indicated, then appropriate covariate values must be supplied as a list in argument \code{covar}.  The elements of \code{covar} should be given the same name as the effect type to which they correspond (e.g., \code{CovSnd}, \code{CovRec}, etc.); any other elements will be ignored.  The format of a given covariate element depends both on the effect type and on the number of covariates specified.  The basic cases are as follows:
\itemize{
\item Single covariate, time invariant: For \code{CovSnd}, \code{CovRec}, or \code{CovInt}, a vector or single-column matrix/array.  For \code{CovEvent}, an \code{n} by \code{n} matrix or array.
\item Multiple covariates, time invariant: For \code{CovSnd}, \code{CovRec}, or \code{CovInt}, a two-dimensional \code{n} by \code{p} matrix/array whose columns contain the respective covariates.  For \code{CovEvent}, a \code{p} by \code{n} by \code{n} array, whose first dimension indexes the covariate matrices.
\item Single or multiple covariates, time varying: For \code{CovSnd}, \code{CovRec}, or \code{CovInt}, an \code{m} by \code{p} by {n} array whose respective dimensions index time (i.e., event number), covariate, and actor.  For \code{CovEvent}, a \code{m} by \code{p} by \code{n} by \code{n} array, whose dimensions are analogous to the previous case.
}
Note that \dQuote{time varying} covariates may only change values when events transpire; thus, they should be regarded as temporally endogenous.  (See the reference below for details.) 

If called with \code{edgelist==NULL}, \code{rem.dyad} will produce a \dQuote{model skeleton} object containing the effects and other information, but no model fit.  (The seed coefficients, if given, are entered as the coefficients in the model, or else an uninteresting default set is used.)  The main purpose for this object is to set up an \emph{ab initio} simulation, as described below: once the skeleton is created, the \code{\link{simulate}} method can be used to generate draws from that model (without fitting to a data set).

A \code{\link{simulate}} method is provided for \code{rem.dyad} objects, which allows simulation of new event sequences from a fitted or skeleton model.  By default, a new sequence of length equal to the original sequence to which the model object was fitted is simulated (if applicable), but other lengths may be chosen using \code{nsim}.  Although the coefficients in the model object are used by default, this may also be altered by specifying \code{coef}.  Note that any covariates used must be passed to the simulate command via \code{covar} (using the same format as in the original model); this is in part because \code{rem.dyad} objects do not currently save their input data, and in part because dynamic covariates must always be the length of the simulated sequence (and hence must be factored when a non-default \code{nsim} value is used).  For models fit using \code{ordinal=TRUE}, the overall pacing of events will be arbitrary (more specifically, the simulation will tacitly assume that each event has a unit base hazard), but the relative timing is not.  See below for examples of both simulation using a fitted model object and \emph{ab initio} simulation without fitting a model to data.

For simulation, it is possible to fix the first portion of the event history by passing an event list matrix to the \code{edgelist} argument; this must be compatible with the target model (i.e., the vertex IDs must match), and it cannot contain \code{NA} values.  (Thus, if starting with an exact timing seqence with a last line containing \code{NA}s, this must be removed.)  If the input event list contains \code{m} events, then these are assumed to supply the first \code{m} events of the target sequence; if \code{m>nsim}, then any excess events are discarded.  By default, the input events are taken as fixed.  However, specifying \code{redraw.timing=TRUE} will lead the event timings to be redrawn, and \code{redraw.events} will lead the sender/reciver pairs to be redrawn.  This allows e.g. for an observed ordinal time sequence to be given a simulated exact time realization, by setting \code{nsim} to the event list length and setting \code{redraw.timing=TRUE}.  The more obvious use case is to simply extend an observed sequence, in which case one should use \code{nsim} greater than the input sequence length (i.e., the input length plus the number of new events to generate) and leave the \code{redraw} paraeters set to \code{FALSE}.  
}
\value{
For \code{rem.dyad}, an object of class \code{rem.dyad}.  For the \code{simulate} method, an event list.
%  ~Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
}
\references{ Butts, C.T.  (2008).  \dQuote{A Relational Event Framework for Social Action.}  \emph{Sociological Methodology}, 38(1). }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{rem}} }
\examples{
\dontrun{
#Generate some simple sample data based on fixed effects
roweff<-rnorm(10)                                     #Build rate matrix
roweff<-roweff-roweff[1]                   #Adjust for later convenience
coleff<-rnorm(10)
coleff<-coleff-coleff[1]
lambda<-exp(outer(roweff,coleff,"+"))
diag(lambda)<-0
ratesum<-sum(lambda)
esnd<-as.vector(row(lambda))                  #List of senders/receivers
erec<-as.vector(col(lambda))
time<-0
edgelist<-vector()
while(time<15){                   # Observe the system for 15 time units
  drawsr<-sample(1:100,1,prob=as.vector(lambda))        #Draw from model
  time<-time+rexp(1,ratesum)
  if(time<=15)                                             #Censor at 15
    edgelist<-rbind(edgelist,c(time,esnd[drawsr],erec[drawsr]))
  else
    edgelist<-rbind(edgelist,c(15,NA,NA))
}

#Fit the model, ordinal BPM
effects<-c("FESnd","FERec")
fit.ord<-rem.dyad(edgelist,10,effects=effects,hessian=TRUE)
summary(fit.ord)
par(mfrow=c(1,2))                                #Check the coefficients
plot(roweff[-1],fit.ord$coef[1:9],asp=1)
abline(0,1)
plot(coleff[-1],fit.ord$coef[10:18],asp=1)
abline(0,1)

#Now, find the temporal BPM
fit.time<-rem.dyad(edgelist,10,effects=effects,ordinal=FALSE,hessian=TRUE)
summary(fit.time)
plot(fit.ord$coef,fit.time$coef,asp=1)                  #Similar results
abline(0,1)

#Finally, try the BSIR method (note: a much larger expansion factor
#is recommended in practice)
fit.bsir<-rem.dyad(edgelist,10,effects=effects,fit.method="BSIR",
    sir.draws=100,sir.expand=5)
summary(fit.bsir)
par(mfrow=c(3,3))   #Examine the approximate posterior marginals
for(i in 1:9){
  hist(fit.bsir$post[,i],main=names(fit.bsir$coef)[i],prob=TRUE)
  abline(v=roweff[i+1],col=2,lwd=3)
}
for(i in 10:18){
  hist(fit.bsir$post[,i],main=names(fit.bsir$coef)[i],prob=TRUE)
  abline(v=coleff[i-8],col=2,lwd=3)
}

#Simulate an event sequence from the temporal model
sim<-simulate(fit.time,nsim=50000) #Simulate 50000 events
head(sim)                          #Show the event list
par(mfrow=c(1,2))                  #Check the behavior
esnd<-exp(c(0,fit.time$coef[1:9]))
esnd<-esnd/sum(esnd)*5e4           #Expected sending count
erec<-exp(c(0,fit.time$coef[10:18]))
erec<-erec/sum(erec)*5e4           #Expected sending count
plot(esnd,tabulate(sim[,2]),xlab="Expected Out-events",ylab="Out-events")
abline(0,1,col=2)
plot(erec,tabulate(sim[,3]),xlab="Expected In-events",ylab="In-events")
abline(0,1,col=2)

#Keep the first 10 events of the simulated sequence, and produce 10 more
sim.pre<-sim[1:10,]
sim2<-simulate(fit.time,nsim=20,edgelist=sim.pre)
sim.pre                            #See the first 10 events
sim2                               #First 10 events preserved
all(sim2[1:10,]==sim.pre)          #All TRUE

#Repeat, but redrawing part of the input sequence
sim2.t<-simulate(fit.time,nsim=20,edgelist=sim.pre,redraw.timing=TRUE)
sim2.e<-simulate(fit.time,nsim=20,edgelist=sim.pre,redraw.events=TRUE)
sim2.t                             #Events kept, timings not
sim2.t[1:10,]==sim.pre             #Second two columns TRUE
sim2.e                             #Timing kept, events not
sim2.e[1:10,]==sim.pre             #(Note: some events may repeat by chance!)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
