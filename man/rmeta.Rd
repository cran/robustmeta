\name{rmeta}
\alias{rmeta}
\title{Robust estimation for meta-analysis with influential outlying studies}
\usage{
rmeta(y, v, model="RE", gamma=0.01)
}
\arguments{
\item{y}{A vector of the outcome measure estimates (e.g., MD, SMD, log OR, log RR, log HR, RD)}

\item{v}{A vector of the variance estimate of \code{y}}

\item{model}{Type of the pooling model; \code{"FE"}: Fixed-effect model or \code{"RE"}: Random-effects model; Default is \code{"RE"}}

\item{gamma}{Unit of grid search to explore the optimal value of tuning parameter \code{alpha} on (0,1); Default is 0.01}
}
\value{
Results of the robust inference for meta-analysis.
\itemize{
\item \code{mu}: Estimate of the common effect (for the fixed-effect model) or the grand mean (for the random-effects model).
\item \code{se}: Standard error estimate of \code{mu}.
\item \code{CI}: 95 percent confidence interval of \code{mu}.
\item \code{P}: P-value of the hypothesis test of \code{mu=0}.
\item \code{alpha}: Selected \code{alpha} by the Hyvarinen score.
\item \code{W}: Contribution rates of individual studies (\code{ui}: contribution rates of the conventional methods, \code{wi}: contribution rates of the robust methods). 
}
}
\description{
Implementing the robust inference for meta-analysis involving influential outlying studies based on the density power divergence.
}
\references{
Noma, H., Sugasawa, S. and Furukawa, T. A. (2022).
Robust inference methods for meta-analysis involving influential outlying studies. In Preparation.

Basu, A., Harris, I. R., Hjort, N. L., Jones, M. C. (1998).
Robust and efficient estimation by minimizing a density power divergence.
\emph{Biometrika.}
\strong{85}: 549-559.

Sugasawa, S. and Yonekura, S. (2021).
On selection criteria for the tuning parameter in robust divergence.
\emph{Entropy.}
\strong{23}: 1147.
}
\examples{
require(metafor)
data(clbp)
edat1 <- escalc(measure="SMD",m1i=m1,m2i=m2,sd1i=s1,sd2i=s2,n1i=n1,n2i=n2,data=clbp)
DL1 <- rma(yi, vi, data=edat1, method="DL")
print(DL1)         # ordinary DerSimonian-Laird method
plot(DL1)   # plots of influential statistics, etc.

###

y <- as.numeric(edat1$yi)		# definition of summary statistics
v <- edat1$vi

rmeta(y,v)                 # robust inference based on the random-effects model
rmeta(y,v,model="FE")      # robust inference based on the fixed-effect model
}
