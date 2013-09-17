\name{BADER}
\alias{BADER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Analysis of RNA Sequencing Data
}
\description{
This function estimates the posterior distribution of various parameters regarding RNA Sequencing data.
The most interesting parameter is the probability of differential expression (DE) between two groups
A and B. But also estimates for the log mean and the log dispersion parameter of the underlying
poisson - log-normal model can be returned.
}
\usage{
BADER(x, design, sizeFactors = TRUE, start = NULL, burn = 1000, reps = 10000, printEvery = 100, saveEvery = 1, t0 = 10, mode = "minimal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{m x n matrix: Every column should contain count data for a subject with m genes or tags.
}
  \item{design}{Factor specifying the samples' treatment groups.
    The first level of 'design' corresponds to the treatment group A,
    the second level to treatment group B }
  \item{sizeFactors}{
boolean: Whether size factors should be estimated (TRUE) or all set to 1 (FALSE)
}
  \item{start}{
list containing start values for MCMC sampler
}
  \item{burn}{Number of burning in steps}
  \item{reps}{Number of repetions}
  \item{printEvery}{After every \code{printEvery} iteration the current step is being printed}
  \item{saveEvery}{Every \code{saveEvery}-th step is saved for inference}
  \item{t0}{Warming up time for Metropolis-Hastings}
  \item{mode}{How much data should be returned? Returning all posterior distributions requires
        large memory.
        \itemize{
        \item \code{mode} = 0: Only posterior means are returned for every parameter
        \item \code{mode} = 1: Full posterior distribution for log fold change parameter is returned
        \item \code{mode} = 2: Posterior distibutions for the following parameters are returned: log fold change
        log mean and log dispersion
        }
  }
}
\value{
A list with posterior distributions / posterior means
}

\author{
Andreas Neudecker
}
\examples{
set.seed(21)

## log mean expression
muA <- rnorm(100,4,1)
gam <- c(rnorm(10,0,2),rep(0,90))
muB <- muA + gam

## log dispersion
alphaA <- alphaB <- rnorm(100,-2,1)

## count tables for treatment group a and b
kA <- t(matrix(rnbinom(300,mu=exp(muA),size=exp(-alphaA)),nrow=3,byrow=TRUE))
kB <- t(matrix(rnbinom(300,mu=exp(muB),size=exp(-alphaB)),nrow=3,byrow=TRUE))

x <- cbind(kA,kB)
design <- factor(c("A","A","A","B","B","B"))

results <- BADER(x,design,burn=1000,reps=2000)

\dontrun{
plot(results$diffProb,xlab="Index",ylab="posterior DE prob.")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
