\name{Optim.Boot}
\alias{Optim.Boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run bootsrap on optimized resistance surfaces}
\description{
This function will allow you to run the bootstrap implemented as part of `Resistance.Optimization` independently
}
\usage{
Optim.Boot(boot.iters = 10000, Optim.input, Sample_Proportion = 0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Optim.input}{
This is the object created from running 'Optim.prep'
}
 \item{boot.iters}{
Number of bootstrap iterations to conduct. Default = 10,000
}
 \item{Sample_Proportion}{
The proportion of samples to subsample in each bootstrap iteration. Default = 0.75}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function will generate two .csv files.
1) Boot_ModelSelection_Freq.csv will contain the summary of how often a given was the top model
2) Boot_AvgRank_summary.csv will contain the average model rank of each model
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Bill Peterman
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}



