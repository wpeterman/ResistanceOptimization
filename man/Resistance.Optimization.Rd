\name{Resistance.Optimization}
\alias{Resistance.Optimization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Primary function to optimize resistance surfaces
}
\description{
This is the main optimization function. Simply provide the object created from running the 'Optim.prep' function. This function will optimize each resistance surface by first implementing a constrained grid search, and then using nlm.
}
\usage{
Resistance.Optimization(Optim.input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Optim.input}{
This is the object created from running 'Optim.prep'
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function will generate numerous files.
In the results directory, 6 .csv files will be created (assuming that bootstrapping was conducted). 
1) All_Optimization_Results gives the model transformation parameters and AIC score for every permutation tested. It includes the constrained search parameters (shape.start and max.start) as well as the final optimized values (Opt.Shape and Opt.Max)
2) TopModel_Optimization_Parameters gives only the parameters for the transformation of each surface with the lowest AIC. 
3) TopModel_Optimization_Results gives the AIC table for the optimized resistance surfaces
4) Coefficient_Table contains the linear mixed effects model parameter coefficient estimates
5) Boot_AvgRank_summary gives the average model rank, sd of model rank, avg weight, sd of weight, and upper and lower confidence intervals from the bootstrap iterations
6) Boot_ModelSelection_Freq provides the percent frequency that each resistance surface was selected as the top model

In the Plots folder, there is a PDF file for each optimized surface showing the transformed response curve.

In the Final_CS_Surfaces folder, are the .asc files and Circuitscape output for the final optimization of each surface. The cumulative resistance .asc file is also generated.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Bill Peterman
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Resistance.Optimization(Optim.input=Optim.input)
}


