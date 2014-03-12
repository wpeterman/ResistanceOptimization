\name{Diagnostic.Plots}
\alias{Diagnostic.Plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create diagnostic plots
}
\description{
This function will generate mixed effect model diagnostic plots following optimization with 'Resistance.Optimization'}
\usage{
Diagnostic.Plots(resist.matrix.path, genetic.dist.vec, XLAB = "Estimated resistance", YLAB = "Genetic distance", plot.dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resist.matrix.path}{
Provide the full path to the square distance matrix produced by Circuitscape (file ending with "_resistances.out")
Example: paste0(Optim.input$Results.dir,"Final_CS_Surfaces/XXX_resistances.out"
}
  \item{genetic.dist.vec}{
Vector of pairwise genetic distance (can supply as "Optim.input$Response.vec" after executing Optim.prep function)
}
  \item{XLAB}{
Label for x-axis (Defaults to "Estimated resistance")
}
  \item{YLAB}{
Label for y-axis (Defaults to "Genetic distance")
}
  \item{plot.dir}{
Directory to output PDF of diagnostic plots (can use "paste0(Optim.input$Results.dir,"Plots/"))" to output PDF in same directory as other response plots)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function will produce a single PDF called DiagnotisticPlots containing four plots.
1) The relationship between the optimized resistance value and your pairwise distance provided as a respsonse (genetic distance)
2) Scatterplot of the mixed effect model residuals in relation to resistance distance
3) Distribution of residuals
4) Q-Q plot 
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
# Diagnostic.Plots(resist.matrix.path=paste0(Optim.input$Results.dir,"Final_CS_Surfaces/resist1_resistances.out"),genetic.dist.vec=Optim.input$Response.vec,XLAB="Transformed resistance",YLAB="Pairwise distance",plot.dir=paste0(Optim.input$Results.dir,"Plots/"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

