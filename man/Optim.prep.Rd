\name{Optim.prep}
\alias{Optim.prep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare elements for optimization
}
\description{
This function prepares and compiles the necessary components to execute subsequent optimization functions, and must be run first. To run the example data, download the data from GitHub (url)
}
\usage{
Optim.prep(Response, n.Pops, ASCII.dir, CS_Point.File, CS.exe, Neighbor.Connect, Constrained.Max, Initial.shape, Bootstrap, boot.iters, Sample_Proportion)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Response}{
This must be a square distance matrix without row or column names
}
  \item{n.Pops}{
The number of populations in your study
}
  \item{ASCII.dir}{
The directory containing the continuous ASCII (.asc) files to be optimized. It is recommended that these files be stored their own directory.
}
  \item{CS_Point.File}{
The path to the Circuitscape formatted point file. See Circuitscape documentation for help.
}
  \item{CS.exe}{
The path to the CIRCUITSCAPE executable file (cs_run.exe). 
For example: 
'"C:/Program Files/Circuitscape/cs_run.exe"'

***NOTE: Double quotation used***
}
  \item{Neighbor.Connect}{
Select 4 or 8 to designate the connection scheme to use in CIRCUITSCAPE (Default = 8)
}
  \item{Constrained.Max}{
The fixed maximum resistance value to use during the initial constrained grid search (Default = 100)
}
  \item{Initial.shape}{
The initial shape parameters to test in constrained grid search (Default = c(seq(0.2,1,by=0.2),seq(1.25,10.75,by=0.75)))
}
  \item{Bootstrap}{
Logical,if TRUE, specify the number of bootstrap iterations for 'boot.iters'
}
  \item{boot.iters}{
Number of bootstrap iterations (Default = 10,000)
}
  \item{Sample_Proportion}{
Proportion of samples to be included in each bootstrap iteration (Default = 0.75)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Results in an object with the necessary components to run other the other functions
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Example.dir <- "C:/Example/"

# Optim.input<-Optim.prep(Response=read.csv(paste0(Example.dir,"RESPONSE_mat.csv"),header=F),n.Pops=64,ASCII.dir=Example.dir,CS_Point.File=paste0(Example.dir,"samples64.txt"),CS.exe='"C:/Program Files/Circuitscape/4.0/Circuitscape/cs_run.exe"',Neighbor.Connect=8,Results.dir=paste0(Example.dir,"Results/"),Constrained.Max=100,Initial.shape=c(seq(0.1,1,0.4),seq(1.5,10,1.5)),Bootstrap=TRUE,boot.iters=10000,Sample_Proportion=0.75)
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
