ResistanceOptimization
======================

R package to optimize continuous resistance surfaces

To install this package, execute the following commands in R:

install.packages("devtools") # Installs the 'devtools' package
library(devtools) # Loads devtools

install_github("wpeterman/ResistanceOptimization")
require(ResistanceOptimization)

In order to use this package, you must have CIRCUITSCAPE installed
I recommend the 4.0 Beta release:
https://github.com/Circuitscape/Circuitscape/releases

You can download example data by going to:
https://github.com/wpeterman/ResistanceOptimization/tree/master/data

Right click on each file, 'Save link as...' and select a folder to save each file in.

If all files were saved to a folder called 'Example' located at 'C:/Example/', and the Circuitscape executable file (cs_run.exe) is located in '"C:/Program Files/Circuitscape/cs_run.exe"', the following lines of code would execute the functions.

Example.dir <- "C:/Example/"

# Prepare data for optimization function
Optim.input <-Optim.prep(
Response=read.csv(paste0(Example.dir,"RESPONSE_mat.csv"),header=F),
n.Pops=64,
ASCII.dir=Example.dir,
CS_Point.File=paste0(Example.dir,"samples64.txt"),
CS.exe='"C:/Program Files/Circuitscape/4.0/Circuitscape/cs_run.exe"',
Neighbor.Connect=8,Results.dir=paste0(Example.dir,"Results/"),
Constrained.Max=100,
Initial.shape=c(seq(0.1,1,0.4),seq(1.5,10,1.5)),
Bootstrap=TRUE,
boot.iters=10000,
Sample_Proportion=0.75)

# Run optimization function
Resistance.Optimization(Optim.input=Optim.input)

# Create Diagnostic plots
Diagnostic.Plots(resist.matrix.path=paste0(Optim.input$Results.dir,"Final_CS_Surfaces/resist1_resistances.out"),
genetic.dist.vec=Optim.input$Response.vec,
XLAB="Transformed resistance",
YLAB="Pairwise distance",
plot.dir=paste0(Optim.input$Results.dir,"Plots/"))