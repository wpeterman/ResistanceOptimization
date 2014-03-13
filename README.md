ResistanceOptimization
======================

R package to optimize continuous resistance surfaces

To install this package, execute the following commands in R:

install.packages("devtools") # Installs the 'devtools' package
library(devtools) # Loads devtools

install_github("wpeterman/ResistanceOptimization") # Downloads package
require(ResistanceOptimization) # Installs package and the other required packages needed
###########################################
In order to use this package, you must have CIRCUITSCAPE installed.
The 4.0-Beta release (or 4.0 version, when released) are required:
Beta: https://github.com/Circuitscape/Circuitscape/releases
Official releases: https://code.google.com/p/circuitscape/downloads/list

You can download example data by going to:
https://github.com/wpeterman/ResistanceOptimization/tree/master/data

Right click on ‘Example.zip’, 'Save link as...' to download the .zip file. Extract the zip file.

If all files are extracted to a folder called 'Example' located at 'C:/Example/', and the Circuitscape executable file (cs_run.exe) is located in '"C:/Program Files/Circuitscape/cs_run.exe"', the following lines of code would execute the functions.
############################################

Example.dir <- "C:/Example/"

# Prepare data for optimization function
Optim.input <-Optim.prep(
Response=read.csv(paste0(Example.dir,"RESPONSE_mat.csv"),header=F),
n.Pops=64,
ASCII.dir=Example.dir,
CS_Point.File=paste0(Example.dir,"samples64.txt"),
CS.exe='"C:/Program Files/Circuitscape/cs_run.exe"',
Neighbor.Connect=8,
Constrained.Max=100,
Initial.shape=c(seq(0.1,1,0.4),seq(1.5,10,1.5)),
Bootstrap=TRUE,
boot.iters=10000,
Sample_Proportion=0.75)

# Run optimization function
Resistance.Optimization(Optim.input=Optim.input)

# Create Diagnostic plots
Diagnostic.Plots(resist.matrix.path=paste0(Optim.input$Results.dir,"Final_CS_Surfaces/Resist_surf1_resistances.out
"),
genetic.dist.vec=Optim.input$Response.vec,
XLAB="Transformed resistance",
YLAB="Pairwise distance",
plot.dir=paste0(Optim.input$Results.dir,"Plots/"))