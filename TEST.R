library(devtools)
install_github("wpeterman/ResistanceOptimization")

require(ResistanceOptimization) # Installs package and the other required packages needed

Example.dir <- "C:/Example/"

# Prepare data for optimization function
Optim.input <-Optim.prep(
  Response=read.csv(paste0(Example.dir,"RESPONSE_mat.csv"),header=F),
  n.Pops=64,
  ASCII.dir=Example.dir,
  CS_Point.File=paste0(Example.dir,"samples64.txt"),
  CS.exe='"C:/Program Files/Circuitscape/4.0/Circuitscape/cs_run.exe"',
  Neighbor.Connect=8,
  Constrained.Max=10,
  Initial.shape=c(seq(0.1,1,0.4),seq(1.5,10,1.5)),
  Bootstrap=TRUE,
  boot.iters=100,
  Sample_Proportion=0.75)

# Run optimization function
Resistance.Optimization(Optim.input=Optim.input)

# Create Diagnostic plots
Diagnostic.Plots(resist.matrix.path=paste0(Optim.input$Results.dir,"Final_CS_Surfaces/Resist_surf1_resistances.out"),
  genetic.dist.vec=Optim.input$Response.vec,
  XLAB="Transformed resistance",
  YLAB="Pairwise distance",
  plot.dir=paste0(Optim.input$Results.dir,"Plots/"))
