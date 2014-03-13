######################################
# Make resposnse curves
#####################################
Response.Figs<- function(Optim.input){
  Top.params <- read.csv(file=paste0(Optim.input$Results.dir,"/TopModel_Optimization_Parameters.csv"),header=T)
  dir.create(file.path(Optim.input$Results.dir, "Plots"))
  
PDF.dir <- paste0(Optim.input$Results.dir,"Plots/")
for (i in 1:nrow(Top.params)){
  ASCII.file <- list.files(Optim.input$ASCII.dir,pattern=paste0(Top.params[i,1],".asc"),full.names=TRUE)
  if(Top.params[i,5]<1e-5 | Top.params[i,6]<1e-5) {
    cat(paste0("Plotting of ", Top.params[i,1]," could not be completed due to extremely small parameter estimates."),"\n")
    next # Use 'next' command to avoid testing invalid combinations    
  } else {
  PLOT.response(PARM=Top.params[i,c(5,6)],Resistance=ASCII.file,equation=Top.params[i,2],AIC=Top.params[i,7], OutputFolder=PDF.dir)
    }
  }
}
#######################################
# PLOT RESPONSE CURVES
#######################################
PLOT.response <- function(PARM,Resistance,equation,AIC, OutputFolder){
  PDF.dir=OutputFolder
  RAST <- raster(Resistance)
  Mn=cellStats(RAST,stat='min')
  Mx=cellStats(RAST,stat='max') 
  
  # Make vector of data
  dat.o <- seq(from=Mn,to=Mx,length.out=1000)
  dat.t <- SCALE.vector(data=dat.o,0,10)
  
  # Apply specified transformation
  if(equation=="Inverse-Reverse Monomolecular"){
    SIGN=-1
    Trans.vec <- SCALE.vector(SIGN*PARM[[2]]*(1-exp(dat.t/PARM[[1]]))+SIGN,MIN=1,MAX=PARM[[2]]+1) # Inverse-Reverse Monomolecular
    TITLE<-"Inverse-Reverse Monomolecular"
  } else if(equation=="Inverse Monomolecular"){
    SIGN=-1
    Trans.vec <- SCALE.vector(SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN,MIN=1,MAX=PARM[[2]]+1) # Inverse Monomolecular
    TITLE<-"Inverse Monomolecular"
  } else if(equation=="Monomolecular"){
    SIGN=1
    Trans.vec <- SCALE.vector(SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN,MIN=1,MAX=PARM[[2]]+1) # Monomolecular
     TITLE<-"Monomolecular"
 } else if(equation=="Reverse Monomolecular"){
    SIGN=1
    Trans.vec <- SIGN*PARM[[2]]*(1-exp(dat.t/PARM[[1]]))+SIGN # Reverse Monomolecular
    Trans.vec <- SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec)))
    TITLE<-"Reverse Monomolecular"

  } else if (equation=="Inverse Ricker") {
    SIGN=-1
    Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN # Inverse Ricker
    Trans.vec <- SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec)))
    TITLE<-"Inverse Ricker"
 }  else  {
    SIGN=1
    Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN #  Ricker
    TITLE<-"Ricker"
} 

  
  
  pdf(paste0(PDF.dir,"Optimized_Response_",RAST@data@names,".pdf"))
  plot(dat.o,Trans.vec,main=paste0("Optimization equation: " ,TITLE), xlab=paste0("Original data: ", RAST@data@names), ylab="Resistance value",type="l", lty=1, lwd=2)
  legend("right",legend=paste0("AIC = ",round(AIC,digits=3)),bty = "n")

dev.off()
  
}

#######################################
# Function to Loop through surfaces
######################################
Resistance.Optimization<-function(Optim.input){
  # Install necessary packages
  libs=c("raster", "lme4", "plyr")
  CheckInstallPackage(packages=libs)
  
  # Load libraries
  require(raster)
  require(lme4)
  require(plyr)
  #####################
  ASCII.files=Optim.input$ASCII.files
  ASCII.names=Optim.input$ASCII.names
  CS.exe=Optim.input$CS.exe
  CS_Point.File=Optim.input$CS_Point.File
  RESPONSE=Optim.input$Response.vec
  EXPORT.DIR=Optim.input$Write.dir
  EQUATION <- c("Ricker","Monomolecular") # Equations that can be run: Ricker, mono, logistic, exp, log, power
  DIRECTION <- c("neg", "pos") # pos, neg
  REVERSE <- c("true","false") # true, false
  count <- 0 # Counter to index iterations
  RESULTS <- list()  # Store optimization results in a list
for (i in 1:length(ASCII.files)){
  for (EQ in EQUATION){
    for (DIR in DIRECTION){
      for (REV in REVERSE){
        count<-count+1 # Count iterations
        EQUATION.name<-EQ.func2(reverse=REV,direction=DIR,EQ=EQ)
        
        if(EQUATION.name=="Reverse Ricker" | EQUATION.name=="Inverse-Reverse Ricker") {
          next # Use 'next' command to avoid testing invalid combinations    
        } else {
          
          RESISTANCE.SURFACE <- SCALE(raster(ASCII.files[i]),MIN=0,MAX=10) # Read in ASCII file, scale 0-10
          names(RESISTANCE.SURFACE)<-ASCII.names[i]
          
          # Smart start search of parameter space  
          ts <- SmartStart.Optimization_func(test.shape=Optim.input$Initial.shape,Resistance=RESISTANCE.SURFACE,equation=EQUATION.name, Max=Optim.input$Constrained.Max,Optim.input=Optim.input)
          
          # Get best random start values
          start.vals <- c(ts[[which.min(ts[,3]),1]],ts[[which.min(ts[,3]),2]])
          
          cat("\n", paste0("Optimizing both shape and maximum value ", count, "/",(length(ASCII.files)*length(EQUATION)*length(DIRECTION)*length(REVERSE))),paste0(": ",ASCII.names[i],"_",EQUATION.name), "\n")
          
          # Informed start values; these are the optimized values from the single parameter optimization
          shape_max.optim <-nlm(Resistance.Optimization_func, log(start.vals), Resistance=RESISTANCE.SURFACE, equation=EQUATION.name,get.best="false",Optim.input=Optim.input)
          
          results<-cbind(ASCII.names[i],EQUATION.name,t((start.vals)),t(exp(shape_max.optim$estimate)),shape_max.optim$minimum)
          colnames(results) <- c("Res.surf","Eq","shape.start", "max.start","Opt.Shape","Opt.Max","AIC")
          RESULTS[[count]] <- results
        } # Close "Skip Ricker" 
      } # Close reverse loop
    } # Close direction loop
  } # Close equation loop
} # Close ascii loop

# Fit distance only and null models
# Get Euclidean distances between points
ID<-Optim.input$ID
RESPONSE<-Optim.input$Response.vec
pC <-read.table(file=CS_Point.File,header=F)
pC<-as.matrix(pC[c("V2","V3")])
d<-pointDistance(pC,longlat=FALSE)
DIST <- d[lower.tri(d)]
dat<-cbind(ID,DIST,RESPONSE)

# Assign value to layer
# LAYER<-assign("LAYER",value=data$cs.matrix)
mod.dist <- lFormula(RESPONSE ~ DIST + (1|pop1), data=dat,REML=FALSE)
mod.dist$reTrms$Zt <- Optim.input$ZZ
dfun <- do.call(mkLmerDevfun,mod.dist)
opt <- optimizeLmer(dfun)
AIC.dist <- AIC(mkMerMod(environment(dfun), opt, mod.dist$reTrms,fr = mod.dist$fr))
Distance<-cbind("Distance","NA",AIC.dist)

# NULL
mod.null <- lFormula(RESPONSE ~ 1 + (1|pop1), data=dat,REML=FALSE)
mod.null$reTrms$Zt <- Optim.input$ZZ
dfun <- do.call(mkLmerDevfun,mod.null)
opt <- optimizeLmer(dfun)
AIC.null <- AIC(mkMerMod(environment(dfun), opt, mod.null$reTrms,fr = mod.null$fr))
Null<-cbind("Null","NA",AIC.null)

Dist_Null<-rbind(Distance,Null)
Dist_Null.df<-as.data.frame(Dist_Null);colnames(Dist_Null.df)<-c("Res.surf","Equation","AIC.stat")

# Process and export results
RESULTS.df<-as.data.frame(do.call(rbind,RESULTS))
# RESULTS.df<-as.data.frame(rbind(RESULTS.df,Dist_Null.df))
RESULTS.df$AIC<-as.numeric(as.character(RESULTS.df$AIC))
(RESULTS.df<-RESULTS.df[order(RESULTS.df$AIC,decreasing=F),])
write.table(RESULTS.df,file=paste0(Optim.input$Results.dir,"/All_Optimization_Results.csv"),sep=",",row.names=F,col.names=T)

#####################################
Optim_Results<-read.csv(paste0(Optim.input$Results.dir,"/All_Optimization_Results.csv"))
# Find optimized model based on AIC
Optimized.Model <- ddply(Optim_Results, .(Res.surf), summarise, Equation= Eq[which.min(AIC)],AIC.stat= min(AIC))
Optimized.Model<-rbind(Optimized.Model,Dist_Null.df)
Optimized.Model$AIC.stat<-as.numeric(as.character(Optimized.Model$AIC.stat))
Optimized.Model$AICc <- Optimized.Model$AIC.stat + ((2*2*(2+1))/(Optim.input$n.Pops-2-1))
Optimized.Model$delta.AICc<-abs(min(Optimized.Model$AICc)-Optimized.Model$AICc)
Optimized.Model$likelihood<-exp(-.5*Optimized.Model$delta.AICc)
Optimized.Model$weight<-Optimized.Model$likelihood/sum(Optimized.Model$likelihood)
(Optimized.Model_Full<-arrange(df=Optimized.Model,delta.AICc))
write.table(Optimized.Model_Full,file=paste0(Optim.input$Results.dir,"/TopModel_Optimization_Results.csv"),sep=",",row.names=F,col.names=T)
############################################
############################################
# Get parameters to run optimal models and then conduct bootstrap
MATCH<-Optimized.Model_Full$Res.surf!="Null"& Optimized.Model_Full$Res.surf!="Distance"
Top.Models <-match(Optimized.Model_Full$Res.surf[MATCH],Optim_Results$Res.surf)

# Write best model parameters to file
Top.params <-Optim_Results[Top.Models,]
write.table(Top.params,file=paste0(Optim.input$Results.dir,"/TopModel_Optimization_Parameters.csv"),sep=",",row.names=F,col.names=T)

  # Where will final CS results be written to?
  dir.create(file.path(Optim.input$Results.dir, "Final_CS_Surfaces"))
  Optim.input2<-Optim.input
  Optim.input2$Write.dir <-paste0(Optim.input$Results.dir,"/Final_CS_Surfaces/")  

# Run through each optimal surface and make final CS resistance matrices
for (i in 1:length(Top.Models)){
  param <-Optim_Results[Top.Models[i],c(1,2,5,6)]
  rast <-SCALE(raster(paste0(Optim.input$ASCII.dir,param[[1]],".asc")),0,10)
  names(rast)<-Optim.input$ASCII.names[i]
  Resistance.Optimization_func(PARM=log(c(param$Opt.Shape,param$Opt.Max)),Resistance=rast,equation=param$Eq,get.best="true",Optim.input=Optim.input2)
}
    # Boostrap
  if(Optim.input$Bootstrap==TRUE){
  Optim.Boot(boot.iters=Optim.input$boot.iters,   #  Number of boostrap iterations (Default = 10 000)
             Optim.input=Optim.input, # Optim input object created from running 'Optim.prep"
             Sample_Proportion=0.75 # Proportion of samples to be included in each bootstrap iteration (Default = 0.75)
  )
  }
  
  # Create response figures
  Response.Figs(Optim.input) # Input object created from running 'Optim.prep' in step 1 above
  
  # Generate parameter estimates
  Coeff.Table(resist.dir=paste0(Optim.input$Results.dir,"Final_CS_Surfaces/"),genetic.dist.vec=Optim.input$Response.vec,Optim.input=Optim.input)
  unlink(paste0(Optim.input$Results.dir,"tmp"),recursive=TRUE)
}

###############################################
# Bootstrap AIC Model Fits
###############################################
Optim.Boot<-function(boot.iters=10000, Optim.input,Sample_Proportion=0.75){
# Where are the CS results stored?
cs<-list.files(path=paste0(Optim.input$Results.dir,"/Final_CS_Surfaces/"),pattern= "resistances.out",full.names=TRUE)
NAMES<-gsub("_resistances.out","",x=list.files(path=paste0(Optim.input$Results.dir,"/Final_CS_Surfaces/"),pattern= "resistances.out"))

# What proportion of data should be sampled?
prop.sample <- Sample_Proportion
SAMPLE<-ceiling(Optim.input$n.Pops*prop.sample) # Number of samples to be drawn
# Read in response matrix
genetic.DIST.mat <- Optim.input$Response.mat

# Get Euclidean distances between points
pC <-read.table(file=Optim.input$CS_Point.File,header=F)
pC<-as.matrix(pC[c("V2","V3")])
d<-pointDistance(pC,longlat=FALSE)

# Number of bootstrap iterations to perform
iterations<-boot.iters
################################################
# Create population ID vectors
ID<-To.From.ID(POPS=SAMPLE)
ZZ<-ZZ.mat(ID)
###############################################
# Run bootstrap
BOOT.AIC <-list()
  progress_bar_text <- create_progress_bar("text")
  progress_bar_text$init(iterations)
for(iter in 1:iterations){
  progress_bar_text$step()
  
#   cat(paste0("Bootstrap iteration ",iter,"/",iterations, "\n")) 
  
  # BOOT.AIC<-foreach(iter=1:iterations,.packages=c("ecodist","lme4","plyr","reshape", "MuMIn")) %dopar% {
  # Array to store results in
  MLPE_Results<-matrix(nrow=length(NAMES),ncol=2);colnames(MLPE_Results)<-list("Model","AICc")
  
  SUBSET <- sample(1:Optim.input$n.Pops,size=SAMPLE,replace=FALSE)
  #   g.DIST<-lower(genetic.DIST.mat[SUBSET,SUBSET])  
  gdm<-genetic.DIST.mat[SUBSET,SUBSET]
  g.DIST<-gdm[lower.tri(gdm)]
  #   g.DIST<-1/(1-(gdm[lower.tri(gdm)]))
  dist<-d[SUBSET,SUBSET]
  DIST<-(dist[lower.tri(dist)])
  
  for (i in 1:length(NAMES)){
    cs.matrix<-read.matrix2(cs[i])
    Mod.Name<-NAMES[i]
    # Make random selection without replacement  
    cs.matrix<-cs.matrix[SUBSET,SUBSET]
    cs.matrix<-scale(as.numeric(cs.matrix[lower.tri(cs.matrix)]),center=TRUE,scale=TRUE)
    
    data<-cbind(ID,cs.matrix, g.DIST)
    
    # Assign value to layer
    LAYER<-assign("LAYER",value=data$cs.matrix)
    
    # Fit model
    mod <- lFormula(g.DIST ~ LAYER + (1|pop1), data=data,REML=FALSE)
    mod$reTrms$Zt <- ZZ
    dfun <- do.call(mkLmerDevfun,mod)
    opt <- optimizeLmer(dfun)
    AIC.stat <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
    
    results<-cbind(Mod.Name,AIC.stat)
    MLPE_Results[i,]<-results 
    
  }
  
  # Fit distance only and null models
  mod.dist <- lFormula(g.DIST ~ DIST + (1|pop1), data=data,REML=FALSE)
  mod.dist$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod.dist)
  opt <- optimizeLmer(dfun)
  AIC.dist <- AIC(mkMerMod(environment(dfun), opt, mod.dist$reTrms,fr = mod.dist$fr))
  Distance<-cbind("Distance",AIC.dist)
  
  # NULL
  mod.null <- lFormula(g.DIST ~ 1 + (1|pop1), data=data,REML=FALSE)
  mod.null$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod.null)
  opt <- optimizeLmer(dfun)
  AIC.null <- AIC(mkMerMod(environment(dfun), opt, mod.null$reTrms,fr = mod.null$fr))
  Null<-cbind("Null",AIC.null)
  
  MLPE_Results<-rbind(MLPE_Results,Distance,Null)
  MLPE_Results.df<-as.data.frame(MLPE_Results)
    
  MLPE_Results.df$AICc<-as.numeric(as.character(MLPE_Results.df$AICc))
  
  # Find optimized model based on chosen summary stat
  Optimized.Model <- ddply(MLPE_Results.df, .(Model), summarise, Model= Model[which.min(AICc)],AIC= min(AICc))
  Optimized.Model$AICc <- Optimized.Model$AIC + ((2*2*(2+1))/(Optim.input$n.Pops-2-1))
  Optimized.Model$delta.AICc<-abs(min(Optimized.Model$AICc)-Optimized.Model$AICc)
  Optimized.Model$likelihood<-exp(-.5*Optimized.Model$delta.AICc)
  Optimized.Model$weight<-Optimized.Model$likelihood/sum(Optimized.Model$likelihood)
  Optimized.Model_Full<-arrange(df=Optimized.Model,delta.AICc)
  BOOT.AIC[[iter]] <- Optimized.Model_Full
}


# Extract data from list
BOOT.BEST.summary.Fst<-ldply(BOOT.AIC, .fun=function(x) cbind(as.character(x$Model[1]),min(x$AICc),max(x$weight)))
colnames(BOOT.BEST.summary.Fst)<-c("Surface","AICc","weight")
# BOOT.BEST.summary

Freq.TopModel<-(table(BOOT.BEST.summary.Fst$Surface)/iterations)*100
Freq.TopModel

write.table(Freq.TopModel,file=paste0(Optim.input$Results.dir,"Boot_ModelSelection_Freq.csv"),sep=",",row.names=F,col.names=T)

# Get ranks of each layer for each iteration
BOOT.rank<-ldply(BOOT.AIC, .fun=function(x) cbind(as.character(x$Model),rank(x$AICc), x$weight))
colnames(BOOT.rank)<-c("Surface","Rank", "weight")
BOOT.rank$Rank<-as.numeric(as.character(BOOT.rank$Rank))
BOOT.rank$weight<-as.numeric(as.character(BOOT.rank$weight))
  BOOT.LAYER.summary<-ddply(BOOT.rank,.(Surface),summarise,Avg.Rank=mean(Rank), SD.Rank=sd(Rank), Avg.weight=mean(weight), SD.weight=sd(weight),LCI.weight=quantile(weight,probs=c(0.025)),UCI.weight=quantile(weight,probs=c(0.975)))
# BOOT.LAYER.summary
(BOOT.LAYER.summary<-BOOT.LAYER.summary[order(BOOT.LAYER.summary[,"Avg.Rank"],decreasing=F), ]) # Sort layers by rank

write.table(BOOT.LAYER.summary,file=paste0(Optim.input$Results.dir,"Boot_AvgRank_summary.csv"),sep=",",row.names=F,col.names=T)
}

#######################################
# ORIGINAL OPTIMIZATION FUNCTION
#######################################
Resistance.Optimization_func<-function(PARM,Resistance,equation, Optim.input,get.best) {
  t1<-Sys.time()
  
  #   File.name <- paste0("exp",TRAN,"_MAX", MAX)
  if (get.best=="false"){
    File.name <- "optim_iter"
    MAP="write_cum_cur_map_only = False"
    CURRENT.MAP="write_cur_maps = False"

  } else {
    File.name <- Resistance@data@names
    MAP="write_cum_cur_map_only = True"
    CURRENT.MAP="write_cur_maps = 1"
  }
  EXPORT.DIR<-Optim.input$Write.dir
  RESIST<-Resistance
#   SIGN <- ifelse(direction=="pos",1,-1)

  # Read in resistance surface to be optimized
  SHAPE <- exp(PARM[1])
  Max.SCALE <- exp(PARM[2])
  SHAPE<-ifelse(exp(PARM[1])>10000,10000,exp(PARM[1])) # Upper boundaries on parameters
  SHAPE<-ifelse(exp(PARM[1])<1e-4,1e-4,exp(PARM[1])) # Upper boundaries on parameters
  Max.SCALE<-ifelse(exp(PARM[2])>100e6,100e6,exp(PARM[2]))
  Max.SCALE<-ifelse(exp(PARM[2])<1e-6,1e-6,exp(PARM[2]))
  cat("\n", Resistance@data@names, as.character(equation),paste0("| Shape = ",SHAPE,"; "),paste0("Maximum scale = ",Max.SCALE),"\n")  
  
  # Apply specified transformation
  if(equation==1|equation=="Inverse-Reverse Monomolecular"){
    SIGN=-1
    R1 <- SIGN*Max.SCALE*(1-exp(RESIST/SHAPE))+SIGN # Inverse-Reverse Monomolecular
    RESIST.t <- SCALE(R1,MIN=abs(cellStats(R1,stat='max')),MAX=abs(cellStats(R1,stat='min')))
    #     EQ <- "Inverse-Reverse Monomolecular"
    
  } else if(equation==2|equation=="Reverse Monomolecular"){
    SIGN=1
    R1 <- SIGN*Max.SCALE*(1-exp(RESIST/SHAPE))+SIGN # Reverse Monomolecular
    RESIST.t <- SCALE(R1,MIN=abs(cellStats(R1,stat='max')),MAX=abs(cellStats(R1,stat='min')))
    #     EQ <- "Reverse Monomolecular"        
    
  } else if(equation==3|equation=="Monomolecular"){
    SIGN=1
    RESIST.t <- SIGN*Max.SCALE*(1-exp(-1*RESIST/SHAPE))+SIGN # Monomolecular
    #     EQ <- "Monomolecular"
    
  } else if (equation==4|equation=="Inverse Monomolecular") {
    SIGN=-1
    R1 <- SIGN*Max.SCALE*(1-exp(-1*RESIST/SHAPE))+SIGN # Inverse Monomolecular
    RESIST.t <- SCALE(R1,MIN=abs(cellStats(R1,stat='max')),MAX=abs(cellStats(R1,stat='min')))
    #     EQ <- "Inverse Monomolecular"  	    
    
  } else if (equation==5|equation=="Inverse Ricker") {
    SIGN=-1
    R1 <- SIGN*(Max.SCALE*RESIST*exp(-1*RESIST/SHAPE))+SIGN # Inverse Ricker
    RESIST.t <- SCALE(R1,MIN=abs(cellStats(R1,stat='max')),MAX=abs(cellStats(R1,stat='min')))
    #     EQ <- "Inverse Ricker"  
    
  } else if (equation==6|equation=="Ricker") {
    SIGN=1
    RESIST.t <- SIGN*(Max.SCALE*RESIST*exp(-1*RESIST/SHAPE))+SIGN #  Ricker
    #     EQ <- "Ricker"
    
  } else if (equation==7|equation=="Inverse Gaussian") {
    RESIST.t <- Max.SCALE - Max.SCALE*exp((-1*(RESIST-OPT)^2)/(2*SD^2))+SIGN #  Inverse Gaussian
    #     EQ <- "Inverse Gaussian"  
    
  } else if (equation==8|equation=="Gaussian") {
    RESIST.t <- Max.SCALE*exp((-1*(RESIST-OPT)^2)/(2*SD^2))+SIGN #  Gaussian
    #     EQ <- "Gaussian"      
  } else {
    RESIST.t <- (RESIST*0)+1 #  Distance
    #     EQ <- "Distance"    
  } # End if-else
    
  RESIST.t[is.na(RESIST.t[])]<-1 
    
  if(cellStats(RESIST.t,"max")>100e6)  RESIST<-SCALE(RESIST.t,1,100e6) # Rescale surface in case resistance are too high
  RESIST.t <- reclassify(RESIST.t, c(-Inf,1e-06, 1e-06,100e6,Inf,100e6))
  writeRaster(x=RESIST.t,filename=paste0(EXPORT.DIR,File.name,".asc"), overwrite=TRUE)
  
  
  # Write Circuitscape Batch files
  BATCH<-paste0(EXPORT.DIR,File.name,".ini")        
  OUT<-paste0(paste0("output_file = ",EXPORT.DIR), File.name,".out")
  HABITAT<-paste0("habitat_file = ",paste0(EXPORT.DIR,File.name,".asc"))
  LOCATION.FILE <- paste0("point_file = ", Optim.input$CS_Point.File)
  ifelse(Optim.input$Neighbor.Connect==4,connect<-"True",connect<-"False")
  CONNECTION=paste0("connect_four_neighbors_only=",connect)
  
#   if(CS.version=='3.5.8'){
#     write.CS_3.5.8(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,CONNECTION=CONNECTION)
#   } else {
    write.CS_4.0(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,CONNECTION=CONNECTION,CURRENT.MAP=CURRENT.MAP,MAP=MAP)    
#   }
##########################################################################################
  # Run Circuitscape
  CS.exe=Optim.input$CS.exe
  
  # Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
  hidden = TRUE
  
  # CS.Batch<- list.files(path=EXPORT.DIR, pattern = "\\.ini$") # Make list of all files with '.asc' extension
  CS.ini <- paste0(EXPORT.DIR,File.name,".ini")
  CS.Run.output<-system(paste(CS.exe, CS.ini), hidden) 
  
  #########################################
  # Run mixed effect model on each Circuitscape effective resistance
  
  CS.results<-paste0(EXPORT.DIR,File.name,"_resistances.out")
  
  # Get AIC statistic for transformed-scaled resistance surface
  cs.matrix<-scale(read.matrix(CS.results),center=TRUE,scale=TRUE)
  ID=Optim.input$ID
  RESPONSE=Optim.input$Response.vec
  data<-cbind(ID,cs.matrix,RESPONSE)
  
  # Assign value to layer
  LAYER<-assign("LAYER",value=data$cs.matrix)
  
  # Fit model
  mod <- lFormula(RESPONSE ~ LAYER + (1|pop1), data=data,REML=FALSE)
  mod$reTrms$Zt <-  Optim.input$ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  AIC.stat <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  
  t2 <-Sys.time()
  cat(paste0("\t", "Iteration took ", round(t2-t1,digits=2), " seconds to complete"),"\n")
  cat(paste0("\t", "AIC = ",round(AIC.stat,3)),"\n")
  
  AIC.stat # Function to be minimized    
  
}


#############################
# SMART START FUNCTION
#############################
SmartStart.Optimization_func<-function(test.shape, Resistance,equation,Max,Optim.input) {
  SmartStart_results <- matrix(nrow=length(test.shape),ncol=3); colnames(SmartStart_results)<-c("Shape", "Max", "AIC")
 EXPORT.DIR<-Optim.input$Write.dir
 RESIST<-Resistance
#   SIGN <- ifelse(direction=="pos",1,-1)
    
  for(i in 1:length(test.shape)){
    t1<-Sys.time()
    
    cat("\n", "Testing parameter space:",paste0("Shape = ",test.shape[i]),";",paste0("Maximum = ",Max), "\n")
    cat(paste0("\t" ,"Test iteration: ",i, "/",length(test.shape)), "--->",paste0(Resistance@data@names,"_",equation),"\n")
    
# Read in resistance surface to be optimized
SHAPE <- test.shape[i]
Max.SCALE <- Max

# Apply specified transformation
if(equation==1|equation=="Inverse-Reverse Monomolecular"){
  SIGN=-1
  R1 <- SIGN*Max.SCALE*(1-exp(RESIST/SHAPE))+SIGN # Inverse-Reverse Monomolecular
  RESIST.t <- SCALE(R1,MIN=abs(cellStats(R1,stat='max')),MAX=abs(cellStats(R1,stat='min')))
  #     EQ <- "Inverse-Reverse Monomolecular"
  
} else if(equation==2|equation=="Reverse Monomolecular"){
  SIGN=1
  R1 <- SIGN*Max.SCALE*(1-exp(RESIST/SHAPE))+SIGN # Reverse Monomolecular
  RESIST.t <- SCALE(R1,MIN=abs(cellStats(R1,stat='max')),MAX=abs(cellStats(R1,stat='min')))
  #     EQ <- "Reverse Monomolecular"        
  
} else if(equation==3|equation=="Monomolecular"){
  SIGN=1
  RESIST.t <- SIGN*Max.SCALE*(1-exp(-1*RESIST/SHAPE))+SIGN # Monomolecular
  #     EQ <- "Monomolecular"
  
} else if (equation==4|equation=="Inverse Monomolecular") {
  SIGN=-1
  R1 <- SIGN*Max.SCALE*(1-exp(-1*RESIST/SHAPE))+SIGN # Inverse Monomolecular
  RESIST.t <- SCALE(R1,MIN=abs(cellStats(R1,stat='max')),MAX=abs(cellStats(R1,stat='min')))
  #     EQ <- "Inverse Monomolecular"        
  
} else if (equation==5|equation=="Inverse Ricker") {
  SIGN=-1
  R1 <- SIGN*(Max.SCALE*RESIST*exp(-1*RESIST/SHAPE))+SIGN # Inverse Ricker
  RESIST.t <- SCALE(R1,MIN=abs(cellStats(R1,stat='max')),MAX=abs(cellStats(R1,stat='min')))
  #     EQ <- "Inverse Ricker"  
  
} else if (equation==6|equation=="Ricker") {
  SIGN=1
  RESIST.t <- SIGN*(Max.SCALE*RESIST*exp(-1*RESIST/SHAPE))+SIGN #  Ricker
  #     EQ <- "Ricker"
  
} else if (equation==7|equation=="Inverse Gaussian") {
  RESIST.t <- Max.SCALE - Max.SCALE*exp((-1*(RESIST-OPT)^2)/(2*SD^2))+SIGN #  Inverse Gaussian
  #     EQ <- "Inverse Gaussian"  
  
} else if (equation==8|equation=="Gaussian") {
  RESIST.t <- Max.SCALE*exp((-1*(RESIST-OPT)^2)/(2*SD^2))+SIGN #  Gaussian
  #     EQ <- "Gaussian"      
} else {
  RESIST.t <- (RESIST*0)+1 #  Distance
  #     EQ <- "Distance"    
} # End if-else

RESIST.t[is.na(RESIST.t[])]<-1 
    
    File.name <- "optim_iter"
    MAP<-"write_cum_cur_map_only = False"
    CURRENT.MAP<-"write_cur_maps = False"
  
    writeRaster(x=RESIST.t,filename=paste0(EXPORT.DIR,File.name,".asc"), overwrite=TRUE)    
    
    # Write Circuitscape Batch files
    BATCH<-paste0(EXPORT.DIR,File.name,".ini")        
    OUT<-paste0(paste0("output_file = ",EXPORT.DIR), File.name,".out")
    HABITAT<-paste0("habitat_file = ",paste0(EXPORT.DIR,File.name,".asc"))
    LOCATION.FILE <- paste0("point_file = ", Optim.input$CS_Point.File)
    ifelse(Optim.input$Neighbor.Connect==4,connect<-"True",connect<-"False")
    CONNECTION=paste0("connect_four_neighbors_only=",connect)
    
#     if(CS.version=='3.5.8'){
#       write.CS_3.5.8(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,CONNECTION=CONNECTION)
#     } else {
      write.CS_4.0(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,CONNECTION=CONNECTION,CURRENT.MAP=CURRENT.MAP,MAP=MAP)    
#     }
##########################################################################################
    # Run Circuitscape
  CS.exe=Optim.input$CS.exe
  
    # Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
    hidden = TRUE
    
    # CS.Batch<- list.files(path=EXPORT.DIR, pattern = "\\.ini$") # Make list of all files with '.asc' extension
    CS.ini <- paste0(EXPORT.DIR,File.name,".ini")
    CS.Run.output<-system(paste(CS.exe, CS.ini), hidden) 
    
    #########################################
    # Run mixed effect model on each Circuitscape effective resistance
    
    CS.results<-paste0(EXPORT.DIR,File.name,"_resistances.out")
    
    # Get AIC statistic for transformed-scaled resistance surface
    cs.matrix<-scale(read.matrix(CS.results),center=TRUE,scale=TRUE)
    ID=Optim.input$ID
    RESPONSE=Optim.input$Response.vec
    data<-cbind(ID,cs.matrix,RESPONSE)
    
    # Assign value to layer
    LAYER<-assign("LAYER",value=data$cs.matrix)
    
    # Fit model
    mod <- lFormula(RESPONSE ~ LAYER + (1|pop1), data=data,REML=FALSE)
    mod$reTrms$Zt <- Optim.input$ZZ
    dfun <- do.call(mkLmerDevfun,mod)
    opt <- optimizeLmer(dfun)
    AIC.stat <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
    
    t2 <-Sys.time()
    cat(paste0("\t", "Iteration took ", round(t2-t1,digits=2), " seconds to complete"),"\n")
    cat(paste0("\t", "AIC = ",round(AIC.stat,3)),"\n")
    
    SmartStart_results[i,] <- cbind(SHAPE,Max.SCALE,AIC.stat)
    #   AIC.stat # Function to be minimized    
    
    }
  (SmartStart_results)
  }

##############################################################
Diagnostic.Plots<-function(resist.matrix.path, genetic.dist.vec, XLAB="Estimated resistance",YLAB="Genetic distance",plot.dir){
  RESPONSE=genetic.dist.vec
  mm<-read.table(resist.matrix.path)[-1,-1]
  m<-length(mm)
  ID<-To.From.ID(POPS=m)
  ZZ<-ZZ.mat(ID=ID)
  cs.matrix<-scale(mm[lower.tri(mm)],center=TRUE,scale=TRUE)
  cs.unscale<-mm[lower.tri(mm)]
  dat<-cbind(ID,cs.matrix,RESPONSE)
  
  # Assign value to layer
  LAYER<-assign("LAYER",value=dat$cs.matrix)
  
  # Fit model
  mod <- lFormula(RESPONSE ~ LAYER + (1|pop1), data=dat,REML=TRUE)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  Mod <- (mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  #######
  # Make diagnostic plots
  #   par(mfrow=c(2,2))
  pdf(file = paste0(plot.dir,"DiagnosticPlots.pdf"))
  par(mfrow=c(2,2),
      oma = c(0,4,0,0) + 0.1,
      mar = c(4,4,1,1) + 0.1)
  plot(genetic.dist.vec~cs.unscale,xlab=XLAB,ylab=YLAB)
  abline(lm(genetic.dist.vec~cs.unscale))
  plot(residuals(Mod)~cs.unscale,xlab=XLAB,ylab="Residuals")
  abline(lm(residuals(Mod)~cs.unscale))
  hist(residuals(Mod),xlab="Residuals",main="")
  qqnorm(resid(Mod),main="")
  qqline(resid(Mod))
  dev.off()
  par(mfrow=c(1,1))  
}
##############################################################
# Run Mixed effects models, recovery parameter estimates
Coeff.Table <- function(resist.dir, genetic.dist.vec,Optim.input){ 
  RESPONSE=genetic.dist.vec
  resist.mat<-list.files(resist.dir,pattern="*_resistances.out",full.names=TRUE)
  resist.names<-gsub(pattern="_resistances.out","",x=list.files(resist.dir,pattern="*_resistances.out"))
  COEF.Table<-array()
  for(i in 1:length(resist.mat)){
    m<-length(read.table(resist.mat[i])[-1,-1])
    mm<-read.table(resist.mat[i])[-1,-1]
    ID<-To.From.ID(POPS=m)
    ZZ<-ZZ.mat(ID=ID)
    cs.matrix<-scale(mm[lower.tri(mm)],center=TRUE,scale=TRUE)
    
    dat<-cbind(ID,cs.matrix,RESPONSE)
    
    # Assign value to layer
    LAYER<-assign(resist.names[i],value=dat$cs.matrix)
    
    # Fit model
    mod <- lFormula(RESPONSE ~ LAYER + (1|pop1), data=dat,REML=TRUE)
    mod$reTrms$Zt <- ZZ
    dfun <- do.call(mkLmerDevfun,mod)
    opt <- optimizeLmer(dfun)
    Mod.Summary <- summary(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
    COEF<-Mod.Summary$coefficients
    row.names(COEF)<-c("Intercept",resist.names[i])
    COEF.Table<-rbind(COEF.Table, COEF)
  }
  COEF.Table<-COEF.Table[-1,]
  write.table(COEF.Table,file=paste0(Optim.input$Results.dir,"Coefficient_Table.csv"),sep=",",row.names=TRUE,col.names=NA)
}
##############################################################

##############################################################
############ OTHER NECESSARY FUNCTIONS  #####################
#############################################################

# FUNCTIONS
read.matrix<-function(cs.matrix){  m<-read.table(cs.matrix)[-1,-1]
                                   m[lower.tri(m)]}

read.matrix2<-function(cs.matrix){  m<-read.table(cs.matrix)[-1,-1]
                                   }
# Make to-from population list
To.From.ID<-function(POPS){
  tmp <- matrix(nrow=POPS,ncol=POPS)
  dimnames(tmp) <- list( 1:POPS, 1:POPS)  
  tmp2 <- as.data.frame( which( row(tmp) < col(tmp), arr.ind=TRUE))  
  tmp2[[2]] <-dimnames(tmp)[[2]][tmp2$col]
  tmp2[[1]] <-dimnames(tmp)[[2]][tmp2$row]
  colnames(tmp2)<-c("pop1","pop2")
  as.numeric(tmp2$pop1);as.numeric(tmp2$pop2)
  ID<-arrange(tmp2,as.numeric(pop1),as.numeric(pop2))
  p1<-ID[POPS-1,1]; p2<-ID[POPS-1,2]
  ID[POPS-1,1]<-p2; ID[POPS-1,2]<-p1
  ID$pop1 <- factor(ID$pop1)
  ID$pop2 <- factor(ID$pop2)
  return(ID)
}

# Create ZZ matrix for mixed effects model
ZZ.mat <- function(ID) {
Zl <- lapply(c("pop1","pop2"), function(nm) Matrix:::fac2sparse(ID[[nm]],"d", drop=FALSE))
ZZ <- Reduce("+", Zl[-1], Zl[[1]])
return(ZZ)
}

# Rescale function
SCALE.vector <-function(data,MIN,MAX){Mn=min(data)
                               Mx=max(data)
                               (MAX-MIN)/(Mx-Mn)*(data-Mx)+MAX}

# Define scaling function
# This will rescale from 1 to specified MAX
SCALE <-function(data,MIN,MAX){Mn=cellStats(data,stat='min')
                               Mx=cellStats(data,stat='max')
                               (MAX-MIN)/(Mx-Mn)*(data-Mx)+MAX}

# Set random start values
# This function will ensure a more equitable search of shape values less than and greater than one
rand.scale.func <-function(){
  x=rnorm(1,0,.25)
  ifelse(x<=0,runif(1, .05, .95),runif(1, 1.01, 5))
}

rand.start <- function(){
  x=rand.scale.func()
  list(TRAN=log(x),
       MAX=log(runif(1, 25, 500)))
}
 
##################
EQ.func <-function(reverse,direction,EQ){  if(reverse=="true" && direction=="neg") {
  EQ <- paste0("rev.inv_",EQ)  
} else if(reverse=="true" && direction=="pos"){
  EQ <- paste0("rev_",EQ)
} else if(reverse=="false" && direction=="pos"){
  EQ <- paste0(EQ)
} else {
  EQ <- paste0("inv_",EQ)
}
 (EQ)                                       
}


EQ.func2 <-function(reverse,direction,EQ){  if(reverse=="true" && direction=="neg") {
  EQ <- paste0("Inverse-Reverse ",EQ)  
} else if(reverse=="true" && direction=="pos"){
  EQ <- paste0("Reverse ",EQ)
} else if(reverse=="false" && direction=="pos"){
  EQ <- paste0(EQ)
} else {
  EQ <- paste0("Inverse ",EQ)
}
 (EQ)                                       
}

# Function to write .ini file for Circuitscape 
write.CS_3.5.8 <- function(BATCH,OUT,HABITAT,LOCATION.FILE,CONNECTION,MAP="write_cum_cur_map_only=False"){
sink(BATCH)
cat("[Options for advanced mode]
ground_file_is_resistances=True
source_file=(Browseforacurrentsourcefile)
remove_src_or_gnd=keepall
ground_file=(Browseforagroundpointfile)
use_unit_currents=False
use_direct_grounds=False

[Calculation options]
low_memory_mode=False
solver=cg+amg
print_timings=True

[Options for pairwise and one-to-all and all-to-one modes]
included_pairs_file=None
point_file_contains_polygons=False
use_included_pairs=False")
cat("\n")
cat(LOCATION.FILE)
cat("\n")

cat("
[Output options]")
cat("\n")
cat(MAP)
cat("\n")
cat("log_transform_maps=False
set_focal_node_currents_to_zero=False
write_max_cur_maps=False
write_volt_maps=False
set_null_currents_to_nodata=True
set_null_voltages_to_nodata=True
compress_grids=False
write_cur_maps=False")
cat("\n")
cat(OUT)
cat("\n")
cat("\n")
cat("[Shortcircuit regions(aka polygons)]
use_polygons=False
polygon_file=(Browse for a short-circuit region file)

[Connection scheme for raster habitat data]")
cat(CONNECTION)
cat("\n")
cat("connect_using_avg_resistances=True

[Habitat raster or graph]
habitat_map_is_resistances=True")
cat("\n")
cat(HABITAT)
cat("\n")
cat("\n")
cat("[Options for one-to-all and all-to-one modes]
use_variable_source_strengths=False
variable_source_file=None

[Version]
version = 3.5.8

[Maskfile]
use_mask=False
mask_file=None

[Circuitscape mode]
data_type=raster
scenario=pairwise")
sink()
}


write.CS_4.0 <- function(BATCH,OUT,HABITAT,LOCATION.FILE,CONNECTION,CURRENT.MAP="write_cur_maps = False"
,MAP="write_cum_cur_map_only = False",PARALLELIZE="parallelize = False",CORES="max_parallel = 0"){
sink(BATCH)
cat("[Options for advanced mode]
ground_file_is_resistances = True
remove_src_or_gnd = rmvsrc
ground_file = (Browse for a raster mask file)
use_unit_currents = False
source_file = (Browse for a raster mask file)
use_direct_grounds = False

[Mask file]
mask_file = (Browse for a raster mask file)
use_mask = False

[Calculation options]
low_memory_mode = False")
cat("\n")
cat(PARALLELIZE)
cat("\n")
cat(CORES)
cat("\n")
cat("
solver = cg+amg
print_timings = False
preemptive_memory_release = False
print_rusages = False

[Short circuit regions (aka polygons)]
polygon_file = (Browse for a short-circuit region file)
use_polygons = False

[Options for one-to-all and all-to-one modes]
use_variable_source_strengths = False
variable_source_file = (Browse for a short-circuit region file)

[Output options]
set_null_currents_to_nodata = True
set_focal_node_currents_to_zero = False
set_null_voltages_to_nodata = True
compress_grids = False
write_volt_maps = False")
cat("\n")
cat(CURRENT.MAP)
cat("\n")
cat(OUT)
cat("\n")
cat(MAP)
cat("\n")
cat("
log_transform_maps = False
write_max_cur_maps = False

[Version]
version = 4.0-beta

[Options for reclassification of habitat data]
reclass_file = (Browse for file with reclassification data)
use_reclass_table = False

[Logging Options]
log_level = INFO
log_file = None
profiler_log_file = None
screenprint_log = False

[Options for pairwise and one-to-all and all-to-one modes]
included_pairs_file = (Browse for a file with pairs to include or exclude)
use_included_pairs = False")
cat("\n")
cat(LOCATION.FILE)
cat("\n")
cat("\n")

cat("
[Connection scheme for raster habitat data]
connect_using_avg_resistances = True")
cat(CONNECTION)
cat("\n")
cat("\n")
cat("
[Habitat raster or graph]
habitat_map_is_resistances = True")
cat("\n")
cat(HABITAT)
cat("\n")
cat("\n")
cat("
[Circuitscape mode]
data_type = raster
scenario = pairwise")
sink()
}

# Function to install and update necessary R packages
package <- function(pkgs, install=TRUE, update=FALSE, quiet=TRUE, verbose=TRUE, ...) {
  myrequire <- function(package, ...) {
    result <- FALSE
    if(quiet) { 
      suppressMessages(suppressWarnings(result <- require(package, ...)))
    } else {
      result <- suppressWarnings(require(package, ...))
    }
    return(result)
  }
  mymessage <- function(msg) {
    if(verbose) {
      message(msg)
    }
  }
  
  installedpkgs <- installed.packages()
  availpkgs <- available.packages(...)[,c('Package','Version')]
  if(nrow(availpkgs) == 0) {
    warning(paste0('There appear to be no packages available from the ',
      'repositories. Perhaps you are not connected to the ',
      'Internet?'))
  }
  # It appears that hyphens (-) will be replaced with dots (.) in version
  # numbers by the packageVersion function
  availpkgs[,'Version'] <- gsub('-', '.', availpkgs[,'Version'])
  results <- data.frame(loaded=rep(FALSE, length(pkgs)),
    installed=rep(FALSE, length(pkgs)),
    loaded.version=rep(as.character(NA), length(pkgs)),
    available.version=rep(as.character(NA), length(pkgs)),
    stringsAsFactors=FALSE)
  row.names(results) <- pkgs
  for(i in pkgs) {
    needInstall <- FALSE
    if(i %in% row.names(installedpkgs)) {
      v <- as.character(packageVersion(i))
      if(i %in% row.names(availpkgs)) {
        if(v != availpkgs[i,'Version']) {
          if(!update) {
            mymessage(paste0('A newer version of ', i, 
              ' is available ', '(current=', v, 
              '; available=',
              availpkgs[i,'Version'], ')'))
          }
          needInstall <- update
        }
        results[i,]$available.version <- availpkgs[i,'Version']
      } else {
        mymessage(paste0(i, ' is not available on the repositories.'))
      }
    } else {
      if(i %in% row.names(availpkgs)) {
        needInstall <- TRUE & install
        results[i,]$available.version <- availpkgs[i,'Version']
      } else {
        warning(paste0(i, ' is not available on the repositories and ',
          'is not installed locally'))
      }
    }
    if(needInstall | !myrequire(i, character.only=TRUE, ...)) {
      install.packages(pkgs=i, quiet=quiet, ...)
      if(!myrequire(i, character.only=TRUE, ...)) {
        warning(paste0('Error loading package: ', i))
      } else {
        results[i,]$installed <- TRUE
        results[i,]$loaded <- TRUE
        results[i,]$loaded.version <- as.character(packageVersion(i))
      }
    } else {
      results[i,]$loaded <- TRUE
      results[i,]$loaded.version <- as.character(packageVersion(i))
    }
  }
  if(verbose) {
    return(results)
  } else {
    invisible(results)
  }
}

##################
# Optimiazation preparation
Optim.prep<-function(Response,n.Pops,ASCII.dir,CS_Point.File,CS.exe,Neighbor.Connect=8,Constrained.Max=100,Initial.shape=c(seq(0.2,1,by=0.2),seq(1.25,10.75,by=0.75)),Bootstrap=FALSE,boot.iters=10000,Sample_Proportion=0.75){
  # Install necessary packages
  libs=c("raster", "lme4", "plyr")
  CheckInstallPackage(packages=libs)
  
  # Load libraries
  require(raster)
  require(lme4)
  require(plyr)
  #####################
  if(is.vector(Response)==TRUE || dim(Response)[1]!=dim(Response)[2]) {warning("Must provide square distance matrix with no column or row names")}
  Response.vec<-Response[lower.tri(Response)]
  ID <- To.From.ID(n.Pops)
  ZZ<-ZZ.mat(ID)
  ASCII.files<-list.files(ASCII.dir,pattern="*.asc",full.names=TRUE)
  ASCII.names<-gsub(pattern="*.asc","",x=(list.files(ASCII.dir,pattern="*.asc")))
  dir.create(file.path(ASCII.dir, "Results"))
  Results.dir<-paste0(ASCII.dir, "/Results/")
  dir.create(file.path(Results.dir, "tmp"))
  Write.dir <-paste0(Results.dir,"/tmp/")  
  
  list(Response.vec=Response.vec,Response.mat=Response,n.Pops=n.Pops,ID=ID,ZZ=ZZ,ASCII.files=ASCII.files,ASCII.names=ASCII.names,ASCII.dir=ASCII.dir,Write.dir=Write.dir,Results.dir=Results.dir,CS_Point.File=CS_Point.File,CS.exe=CS.exe,Neighbor.Connect=Neighbor.Connect,Constrained.Max=Constrained.Max,Initial.shape=Initial.shape,Bootstrap=Bootstrap,boot.iters=boot.iters,Sample_Proportion=Sample_Proportion)
}

#########################################
# Install necessary packages
CheckInstallPackage <- function(packages, repos="http://cran.r-project.org") {
  installed=as.data.frame(installed.packages())
  for(p in packages) {
    if(is.na(charmatch(p, installed[,1]))) { 
      install.packages(p, repos=repos) 
    }
  }
} 
##########################################################################################
