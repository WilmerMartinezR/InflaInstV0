

######################################################################
### Funcion Kernel
######################################################################

KernelA <- function(a,TT){
  tau <- 0:(TT-1)
  nuM <- (TT - tau)^a
  return(kapa = nuM*TT/sum(nuM))
}
# KernelA(4,12)

#Funcion inflacion instantanea
InflaIns <- function(x,kapa){
  return(prod(x^kapa)-1)
}

GuardaOut <- function(FechaFinAll, NombreOut, Grupo, BaseG, path.out){
  IniF <- FechaFinAll[1]
  FinF <- FechaFinAll[length(FechaFinAll)]
  file.out <- paste(NombreOut, "_", Grupo, "_",IniF,"_",FinF, ".txt",sep = "")
  write.table(BaseG, paste(path.out, file.out, sep = "") 
              ,sep = "\t")
}

## Pre-whiting function
# This return the optimum lag
LagOpt_ccf <- function(x,y){
  ModX <- auto.arima(x)
  pwx <- ModX$residuals
  newpwy <- residuals(Arima(y, model = ModX))
  CrossC <- c(ccf(pwx,newpwy,na.action=na.omit))
  CrossC.1 <- data.frame(Lag = CrossC$lag, ccf = CrossC$acf)
  CrossC.2 <- CrossC.1[which.max(abs(CrossC.1$ccf)),]
  return(CrossC.2$Lag)
}

###########################################################################
RECM <- function(x) sqrt(mean(x^2,na.rm=T))
###########################################################################

# The data has to have these columns Fecha, model, pronos and obs per horizon.
Eval_DM <- function(DCCP, Fecha="Fecha",Model="model",
                    IniObs = "obs", IniP = "pron", alpha = 0.05){
  
  ## Calculo de los errores por modelo
  tvar0 <- paste("^",IniObs,sep=""); name0 <- c(grep(tvar0, names(DCCP)))
  tvar1 <- paste("^",IniP,sep=""); name1 <- c(grep(tvar1, names(DCCP)))
  name0 <- name0[1:length(name1)] 
  FNC <- which(colnames(DCCP) == Fecha)
  mNC <- which(colnames(DCCP) == Model)
  
  DCCP.E <- data.frame(DCCP[,c(FNC,mNC)], DCCP[,name0] - DCCP[,name1])
  colnames(DCCP.E) <- c(colnames(DCCP.E)[1:2], paste("E",1:((ncol(DCCP.E)-2)), sep="_"))
  
  
  ## Calcula los RMSE por modelo
  Modelos <- unique(DCCP.E$model)
  RMSE <- NULL
  for(i in 1:length(Modelos)){
    #i=1
    DCCP.2 <- DCCP.E[DCCP.E$model == Modelos[i],]
    DCCP.3 <- DCCP.2[nrow(DCCP.2),]
    tvar1 <- "^E"; name1 <- c(grep(tvar1, names(DCCP.2)))
    DCCP.3[,name1] <- apply(DCCP.2[,name1], 2, RECM)
    RMSE <- rbind(RMSE, DCCP.3[-1])
  }
  
  #### Mejor modelo segun RECM
  tvar1 <- "^E"; name1 <- c(grep(tvar1, names(RMSE)))
  MinRECM <- apply(RMSE[,name1], 2, which.min)
  #RMSE[MinRECM[6],c("model",names(MinRECM)[6])]
  ##################################################################
  ###   Prueba DM
  ##################################################################
  
  tvar1 <- "^E"; name1 <- c(grep(tvar1, names(DCCP.E)))
  
  Modelo <- unique(DCCP.E$model)

  DM <- NULL
  for(i in 1:length(Modelo)){
    # i=91 # itera modelo
    PvalDM <- c()
    for(h in 1:(ncol(DCCP.E)-2)){
      #h=1 # itera horizonte
      mejor <- RMSE[MinRECM[h], ] 
      xx <- DCCP.E[DCCP.E$model == mejor$model, name1[h]]
      yy <- DCCP.E[DCCP.E$model == Modelo[i], name1[h]]
      
      if(mejor$model != Modelo[i]){
        Pval0 <- try(dm.test(na.omit(xx), 
                        na.omit(yy), 
                        h=h, alternative ="two.sided",
                        varestimator = c("bartlett")),TRUE)
        if(class(Pval0)[1] == "try-error") Pval <- 1 else Pval <- Pval0$p.value
                    
      }else Pval <- 1 
      PvalDM <- c(PvalDM, Pval) 
    }
    
    DM0 <- data.frame(RMSE[i,1], t(PvalDM))
    colnames(DM0) <- c(colnames(RMSE)[1], paste("DM",1:(ncol(DCCP.E)-2), sep="_"))
    DM <- rbind(DM,  DM0) 
  }
  
  DM.Filtra <- apply(DM[,2:ncol(DM)], 2, function(x) ifelse(x > alpha, 1,0))
  colnames(DM.Filtra) <- paste("Filtra",1:(ncol(DM.Filtra)), sep="_")
  DM <- cbind(DM, DM.Filtra)
  
  tvar0 <- "^Filtra"; name10 <- c(grep(tvar0, names(DM)))
  ReportF0 <- apply(DM[,name10], 2, sum)
  ReportF <- apply(DM[,name10], 2, sum)/nrow(DM)*100
  
  ### Crea matrix con el listado de los modelos equivalentes, de acuerdo con DM test
  ###  La primera fila corresponde al porcentaje de modelos equivalentes,
  ###   la segunda fila corresponde al modelo con el min RMSE.
  ###  Cada columna representa el horizonte de pronostico
  
  BMods <- data.frame(matrix(NA,nrow = max(ReportF0), ncol=length(name10)))
  BMods[1,] <- round(ReportF,1)
  BMods[2,] <- DM[MinRECM, Model]
  for(j in 1:length(name10)){
    #j=1
    equiMod <- DM[ DM[, c(1,name10[j])][,2] == 1, 1] # Filtra los modelos equivalentes
    Bmin <- which(equiMod == BMods[2,j])
    if(length(equiMod) > 1) BMods[3:(length(equiMod)+1),j] <- equiMod[-Bmin]
  }
  colnames(BMods) <- paste("h",1:length(name10), sep = "")
  
  return(list(MERROR = DCCP.E, RMSE = RMSE, DM = DM, BestMods = BMods))
  
}

###########################################################################

# The data has to have these columns Fecha, model, pronos and obs per horizon.
Eval_DM_V2 <- function(DCCP, Fecha="Fecha",Model="model",
                    IniObs = "obs", IniP = "pron", alpha = 0.05,ModRef = "Total_BAR."){
  
  ## Calculo de los errores por modelo
  tvar0 <- paste("^",IniObs,sep=""); name0 <- c(grep(tvar0, names(DCCP)))
  tvar1 <- paste("^",IniP,sep=""); name1 <- c(grep(tvar1, names(DCCP)))
  name0 <- name0[1:length(name1)] 
  FNC <- which(colnames(DCCP) == Fecha)
  mNC <- which(colnames(DCCP) == Model)
  
  DCCP.E <- data.frame(DCCP[,c(FNC,mNC)], DCCP[,name0] - DCCP[,name1])
  colnames(DCCP.E) <- c(colnames(DCCP.E)[1:2], paste("E",1:((ncol(DCCP.E)-2)), sep="_"))
  
  
  ## Calcula los RMSE por modelo
  Modelos <- unique(DCCP.E$model)
  RMSE <- NULL
  for(i in 1:length(Modelos)){
    #i=1
    DCCP.2 <- DCCP.E[DCCP.E$model == Modelos[i],]
    DCCP.3 <- DCCP.2[nrow(DCCP.2),]
    tvar1 <- "^E"; name1 <- c(grep(tvar1, names(DCCP.2)))
    DCCP.3[,name1] <- apply(DCCP.2[,name1], 2, RECM)
    RMSE <- rbind(RMSE, DCCP.3[-1])
  }
  
  #### Mejor modelo segun RECM
  tvar1 <- "^E"; name1 <- c(grep(tvar1, names(RMSE)))
  MinRECM <- apply(RMSE[,name1], 2, which.min)
  #RMSE[MinRECM[6],c("model",names(MinRECM)[6])]
  ##################################################################
  ###   Prueba DM
  ##################################################################
  
  tvar1 <- "^E"; name1 <- c(grep(tvar1, names(DCCP.E)))
  Modelo <- unique(DCCP.E$model)
  
  DM00 <- NULL
  for(i in 1:length(Modelo)){
    # i=55 # itera modelo
    PvalDM <- c()
    for(h in 1:(ncol(DCCP.E)-2)){
      #h=1 # itera horizonte
      #mejor <- RMSE[MinRECM[h], ] 
      mejor <- RMSE[which(RMSE$model == ModRef), ]
      xx <- DCCP.E[DCCP.E$model == mejor$model, name1[h]]
      yy <- DCCP.E[DCCP.E$model == Modelo[i], name1[h]]
      
      if(mejor$model != Modelo[i]){
        Pval0 <- try(dm.test(na.omit(xx), 
                             na.omit(yy), 
                             h=h, alternative ="two.sided",
                             varestimator = c("bartlett")),TRUE)
        if(class(Pval0)[1] == "try-error") Pval <- 1 else Pval <- Pval0$p.value
        
      }else Pval <- 1 
      PvalDM <- c(PvalDM, Pval) 
    }
    
    DM0 <- data.frame(RMSE[i,1], t(PvalDM))
    colnames(DM0) <- c(colnames(RMSE)[1], paste("DM",1:(ncol(DCCP.E)-2), sep="_"))
    DM00 <- rbind(DM00,  DM0) 
  }
  
  #### Mejor modelo segun DM
  tvar11 <- "^DM"; name11 <- c(grep(tvar11, names(DM00)))

  # modelos con RMSE menor al modelo de referencia
  RMSE2 <- RMSE
  for(j in 2:ncol(RMSE)) RMSE2[,j] <- ifelse(RMSE[,j] < RMSE[RMSE$model == ModRef,j],1,0)
  # modelos con 
  DM00a <- DM00
  for(j in 2:ncol(DM00)) DM00a[,j] <- ifelse(DM00[,j] <= 0.05,1,0)
  
  JointDM_RMSE <- RMSE2
  for(j in 2:ncol(DM00)) JointDM_RMSE[,j] <- ifelse(RMSE2[,j] + DM00a[,j] == 2,3,
                                                ifelse(RMSE2[,j]==1 & DM00a[,j] == 0, 2,1))
  
  tvar12 <- "^E"; name12 <- c(grep(tvar12, names(JointDM_RMSE)))
  MinDM <- apply(JointDM_RMSE[,name12], 2, which.max)
  
  #RMSE[MinDM[6],c("model",names(MinDM)[6])]
  
  DM <- NULL
  for(i in 1:length(Modelo)){
    # i=91 # itera modelo
    PvalDM <- c()
    for(h in 1:(ncol(DCCP.E)-2)){
      #h=1 # itera horizonte
      mejor <- RMSE[MinDM[h], ] 
      xx <- DCCP.E[DCCP.E$model == mejor$model, name1[h]]
      yy <- DCCP.E[DCCP.E$model == Modelo[i], name1[h]]
      
      if(mejor$model != Modelo[i]){
        Pval0 <- try(dm.test(na.omit(xx), 
                             na.omit(yy), 
                             h=h, alternative ="two.sided",
                             varestimator = c("bartlett")),TRUE)
        if(class(Pval0)[1] == "try-error") Pval <- 1 else Pval <- Pval0$p.value
        
      }else Pval <- 1 
      PvalDM <- c(PvalDM, Pval) 
    }
    
    DM0 <- data.frame(RMSE[i,1], t(PvalDM))
    colnames(DM0) <- c(colnames(RMSE)[1], paste("DM",1:(ncol(DCCP.E)-2), sep="_"))
    DM <- rbind(DM,  DM0) 
  }
  
  DM.Filtra <- apply(DM[,2:ncol(DM)], 2, function(x) ifelse(x > alpha, 1,0))
  colnames(DM.Filtra) <- paste("Filtra",1:(ncol(DM.Filtra)), sep="_")
  DM <- cbind(DM, DM.Filtra)
  
  tvar0 <- "^Filtra"; name10 <- c(grep(tvar0, names(DM)))
  ReportF0 <- apply(DM[,name10], 2, sum)
  ReportF <- apply(DM[,name10], 2, sum)/nrow(DM)*100
  
  ### Crea matrix con el listado de los modelos equivalentes, de acuerdo con DM test
  ###  La primera fila corresponde al porcentaje de modelos equivalentes,
  ###   la segunda fila corresponde al modelo con el min RMSE.
  ###  Cada columna representa el horizonte de pronostico
  
  BMods <- data.frame(matrix(NA,nrow = max(ReportF0), ncol=length(name10)))
  BMods[1,] <- round(ReportF,1)
  BMods[2,] <- DM[MinDM, Model]
  for(j in 1:length(name10)){
    #j=1
    equiMod <- DM[ DM[, c(1,name10[j])][,2] == 1, 1] # Filtra los modelos equivalentes
    Bmin <- which(equiMod == BMods[2,j])
    if(length(equiMod) > 1) BMods[3:(length(equiMod)+1),j] <- equiMod[-Bmin]
  }
  colnames(BMods) <- paste("h",1:length(name10), sep = "")
  
  return(list(MERROR = DCCP.E, RMSE = RMSE, DM = DM, BestMods = BMods))
  
}

#######################################################################
#  This function joints results from several models ARX, VECM, BAR
#######################################################################
AlistaDat <- function(ModT, Fecha="FechaFin",Model="model",
                      IniObs = "Obs", IniP = "h"){
  tvar0 <- paste("^",Fecha,sep=""); name0 <- c(grep(tvar0, names(ModT)))
  tvar1 <- paste("^",Model,sep=""); name1 <- c(grep(tvar1, names(ModT)))
  tvar2 <- paste("^",IniObs,sep=""); name2 <- c(grep(tvar2, names(ModT)))
  tvar3 <- paste("^",IniP,sep=""); name3 <- c(grep(tvar3, names(ModT)))
  ModT <- ModT[,c(name0,name1,name2,name3)]
  return(ModT)
}

#######################################################################
## This function calculates the frequencies of the best models
ResumMods <- function(x, Nombre1){
  # x=BestMods$h1[1:17]
  x <- as.character(na.omit(x[-1]))
  #x <- sort(x)
  p1 <- nchar(Nombre1)
  p2 <- regexpr("*[0-9]", x[1])
  AB0 <- substr(x,p1+2,p2-2)
  AB0 <- ifelse(AB0 == "",NA,AB0)
  AB <- data.frame(table(na.omit( AB0 )))
  AB$ModFreq <- paste(AB[,1],"_(",AB[,2],")",sep = "")
  
  AC0 <- substr(x,p2,1000000L)
  AC0 <- as.character(na.omit(ifelse(AC0 == "",NA,AC0)))
  
  if(nrow(AB) == 1){
    pega1 <- paste(AC0[1:(AB$Freq[1])],collapse="_")
    AB$ValAs <- pega1 
  }else if(nrow(AB) == 2){
    pega1 <- paste(AC0[1:(AB$Freq[1])],collapse="_")
    pega2 <- paste(AC0[((AB$Freq[1])+1):cumsum(AB$Freq[1:2])[2]],collapse="")
    AB$ValAs <- c(pega1,pega2) 
  }else if(nrow(AB) == 3){
    pega1 <- paste(AC0[1:(AB$Freq[1])],collapse="_")
    pega2 <- paste(AC0[((AB$Freq[1])+1):cumsum(AB$Freq[1:2])[2]],collapse="")
    pega3 <- paste(AC0[((AB$Freq[2])+1):cumsum(AB$Freq[1:3])[3]],collapse="")
    AB$ValAs <- c(pega1,pega2,pega3) 
  }
  
  AB2 <- paste(AB[,3],collapse="_")
  return(list(ModsFreqAs = AB,ModsFreq = AB2))
}


######################################################################
##   Funciones
######################################################################
Resport1 <- function(outputForest, fechaPand0 = TRUE, fecPand = "Dic_2019", 
                     antes=TRUE, i=1, filterMod=TRUE,
                     GuardaPlot=FALSE, file.out, tab=TRUE){
  
  
  #outputForest = AB2
  if(filterMod){
    outputForestE0 <- outputForest[outputForest$Mod1 == 2,]
  }else outputForestE0 <- outputForest 
  
  outputForestE0 <- outputForestE0[order(outputForestE0$Codigo),]
  outputForestE0 <- outputForestE0 %>% dplyr::filter(Codigo == 0 | Codigo > 5)
  
  
  if(fechaPand0){
    fechaPand <- min(outputForestE0[outputForestE0$fechaf == fecPand,"ordenf"])
    if(antes){
      outputForestE0 <- outputForestE0 %>% dplyr::filter(ordenf <= fechaPand)
    }else outputForestE0 <- outputForestE0 %>% dplyr::filter(ordenf > fechaPand)
    perform.Ys.Full <- CheckResiFores2(outputForestE = outputForestE0)
  }else perform.Ys.Full <- CheckResiFores2(outputForestE = outputForestE0)
  
  # Agregando errores
  # perform.Ys.Full$PerformanceAD
  # Errores desde los agregados
  #  perform.Ys.Full$PerformanceAyAfD
  ForestW <- perform.Ys.Full$ForesAyAfD
  
  tvar0 <- "^Obs"; name0 <- c(grep(tvar0, names(ForestW)))
  tvar1 <- "^h"; name1 <- c(grep(tvar1, names(ForestW)))
  #i=1
  if(GuardaPlot){
    
    png(file.out)
    
    op <- par(mfrow=c(1,2),las=3,mar=c(4,3.5,2,0)+.5, mgp=c(2.2,0.3,0), cex.axis=1.1, cex.lab=1.2)
    
    ## Filtra el pronostico usando el agregado
    ForestW.F <- ForestW[ForestW$Codigo == 0,]
    plot(ForestW.F$ordenf, ForestW.F[,name0[i]], 
         ylim=range(ForestW[,c(name0[i],name1[i])],na.rm = T),pch=20,
         ylab="Inflación", xlab="", axes=F)
    axis(2)
    axis(1, ForestW.F$ordenf, substr(unique(ForestW.F$fechaf),5,8) )
    # pronósticos del total
    lines(ForestW.F$ordenf, ForestW.F[,name1[i]], col="red")
    
    # pronósticos del agregado
    ForestW.AD <- ForestW[ForestW$Codigo == 1,]
    lines(ForestW.AD$ordenf, ForestW.AD[,name1[i]], col="blue",lty=2,lwd=2)
    
    ### Residuals
    tvar2 <- "^Resi"; name2 <- c(grep(tvar2, names(ForestW)))
    i=1
    ForestW.F <- ForestW[ForestW$Codigo == 0,]
    plot(ForestW.F$ordenf, ForestW.F[,name2[i]], 
         ylim=range(ForestW[,c(name2[i])],na.rm = T), col="red",
         ylab="Errores de pronósticos", xlab="", axes=F, type = "o")
    axis(2)
    axis(1, ForestW.F$ordenf, substr(unique(ForestW.F$fechaf),5,8) )
    
    ForestW.AD <- ForestW[ForestW$Codigo == 1,]
    lines(ForestW.AD$ordenf, ForestW.AD[,name2[i]], col="blue",lty=2,lwd=2, type = "o")
    #name2 <- c(grep(tvar2, names(ForestWBest)))
    #lines(ForestW.ADBest$ordenf, ForestW.ADBest[,name2[i]], col="green")
    
    par(op)
    dev.off()
  }
  
  Perform <- merge(perform.Ys.Full$PerformanceAD,  perform.Ys.Full$PerformanceAyAfD,
                   by.x = "h", by.y = "h")
  Perform <- Perform[,-ncol(Perform)]
  colnames(Perform) <- c("h","RMSEA.0","RMSED.0","Acc.0","PvalDM.0","SsFM.0",     
                         "RMSEA.1","RMSED.1","Acc.1","PvalDM.1")
  if(tab) Perform
  
}

# This function calculates the forecast from the disagregates for the total
#  the total should define as 0 for the code value
Resport1A <- function(outputForest){
  
  outputForestE0 <- outputForest 
  
  outputForestE0 <- outputForestE0[order(outputForestE0$Codigo),]
  outputForestE0 <- outputForestE0 %>% dplyr::filter(Codigo == 0 | Codigo > 5)
  
  perform.Ys.Full <- CheckResiFores2(outputForestE = outputForestE0)
  
  ForestW <- perform.Ys.Full$ForesAyAfD
  
  tvar0 <- "^Obs"; name0 <- c(grep(tvar0, names(ForestW)))
  tvar1 <- "^h"; name1 <- c(grep(tvar1, names(ForestW)))
  tvar2 <- "^Forest_AfD"; name2 <- c(grep(tvar2, names(ForestW)))
  #i=1
  
  Out <- ForestW[,c(which(names(ForestW) == "fechaf"),
                    which(names(ForestW) == "ordenf"),
                    name0, name1, name2)]
  return(Out)
}
