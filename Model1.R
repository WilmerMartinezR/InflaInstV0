
#############################################################################
#############################################################################
#     Este programa permite calcular pronosticos a partir de desagregados   #

#     Programa elaborado por: WILMER OSVALDO MARTINEZ RIVERA                #
#                             PhD in STATISTICS                             #

#     BANCO DE LA REPUBLICA                                                 #
#     UNIDAD DE INVESTIGACIONES                                             #

#############################################################################

#library(forecast)
#library(lawstat)
#library(stats)
#library(FinTS)
#library(TSA)
#library(glmnet)
#library(tseries)
#library(randomForest)
### Fitting the best ARIMA model

# Function parameters

# datos: data set
# fecha.ini: Initial date of the analysis, 
# fecha.fin: Last date of the analysis, 
# h: forecast horizon value

# Note: The first colum in datos must be the date

# This function fit the best arima model according to the auto.arima
# function and the AIC criterion. The object data must have at least two
#  columns where the first must be the date of the data.
# The output a data frame where each row correspoding to the results
#  per variable. The first column is the name or code of the variable,
#  the subsequent columns are:
#  initial date, end date, fitted model ARIMA(p,d,q)x(P, D, Q)[s],
#  Pvalues of ShapiroWilk, LB for the lags 5, 10 and 15, and the ARCH
#  test. Forecast for 1:12 steps ahead and a checking variable ResiFail
#  =1 if one of the residuals test fail (with alpha = 0.05), = 0 o.w. 
ModResiForest <- function(datos, fecha.ini, fecha.fin, h){
  
  model.1 <- NULL
  #y1 <- IPCVA
  FIni <- which(datos[,1] == fecha.ini)
  FEnd <- which(datos[,1] == fecha.fin)
  for(ii in 2:ncol(datos)){
    #ii <- 2
    y1 <- datos[FIni:FEnd ,ii]
    
    fit1 <- auto.arima(y1)
    orden1 <- c(fit1$arma)
    orden <- c(orden1[1],orden1[6],orden1[2])
    ordenS <- c(orden1[3],orden1[7],orden1[4], orden1[5])
    fc1 <- data.frame(forecast(fit1, h=h))
    
    ####################################################################
    ### Chequeo de residuales
    ####################################################################
    ######### #              Test Ljung and Box
    
    x <- fit1$residuals
    ###  prueba Ljung_Box para 18 rezagos
    Ljung_Box <- NULL
    for(i in 1:18){
      Ljung_Box1 <- Box.test(x, lag = i, type="Ljung")$p.value
      Ljung_Box <- rbind(Ljung_Box, cbind(i,Ljung_Box1))
    }
    colnames(Ljung_Box) <- c("Rezago","p-valor");
    
    ##############################################################
    ResidualsCheck <- c(NT_SWilk = shapiro.test(x)$p.value,
                        LBlag5 = Ljung_Box[5,2], LBlag10 = Ljung_Box[10,2],
                        LBlag15 = Ljung_Box[15,2], ARCH = ArchTest(x)$p.value) 
    
    model.1 <- rbind(model.1, c(Codigo = names(datos)[ii],
                 fechaI = fecha.ini, fechaF = fecha.fin,
                 orden, ordenS, 
                 round(ResidualsCheck,4),
                 Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                 t(fc1[,1]) ) )
  }
  model.1 <- data.frame(model.1)
  model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
  
  colnames(model.1) <- c("Codigo","fechaI","fechaf",
                         "p","d","q","P","D","Q", "S", 
                         "NT_SWilk","LBlag5","LBlag10","LBlag15",
                         "ARCH", "ResiFail", paste("h",1:h,sep="")) 
  for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
  
  return(model.1)
}

################################################
#   Run example
################################################

# datos <- data.frame(fecha=as.character(1:500), 
#                     arima.sim(list(order=c(1,0,0), ar=0.5), n=500),
#                     arima.sim(list(order=c(1,0,1), ar=0.5, ma=1), n=500))
# colnames(datos) <- c("fecha", "ind1","ind2")
# fecha.ini <- "1" 
# fecha.fin <- "100"
# h <- 12
# output <- ModResiForest(datos, fecha.ini, fecha.fin, h)


#################################################################
### Fitting the best ARIMA model taking into account the 
###    structure of the CPI index at Central Bank
#################################################################
T_ModResiForest <- function(datos, fecha.ini, fecha.fin, h){
  
  model.1 <- NULL
  #y1 <- IPCVA
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  for(ii in 1:nrow(datos)){
    #ii <- 1
    y1 <- t(datos[ii, FIni:FEnd])
    
    fit1 <- auto.arima(y1)
    orden1 <- c(fit1$arma)
    orden <- c(orden1[1],orden1[6],orden1[2])
    ordenS <- c(orden1[3],orden1[7],orden1[4], orden1[5])
    fc1 <- data.frame(forecast(fit1, h=h))
    
    ####################################################################
    ### Chequeo de residuales
    ####################################################################
    ######### #              Test Ljung and Box
    
    x <- fit1$residuals
    ###  prueba Ljung_Box para 18 rezagos
    Ljung_Box <- NULL
    for(i in 1:18){
      Ljung_Box1 <- Box.test(x, lag = i, type="Ljung")$p.value
      Ljung_Box <- rbind(Ljung_Box, cbind(i,Ljung_Box1))
    }
    colnames(Ljung_Box) <- c("Rezago","p-valor");
    
    ##############################################################
    ResidualsCheck <- c(NT_SWilk = shapiro.test(x)$p.value,
                        LBlag5 = Ljung_Box[5,2], LBlag10 = Ljung_Box[10,2],
                        LBlag15 = Ljung_Box[15,2], ARCH = ArchTest(x)$p.value) 
    
    model.1 <- rbind(model.1, c(Codigo = datos[ii,1],
                                fechaI = fecha.ini, fechaF = fecha.fin,
                                orden, ordenS, 
                                round(ResidualsCheck,4),
                                Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                t(fc1[,1]) ) )
  }
  model.1 <- data.frame(model.1)
  model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
  
  colnames(model.1) <- c("Codigo","fechaI","fechaf",
                         "p","d","q","P","D","Q", "S", 
                         "NT_SWilk","LBlag5","LBlag10","LBlag15",
                         "ARCH", "ResiFail", paste("h",1:h,sep="")) 
  for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
  
  return(model.1)
}


#################################################################
### Fitting the best ARIMA model taking into account the 
###    structure of the CPI index at Central Bank
#################################################################
T_ModResiForestARX <- function(datos, fecha.ini, fecha.fin, h, datosX){
  
  model.1 <- NULL
  #y1 <- IPCVA
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  for(ii in 1:nrow(datos)){
    #ii <- 2
    y1.0 <- t(datos[ii, FIni:FEnd])
    Xvar.0 <- t(datosX[ii, FIni:FEnd])
    
    dat0 <- as.matrix(cbind(y1.0, Xvar.0))
    colnames(dat0) <- c("ipc","ipp")
    
    dat1 <- embed(dat0,5)
    
    y1 <- dat1[,1]
    Xvar <- dat1[,c(2,4,6,8,10)] #lags effects of the PPI
    Xvar2 <- dat1[,c(3,5)] #lags effects of the CPI
    
    fit1 <- auto.arima(y1, xreg=cbind(Xvar,Xvar2^2, Xvar2^3))
    #fitX <- auto.arima(Xvar)
    #fitX <- rep(mean(Xvar[(length(Xvar)-4):length(Xvar)]),h)
    #fitX <- matrix(rep(Xvar[nrow(Xvar),],each=h), ncol = 5)
    fitX <- matrix(c(rep(Xvar[nrow(Xvar),],each=h), 
                     rep(Xvar2[nrow(Xvar2),]^2,each=h),
                     rep(Xvar2[nrow(Xvar2),]^3,each=h)),
                   ncol = 9)
    
    
    orden1 <- c(fit1$arma)
    orden <- c(orden1[1],orden1[6],orden1[2])
    ordenS <- c(orden1[3],orden1[7],orden1[4], orden1[5])
    #fcX <- data.frame(forecast(fitX, h=h))
    #colnames(fcX)[1] <- c("Xvar")
    #fc1 <- data.frame(forecast(fit1, h=h, xreg=fcX[,1]))
    fc1 <- data.frame(forecast(fit1, h=h, xreg=fitX))
    
    ####################################################################
    ### Chequeo de residuales
    ####################################################################
    ######### #              Test Ljung and Box
    
    x <- fit1$residuals
    ###  prueba Ljung_Box para 18 rezagos
    Ljung_Box <- NULL
    for(i in 1:18){
      Ljung_Box1 <- Box.test(x, lag = i, type="Ljung")$p.value
      Ljung_Box <- rbind(Ljung_Box, cbind(i,Ljung_Box1))
    }
    colnames(Ljung_Box) <- c("Rezago","p-valor");
    
    ##############################################################
    ResidualsCheck <- c(NT_SWilk = shapiro.test(x)$p.value,
                        LBlag5 = Ljung_Box[5,2], LBlag10 = Ljung_Box[10,2],
                        LBlag15 = Ljung_Box[15,2], ARCH = ArchTest(x)$p.value) 
    
    model.1 <- rbind(model.1, c(Codigo = datos[ii,1],
                                fechaI = fecha.ini, fechaF = fecha.fin,
                                orden, ordenS, 
                                round(ResidualsCheck,4),
                                Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                t(fc1[,1]) ) )
  }
  model.1 <- data.frame(model.1)
  #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
  
  colnames(model.1) <- c("Codigo","fechaI","fechaf",
                         "p","d","q","P","D","Q", "S", 
                         "NT_SWilk","LBlag5","LBlag10","LBlag15",
                         "ARCH", "ResiFail", paste("h",1:h,sep="")) 
  for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
  
  return(model.1)
}


#################################################################
### Fitting the best ARIMA model taking into account the 
###    structure of the CPI index at Central Bank
#################################################################

####################################################################
### FUNCION QUE EXTRAE ORDEN DEL MODELO auto.arima
####################################################################
OrdenArima <- function(model_fit){
  # non_seasonal_ar_order = model_fit$arma[1]
  # non_seasonal_ma_order = model_fit$arma[2]
  # 
  # seasonal_ar_order = model_fit$arma[3]
  # seasonal_ma_order = model_fit$arma[4]
  # 
  # period_of_data = model_fit$arma[5] # 1 for is non-seasonal data
  # 
  # non_seasonal_diff_order = model_fit$arma[6]
  # seasonal_diff_order =  model_fit$arma[7]
  
  p = model_fit$arma[1]
  q = model_fit$arma[2]
  
  P = model_fit$arma[3]
  Q = model_fit$arma[4]
  
  S = model_fit$arma[5] # 1 for is non-seasonal data
  
  d = model_fit$arma[6]
  D =  model_fit$arma[7]
  return(list(orden = c(p,d,q), ordenS = c(P,D,Q), Period = S))
}

####################################################################
### FUNCION Chequeo de residuales
####################################################################

Residuals <- function(x){
  ###  prueba Ljung_Box para 18 rezagos
  Ljung_Box <- NULL
  for(i in 1:18){
    Ljung_Box1 <- Box.test(x, lag = i, type="Ljung")$p.value
    Ljung_Box <- rbind(Ljung_Box, cbind(i,Ljung_Box1))
  }
  colnames(Ljung_Box) <- c("Rezago","p-valor");
  
  ##############################################################
  
  if( length(x) > 18){
    ResidualsCheck <- c(NT_SWilk = shapiro.test(x)$p.value,
                      LBlag5 = Ljung_Box[5,2], LBlag10 = Ljung_Box[10,2],
                      LBlag15 = Ljung_Box[15,2], ARCH = ArchTest(x)$p.value) 
  }else{
    ResidualsCheck <- c(NT_SWilk = NA,
                        LBlag5 = NA, LBlag10 = NA,
                        LBlag15 = NA, ARCH = NA)
  }
  return(ResidualsCheck)
}

#################################################################### 


T_ModResiForestARX_select <- function(datos, fecha.ini, fecha.fin, h, datosX){
  
  model.3 <- NULL
  #y1 <- IPCVA
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  for(ii in 1:nrow(datos)){
    #ii <- 1
    y1.0 <- t(datos[ii, FIni:FEnd])
    y1.Obs <- t(datos[ii, (FEnd+1):(FEnd+h)])
    Xvar.0 <- t(datosX[ii, FIni:FEnd])
    
    dat0 <- as.matrix(cbind(y1.0, Xvar.0))
    colnames(dat0) <- c("ipc","ipp")
    
    dat1 <- embed(dat0,6)
    
    y1 <- dat1[,1]
    
    ## Fitting the best ARIMA model
    
    fit0 <- auto.arima(y1)
    fc0 <- data.frame(forecast(fit0, h=h))
    orden0 <- OrdenArima(fit0)
    orden <- orden0$orden
    ordenS <- c(orden0$ordenS, orden0$Period)
    
    ResidualsCheck <- Residuals(fit0$residuals) 
    model.1 <- c(Codigo = datos[ii,1],
                                fechaI = fecha.ini, fechaF = fecha.fin,
                                orden, ordenS, 
                                round(ResidualsCheck,4),
                                Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                t(fc1[,1]), y1.Obs)
    
    ### Finding a better model
    Vari <- list(c(2,4), c(2,4,6),c(2,4,6,8),
                 c(2,4,6,8,10), c(2,4,6,8,10,12))
    
    for(jv in 1:length(Vari)){
      #jv=1
      v = Vari[[jv]]
      Xvar <- dat1[,v] #lags effects of the PPI
      #Xvar2 <- dat1[,c(3,5)] #lags effects of the CPI
      
  
      fit1 <- try(Arima(y1, order = orden0$orden, 
                    seasonal=list(order=orden0$ordenS, period=orden0$Period),
                    xreg=cbind(Xvar)),TRUE)
      
      fitX <- matrix(rep(Xvar[nrow(Xvar),],each=h), ncol = ncol(Xvar))
      #fitX <- matrix(c(rep(Xvar[nrow(Xvar),],each=h), 
      #                 rep(Xvar2[nrow(Xvar2),]^2,each=h),
      #                 rep(Xvar2[nrow(Xvar2),]^3,each=h)), ncol = 9)
      
    if(class(fit1)[1] == "try-error") fit1 = fit0
  
        fc1 <- data.frame(forecast(fit1, h=h, xreg=fitX))
        
        
        ResidualsCheck <- Residuals(fit1$residuals) 
        
        
        model.1 <- rbind(model.1, c(Codigo = datos[ii,1],
                                    fechaI = fecha.ini, fechaF = fecha.fin,
                                    orden, ordenS, 
                                    round(ResidualsCheck,4),
                                    Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                    t(fc1[,1]), y1.Obs) )
      }
      
      model.1 <- data.frame(model.1)
      #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
      
      colnames(model.1) <- c("Codigo","fechaI","fechaf",
                             "p","d","q","P","D","Q", "S", 
                             "NT_SWilk","LBlag5","LBlag10","LBlag15",
                             "ARCH", "ResiFail", paste("h",1:h,sep=""), paste("Mod",1:h,sep="")) 
      for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
      
      # determina cual modelo calcula el pronostico con menor error
      tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
      tvar1 <- "^Mod"; name1 <- c(grep(tvar1, names(model.1)))
      MinFores <- c()
      for(ij in 1:h) MinFores <- c(MinFores, which.min(abs(model.1[,name0[ij]] - model.1[,name1[ij]])))
      
      model.2 <- model.1[1,]
      for(ij in 1:h) model.2[1,name0[ij]] <- model.1[MinFores[ij],name0[ij]]
      for(ij in 1:h) model.2[1,name1[ij]] <- MinFores[ij]

      model.3 <- rbind(model.3, model.2)      
  }
  
  
  return(model.3)
}

#################################################################### 

# Esta función incluye IPP individual (datosX)
T_ModResiForestARX_selectTAR <- function(datos, fecha.ini, fecha.fin, h, datosX, TARMod=TRUE){
  
  model.3 <- model.Full <- NULL
  #y1 <- IPCVA
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  for(ii in 1:nrow(datos)){
    #ii <- 1
    y1.0 <- t(datos[ii, FIni:FEnd])
    y1.Obs <- t(datos[ii, (FEnd+1):(FEnd+h)])
    Xvar.0 <- t(datosX[ii, FIni:FEnd])
    
    dat0 <- as.matrix(cbind(y1.0, Xvar.0))
    colnames(dat0) <- c("ipc","ipp")
    
    dat1 <- embed(dat0,6)
    
    y1 <- dat1[,1]
    
    ## Fitting the best ARIMA model
    
    fit0 <- auto.arima(y1)
    fc0 <- data.frame(forecast(fit0, h=h))
    orden0 <- OrdenArima(fit0)
    orden <- orden0$orden
    ordenS <- c(orden0$ordenS, orden0$Period)
    
    ResidualsCheck <- Residuals(fit0$residuals) 
    model.1 <- c(Codigo = datos[ii,1],
                 fechaI = fecha.ini, fechaF = fecha.fin,
                 orden, ordenS, 
                 round(ResidualsCheck,4),
                 Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                 t(fc0[,1]), y1.Obs)
    
    ### Finding a better model
    Vari <- list(c(2,4), c(2,4,6),c(2,4,6,8),
                 c(2,4,6,8,10), c(2,4,6,8,10,12))
    
    for(jv in 1:length(Vari)){
      #jv=1
      v = Vari[[jv]]
      Xvar <- dat1[,v] #lags effects of the PPI
      #Xvar2 <- dat1[,c(3,5)] #lags effects of the CPI
      
      
      fit1 <- try(Arima(y1, order = orden0$orden, 
                        seasonal=list(order=orden0$ordenS, period=orden0$Period),
                        xreg=cbind(Xvar)),TRUE)
      
      fitX <- matrix(rep(Xvar[nrow(Xvar),],each=h), ncol = ncol(Xvar))
      #fitX <- matrix(c(rep(Xvar[nrow(Xvar),],each=h), 
      #                 rep(Xvar2[nrow(Xvar2),]^2,each=h),
      #                 rep(Xvar2[nrow(Xvar2),]^3,each=h)), ncol = 9)
      
      if(class(fit1)[1] == "try-error") fit1 = fit0
      
      fc1 <- data.frame(forecast(fit1, h=h, xreg=fitX))
      
      
      ResidualsCheck <- Residuals(fit1$residuals) 
      
      
      model.1 <- rbind(model.1, c(Codigo = datos[ii,1],
                                  fechaI = fecha.ini, fechaF = fecha.fin,
                                  orden, ordenS, 
                                  round(ResidualsCheck,4),
                                  Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                  t(fc1[,1]), y1.Obs) )
    }
    
    #########  Fitting a TAR model
    if(TARMod){
      set.seed(2357125)
      fit1 <- tar(y=y1,p1=4,p2=4,d=3,a=.1,b=.9,print=FALSE)
      fc1 <- predict(fit1, n.ahead=12,n.sim=1000)
  
      orden <- c(fit1$p1, fit1$d, fit1$p2)
      ordenS <- c(10,0,0, round(fit1$thd,4))
      ResidualsCheck <- Residuals(fit1$residuals)
          
      model.1 <- rbind(model.1, c(Codigo = datos[ii,1],
                                  fechaI = fecha.ini, fechaF = fecha.fin,
                                  orden, ordenS, 
                                  round(ResidualsCheck,4),
                                  Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                  t(fc1$fit), y1.Obs) )
    }
    
    model.1 <- data.frame(model.1)
    #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
    
    colnames(model.1) <- c("Codigo","fechaI","fechaf",
                           "p","d","q","P","D","Q", "S", 
                           "NT_SWilk","LBlag5","LBlag10","LBlag15",
                           "ARCH", "ResiFail", paste("h",1:h,sep=""), paste("Mod",1:h,sep="")) 
    for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
    
    model.1 <- rbind(model.1, model.1[1,])
    # toma el promedio de los pronosticos calculados hasta el momento
    tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
    model.1[nrow(model.1), name0] <- apply(model.1[-nrow(model.1), name0],2, mean) 
    # determina cual modelo calcula el pronostico con menor error
    tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
    tvar1 <- "^Mod"; name1 <- c(grep(tvar1, names(model.1)))
    MinFores <- c()
    for(ij in 1:h) MinFores <- c(MinFores, which.min(abs(model.1[,name0[ij]] - model.1[,name1[ij]])))
    
    model.2 <- model.1[1,]
    for(ij in 1:h) model.2[1,name0[ij]] <- model.1[MinFores[ij],name0[ij]]
    for(ij in 1:h) model.2[1,name1[ij]] <- MinFores[ij]
    
    for(ij in 1:nrow(model.1)) model.1[ij,name1] <- ij
    
    model.Full <- rbind(model.Full, model.1)
    model.3 <- rbind(model.3, model.2)      
  }
  
  
  return(list(FullModels = model.Full, BestModel = model.3))
}


#################################################################### 

# Esta función incluye IPP individual (datosX) y PC (datosX2)
T_ModResiForestARX2_selectTAR <- function(datos, fecha.ini, fecha.fin, 
                                          h, datosX,datosX2, TARMod=TRUE){
  
  model.3 <- model.Full <- NULL
  #y1 <- IPCVA
  # datosX = datosIPP; datosX2 = XT1
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  #Xvar.1 <- t(datosX2[, FIni:FEnd])
  for(ii in 1:nrow(datos)){
    #ii <- 2
    y1.0 <- t(datos[ii, FIni:FEnd])
    y1.Obs <- t(datos[ii, (FEnd+1):(FEnd+h)])
    Xvar.0 <- t(datosX[ii, FIni:FEnd])
    Xvar.1 <- t(datosX2[, (FIni-5):(FEnd-5)])
    
    dat0 <- as.matrix(cbind(y1.0, Xvar.0))
    colnames(dat0) <- c("ipc","ipp")
    
    Lags <- 6
    dat1 <- embed(dat0, Lags)
    Xvar.1 <- Xvar.1[-c(1:(Lags-1)),]
    
    y1 <- dat1[,1]
    
    ## Fitting the best ARIMA model
    
    fit0 <- auto.arima(y1)
    fc0 <- data.frame(forecast(fit0, h=h))
    orden0 <- OrdenArima(fit0)
    orden <- orden0$orden
    ordenS <- c(orden0$ordenS, orden0$Period)
    
    ResidualsCheck <- Residuals(fit0$residuals) 
    model.1 <- c(Codigo = datos[ii,1],
                 fechaI = fecha.ini, fechaF = fecha.fin,
                 orden, ordenS, 
                 round(ResidualsCheck,4),
                 Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                 t(fc0[,1]), y1.Obs)
    
    ### Finding a better model
    Vari <- list(c(2,4), c(2,4,6),c(2,4,6,8),
                 c(2,4,6,8,10), c(2,4,6,8,10,12))
    
    for(jv in 1:length(Vari)){
      #jv=1
      v = Vari[[jv]]
      Xvar <- dat1[,v] #lags effects of the PPI
      #Xvar2 <- dat1[,c(3,5)] #lags effects of the CPI
      
      covars <- cbind(Xvar, Xvar.1)
      colnames(covars) <- c(paste("Covar", 1:ncol(covars), sep=""))
      fit1 <- try(Arima(y1, order = orden0$orden, 
                        seasonal=list(order=orden0$ordenS, period=orden0$Period),
                        xreg=covars),TRUE)
      
      fitX <- matrix(rep(Xvar[nrow(Xvar),],each=h), ncol = ncol(Xvar))
      fitX2 <- matrix(rep(Xvar.1[nrow(Xvar.1),],each=h), ncol = ncol(Xvar.1))
      #fitX <- matrix(c(rep(Xvar[nrow(Xvar),],each=h), 
      #                 rep(Xvar2[nrow(Xvar2),]^2,each=h),
      #                 rep(Xvar2[nrow(Xvar2),]^3,each=h)), ncol = 9)
      
      if(class(fit1)[1] == "try-error") fit1 = fit0
      
      covarsF <- cbind(fitX, fitX2)
      colnames(covarsF) <- c(paste("Covar", 1:ncol(covars), sep=""))
      fc1 <- data.frame(forecast(fit1, h=h, xreg=covarsF))
      
      ResidualsCheck <- Residuals(fit1$residuals) 
      
      model.1 <- rbind(model.1, c(Codigo = datos[ii,1],
                                  fechaI = fecha.ini, fechaF = fecha.fin,
                                  orden, ordenS, 
                                  round(ResidualsCheck,4),
                                  Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                  t(fc1[,1]), y1.Obs) )
    }
    
    #########  Fitting a TAR model
    if(TARMod){
      set.seed(2357125)
      fit1 <- tar(y=y1,p1=4,p2=4,d=3,a=.1,b=.9,print=FALSE)
      fc1 <- predict(fit1, n.ahead=12,n.sim=1000)
      
      orden <- c(fit1$p1, fit1$d, fit1$p2)
      ordenS <- c(10,0,0, round(fit1$thd,4))
      ResidualsCheck <- Residuals(fit1$residuals)
      
      model.1 <- rbind(model.1, c(Codigo = datos[ii,1],
                                  fechaI = fecha.ini, fechaF = fecha.fin,
                                  orden, ordenS, 
                                  round(ResidualsCheck,4),
                                  Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                  t(fc1$fit), y1.Obs) )
    }
    
    model.1 <- data.frame(model.1)
    #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
    
    colnames(model.1) <- c("Codigo","fechaI","fechaf",
                           "p","d","q","P","D","Q", "S", 
                           "NT_SWilk","LBlag5","LBlag10","LBlag15",
                           "ARCH", "ResiFail", paste("h",1:h,sep=""), paste("Mod",1:h,sep="")) 
    for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
    
    model.1 <- rbind(model.1, model.1[1,])
    # toma el promedio de los pronosticos calculados hasta el momento
    tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
    model.1[nrow(model.1), name0] <- apply(model.1[-nrow(model.1), name0],2, mean) 
    # determina cual modelo calcula el pronostico con menor error
    tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
    tvar1 <- "^Mod"; name1 <- c(grep(tvar1, names(model.1)))
    MinFores <- c()
    for(ij in 1:h) MinFores <- c(MinFores, which.min(abs(model.1[,name0[ij]] - model.1[,name1[ij]])))
    
    model.2 <- model.1[1,]
    for(ij in 1:h) model.2[1,name0[ij]] <- model.1[MinFores[ij],name0[ij]]
    for(ij in 1:h) model.2[1,name1[ij]] <- MinFores[ij]
    
    for(ij in 1:nrow(model.1)) model.1[ij,name1] <- ij
    
    model.Full <- rbind(model.Full, model.1)
    model.3 <- rbind(model.3, model.2)      
  }
  
  
  return(list(FullModels = model.Full, BestModel = model.3))
}

#################################################################### 

# Esta función incluye IPP individual (datosX) y PC (datosX2)
T_ModResiForestARX3_selectTAR <- function(datos, fecha.ini, fecha.fin, 
                                          h, datosX=NULL, datosX2=NULL,
                                          datosX3=NULL, 
                                          TARMod=TRUE, IPPindv = TRUE,CompSeaso=TRUE){
  
  model.3 <- model.Full <- NULL
  #y1 <- IPCVA
  # datosX = datosIPP; datosX2 = XT1
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  #Xvar.1 <- t(datosX2[, FIni:FEnd])
  for(ii in 1:nrow(datos)){
    #ii <- 1
    y1.0 <- t(datos[ii, FIni:FEnd])
    
    # Esta condicion es para evitar el analisis con rubros nuevos, donde aun 
    #  no tienen observaciones
    #if(!all(is.na(y1.0))){
    
      y1.Obs <- t(datos[ii, (FEnd+1):(FEnd+h)])
      if(length(datosX2) > 0) Xvar.1 <- t(datosX2[, (FIni-5):(FEnd-5)]) else Xvar.1 <- NULL
      if(length(datosX3) > 0) Xvar.2 <- t(datosX3[, (FIni-5):(FEnd-5)]) else Xvar.2 <- NULL
      
      if(IPPindv){ 
        Xvar.0 <- t(datosX[ii, FIni:FEnd])
        dat0 <- as.matrix(cbind(y1.0, Xvar.0))
        colnames(dat0) <- c("ipc","ipp")
        Lags <- 6
        dat1 <- embed(dat0, Lags)
        if(length(Xvar.1) > 0) Xvar.1 <- Xvar.1[-c(1:(Lags-1)),]
        if(length(Xvar.2) > 0) Xvar.2 <- Xvar.2[-c(1:(Lags-1)),]
      }else{
        dat1 <- y1.0 
      }
      y1 <- dat1[,1]
      
      # We add this because the data for auto.arima has to be a ts
      if(CompSeaso){
        y1 <- ts(y1, frequency = 12)
      }else y1 <- y1 
      ## Fitting the best ARIMA model
      
      fit0 <- forecast::auto.arima(y1)
      fc0 <- data.frame(forecast::forecast(fit0, h=h))
      orden0 <- OrdenArima(fit0)
      orden <- orden0$orden
      ordenS <- c(orden0$ordenS, orden0$Period)
      
      ResidualsCheck <- Residuals(fit0$residuals) 
      model.1 <- c(Codigo = datos[ii,1],
                   fechaI = fecha.ini, fechaF = fecha.fin,
                   orden, ordenS, 
                   round(ResidualsCheck,4),
                   Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                   t(fc0[,1]), y1.Obs)
      
      if(IPPindv){
        ### Finding a better model INCLUDING IPP
        Vari <- list(c(2,4), c(2,4,6),c(2,4,6,8),
                     c(2,4,6,8,10), c(2,4,6,8,10,12))
        
        for(jv in 1:length(Vari)){
          #jv=1
          v = Vari[[jv]]
          Xvar <- dat1[,v] #lags effects of the PPI
          #Xvar2 <- dat1[,c(3,5)] #lags effects of the CPI
          
          covars <- cbind(Xvar, Xvar.1, Xvar.2)
          colnames(covars) <- c(paste("Covar", 1:ncol(covars), sep=""))
          fit1 <- try(forecast::Arima(y1, order = orden0$orden, 
                            seasonal=list(order=orden0$ordenS, period=orden0$Period),
                            xreg=covars, method = "ML"),TRUE)
          
          fitX <- matrix(rep(Xvar[nrow(Xvar),],each=h), ncol = ncol(Xvar))
          if(length(Xvar.1) > 0) fitX2 <- matrix(rep(Xvar.1[nrow(Xvar.1),],each=h), ncol = ncol(Xvar.1)) else fitX2 <- NULL
          if(length(Xvar.2) > 0) fitX3 <- matrix(rep(Xvar.2[nrow(Xvar.2),],each=h), ncol = ncol(Xvar.2)) else fitX3 <- NULL
          #fitX <- matrix(c(rep(Xvar[nrow(Xvar),],each=h), 
          #                 rep(Xvar2[nrow(Xvar2),]^2,each=h),
          #                 rep(Xvar2[nrow(Xvar2),]^3,each=h)), ncol = 9)
          
          if(class(fit1)[1] == "try-error") fit1 = fit0
          
          covarsF <- cbind(fitX, fitX2, fitX3)
          colnames(covarsF) <- c(paste("Covar", 1:ncol(covars), sep=""))
          fc1 <- data.frame(forecast::forecast(fit1, h=h, xreg=covarsF))
          
          ResidualsCheck <- Residuals(fit1$residuals) 
          
          model.1 <- rbind(model.1, c(Codigo = datos[ii,1],
                                      fechaI = fecha.ini, fechaF = fecha.fin,
                                      orden, ordenS, 
                                      round(ResidualsCheck,4),
                                      Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                      t(fc1[,1]), y1.Obs) )
        }
      }
      
      if(!IPPindv & length(Xvar.1) > 0 || !IPPindv & length(Xvar.2) > 0){
  
          covars <- cbind(Xvar.1, Xvar.2)
          colnames(covars) <- c(paste("Covar", 1:ncol(covars), sep=""))
          fit1 <- try(forecast::Arima(y1, order = orden0$orden, 
                            seasonal=list(order=orden0$ordenS, period=orden0$Period),
                            xreg=covars, method = "ML"),TRUE)
          
          if(class(fit1)[1] == "try-error"){
            fit1 = fit0
            fc1 <- data.frame(forecast::forecast(fit1, h=h))
          }else{
            #fitX <- matrix(rep(Xvar[nrow(Xvar),],each=h), ncol = ncol(Xvar))
            if(length(Xvar.1) > 0) fitX2 <- matrix(rep(Xvar.1[nrow(Xvar.1),],each=h ), ncol = ncol(Xvar.1)) else fitX2 <- NULL
            if(length(Xvar.2) > 0) fitX3 <- matrix(rep(Xvar.2[nrow(Xvar.2),],each=h ), ncol = ncol(Xvar.2)) else fitX3 <- NULL
            covarsF <- cbind(fitX2, fitX3)
            colnames(covarsF) <- c(paste("Covar", 1:ncol(covars), sep=""))
            fc1 <- data.frame(forecast::forecast(fit1, h=h, xreg=covarsF))  
          }  
          
          ResidualsCheck <- Residuals(fit1$residuals) 
          
          model.1 <- rbind(model.1, c(Codigo = datos[ii,1],
                                      fechaI = fecha.ini, fechaF = fecha.fin,
                                      orden, ordenS, 
                                      round(ResidualsCheck,4),
                                      Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                      t(fc1[,1]), y1.Obs) )
      }    
      
      
      #########  Fitting a TAR model
      if(TARMod){
        
        if( length( na.omit(c(y1)) ) > 20 ){
          set.seed(2357125)
          fit1 <- TSA::tar(y=y1,p1=4,p2=4,d=3,a=.1,b=.9,print=FALSE)
          fc1 <- predict(fit1, n.ahead=12,n.sim=1000)
          
          orden <- c(fit1$p1, fit1$d, fit1$p2)
          ordenS <- c(10,0,0, round(fit1$thd,4))
          ResidualsCheck <- Residuals(fit1$residuals)
          
          model.1 <- rbind(model.1, c(Codigo = datos[ii,1],
                                      fechaI = fecha.ini, fechaF = fecha.fin,
                                      orden, ordenS, 
                                      round(ResidualsCheck,4),
                                      Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                      t(fc1$fit), y1.Obs) )
        }
      }
      
      model.1 <- data.frame(model.1)
      #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
      
      colnames(model.1) <- c("Codigo","fechaI","fechaf",
                             "p","d","q","P","D","Q", "S", 
                             "NT_SWilk","LBlag5","LBlag10","LBlag15",
                             "ARCH", "ResiFail", paste("h",1:h,sep=""), paste("Mod",1:h,sep="")) 
      for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
      
      #ModMean <- FALSE
      #if(ModMean){
        model.1 <- rbind(model.1, model.1[1,])
        # toma el promedio de los pronosticos calculados hasta el momento
        tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
        model.1[nrow(model.1), name0] <- apply(model.1[-nrow(model.1), name0],2, mean) 
        # determina cual modelo calcula el pronostico con menor error
        tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
        tvar1 <- "^Mod"; name1 <- c(grep(tvar1, names(model.1)))
        MinFores <- c()
        for(ij in 1:h) MinFores <- c(MinFores, which.min(abs(model.1[,name0[ij]] - model.1[,name1[ij]])))
        
       # model.2 <- model.1[1,]
      #  for(ij in 1:h) model.2[1,name0[ij]] <- model.1[MinFores[ij],name0[ij]]
      #  for(ij in 1:h) model.2[1,name1[ij]] <- MinFores[ij]
      #}
      
      for(ij in 1:nrow(model.1)) model.1[ij,name1] <- ij
      
      model.Full <- rbind(model.Full, model.1)
      #model.3 <- rbind(model.3, model.2)      
    #}
  }
  
  #return(list(FullModels = model.Full, BestModel = model.3))
  return(list(FullModels = model.Full))
}

#################################################################### 

# Esta función estima pronósticos por el método directo
# Esta función incluye IPC (datosX) y Xs (datosX2) en la
#    parte agregada
# y Xs (datosX2) 6 rezagos de IPP desagr en los casos donde
#  se tiene par gemelo, esto para desagregados
# Todo lo anterior bajo la metodologia ARIMAX

T_ModResiForestARX4_direct <- function(datos, fecha.ini, fecha.fin, 
                                       h, datosX=NULL, datosX2=NULL,
                                       datosIPP=NULL, IPPdat=FALSE,
                                       plotReg=FALSE){
  
  # datosX <- YT1; datosX2 <- XT1
  model.3 <- model.Full <- NULL
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  
  FIniX <- which(names(datosX) == fecha.ini)
  FEndX <- which(names(datosX) == fecha.fin)
  
  FIniX2 <- which(names(datosX2) == fecha.ini)
  FEndX2 <- which(names(datosX2) == fecha.fin)
  
  FIniX3 <- which(names(datosIPP) == fecha.ini)
  FEndX3 <- which(names(datosIPP) == fecha.fin)
  
  for(ii in 1:nrow(datos)){
    #ii <- 1
    y1.0 <- t(datos[ii, FIni:FEnd])
    y1.0a <- y1.0
    y1.Obs <- t(datos[ii, (FEnd+1):(FEnd+h)])
    
    # This part consider the set of variables in a different
    # way for the aggregates and the desagregate information
    
    if(datos[ii,"Codigo"] < 10 & any(length(datosX),length(datosX2),length(IPPdat)) > 0){
      
      if(length(datosX) > 0 | length(datosX2) > 0){
        x <- data.matrix( cbind( t(datosX[, FIniX:FEndX]), 
                                 t(datosX2[, FIniX2:FEndX2]) ) )
      }else x <- NULL
      colnames(x) <- paste("V", 1:ncol(x),sep = "")
    #}else{
      if(IPPdat & datos[ii,"indPP"] == 1){
        Cod <- which(datosIPP[,"Codigo"] == datos[ii,"Codigo"])
        datSame <- embed(t(datosIPP[Cod, FIniX3:FEndX3]),7)
        complFirVals <- matrix(rep(datSame[1,],each=6),byrow=F,
                               nrow = 6,ncol = 7)
        datSame <- rbind(complFirVals, datSame)
        x <- data.matrix(cbind(datSame[,-1] , 
                               t(datosX2[, FIniX2:FEndX2]) ))
        colnames(x) <- paste("V", 1:ncol(x),sep = "")
      }else{
        x <- data.matrix( t(datosX2[, FIniX2:FEndX2]) )
        colnames(x) <- paste("V", 1:ncol(x),sep = "")
      }
    }else x <- NULL
    
    
    ## Fitting the best ARIMA model
    
    fit0 <- auto.arima(y1.0)
    #fc0 <- data.frame(forecast(fit0, h=h))
    orden0 <- OrdenArima(fit0)
    orden <- orden0$orden
    ordenS <- c(orden0$ordenS, orden0$Period)
    
    #ResidualsCheck <- Residuals(fit0$residuals) 
    #model.0 <- c(Codigo = datos[ii,1],
    #             fechaI = fecha.ini, fechaF = fecha.fin,
    #             orden, ordenS, 
    #             round(ResidualsCheck,4),
    #             Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
    #             t(fc0[,1]), y1.Obs)
    
    
    #define response variable
    model.1 <- matrix(NA, nrow=h, ncol = 40)
    
    for(hh in 1:h){
      #hh=1
      
      y <- y1.0[(hh+1):length(y1.0)]
      xtrain <- data.matrix(x[1:(nrow(x)-hh),])
      xtest <- data.matrix(x[(nrow(x)-hh+1):nrow(x),])
      
      ############################
      # Step 2: Fit the ARIMAX model
      
      #model <- glmnet(xtrain, y, alpha = KindRegre)
      fit1 <- try(Arima(y, order = orden0$orden, 
                        seasonal=list(order=orden0$ordenS, period=orden0$Period),
                        xreg=xtrain, method = "ML"),TRUE)
      
      if(class(fit1)[1] == "try-error") fit1 = fit0
      
      ###############################################
      #  Step 3: Forecast
      if(hh==1){
        fc1 <- data.frame(forecast(fit1, h=hh, xreg=t(xtest) ))
      }else fc1 <- data.frame(forecast(fit1, h=hh, xreg=xtest ))
      
      ResidualsCheck <- Residuals(fit1$residuals) 
      
      
      if(plotReg){
        y_Forest <- fc1[,1]
        if(hh == 1){
          plot(y1.Obs, ylim = range(na.omit(y1.Obs), y_Forest,0,6), type = "o", pch=20)
          lines(y_Forest, col=hh, type = "o")
        }else lines(y_Forest, col=hh, type = "o")
      }
      
      model.1[hh,1:(ncol(model.1) -12+hh) ] <- c(Codigo = datos[ii,1],
                                                 fechaI = fecha.ini, fechaF = fecha.fin,
                                                 orden, ordenS, 
                                                 round(ResidualsCheck,4),
                                                 Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                                 t(fc1[,1]), t(y1.Obs)) 
    }
    
    model.1 <- data.frame(model.1)
    #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
    
    colnames(model.1) <- c("Codigo","fechaI","fechaf",
                           "p","d","q","P","D","Q", "S", 
                           "NT_SWilk","LBlag5","LBlag10","LBlag15",
                           "ARCH", "ResiFail", paste("h",1:h,sep=""), paste("Mod",1:h,sep="")) 
    
    for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
    
    #model.1 <- rbind(model.1, model.1[1,])
    
    tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
    tvar1 <- "^Mod"; name1 <- c(grep(tvar1, names(model.1)))
    
    model.2 <- model.1[h,]
    for(ij in 1:h) model.2[1,name0[ij]] <- model.1[ij,name0[ij]]
    
    
    for(ij in 1:nrow(model.1)) model.1[ij,name1] <- ij
    
    #model.Full <- rbind(model.Full, model.1)
    model.3 <- rbind(model.3, model.2)      
    
  }
  
  #return(list(FullModels = model.Full, BestModel = model.3))
  return(list( BestModel = model.3))
}

#####################
#####################
# Esta función estima pronósticos por el método iterativo
# Esta función incluye IPC (datosX) y Xs (datosX2) en la
#    parte agregada
# y Xs (datosX2) 6 rezagos de IPP desagr en los casos donde
#  se tiene par gemelo, esto para desagregados
# Todo lo anterior bajo la metodologia ARIMAX
T_ModResiForestARX4 <- function(datos, fecha.ini, fecha.fin, 
                                h, datosX=NULL, datosX2=NULL,
                                datosIPP=NULL, IPPdat=FALSE, 
                                plotReg=FALSE){
  
  # datosX <- YT1; datosX2 <- XT1
  
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  
  FIniX <- which(names(datosX) == fecha.ini)
  FEndX <- which(names(datosX) == fecha.fin)
  
  FIniX2 <- which(names(datosX2) == fecha.ini)
  FEndX2 <- which(names(datosX2) == fecha.fin)
  
  FIniX3 <- which(names(datosIPP) == fecha.ini)
  FEndX3 <- which(names(datosIPP) == fecha.fin)
  
  model.1 <- NULL
  for(ii in 1:nrow(datos)){
    #ii <- 44
    y1.0 <- t(datos[ii, FIni:FEnd])
    y1.0a <- y1.0
    y1.Obs <- t(datos[ii, (FEnd+1):(FEnd+h)])
    
    # This part consider the set of variables in a different
    # way for the aggregates and the desagregate information
    
    if(datos[ii,"Codigo"] < 10){
      
      if(length(datosX) > 0 | length(datosX2) > 0){
        
        x <- data.matrix( cbind( t(datosX[, FIniX:FEndX]), 
                                 t(datosX2[, FIniX2:FEndX2]) ) )
      }else x <- NULL
      colnames(x) <- paste("V", 1:ncol(x),sep = "")
    }else{
      
      if(IPPdat & datos[ii,"indPP"] == 1){
        Cod <- which(datosIPP[,"Codigo"] == datos[ii,"Codigo"])
        datSame <- embed(t(datosIPP[Cod, FIniX3:FEndX3]),7)
        complFirVals <- matrix(rep(datSame[1,],each=6),byrow=F,
                               nrow = 6,ncol = 7)
        datSame <- rbind(complFirVals, datSame)
        
        x <- data.matrix(cbind(datSame[,-1] , 
                               t(datosX2[, FIniX2:FEndX2]) ))
        colnames(x) <- paste("V", 1:ncol(x),sep = "")
      }else{
        x <- data.matrix( t(datosX2[, FIniX2:FEndX2]) )
        colnames(x) <- paste("V", 1:ncol(x),sep = "")
      }
    }    
    
    ## Fitting the best ARIMA model
    
    fit0 <- auto.arima(y1.0)
    #fc0 <- data.frame(forecast(fit0, h=h))
    orden0 <- OrdenArima(fit0)
    orden <- orden0$orden
    ordenS <- c(orden0$ordenS, orden0$Period)
    
    #ResidualsCheck <- Residuals(fit0$residuals) 
    #model.0 <- c(Codigo = datos[ii,1],
    #             fechaI = fecha.ini, fechaF = fecha.fin,
    #             orden, ordenS, 
    #             round(ResidualsCheck,4),
    #             Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
    #             t(fc0[,1]), y1.Obs)
    
    
    #define response variable
    #model.1 <- matrix(NA, nrow=h, ncol = 40)
    
    #for(hh in 1:h){
    #hh=h
    
    y <- y1.0
    #xtrain <- data.matrix(x[1:(nrow(x)-hh),])
    #xtest <- data.matrix(x[(nrow(x)-hh+1):nrow(x),])
    xtrain <- data.matrix(x)
    xtest <- matrix(rep(x[nrow(x),],each=h), ncol = ncol(x))
    colnames(xtest) <- colnames(xtrain)
    ############################
    # Step 2: Fit the ARIMAX model
    
    #model <- glmnet(xtrain, y, alpha = KindRegre)
    fit1 <- try(Arima(y, order = orden0$orden, 
                      seasonal=list(order=orden0$ordenS, period=orden0$Period),
                      xreg=xtrain, method = "ML"),TRUE)
    
    if(class(fit1)[1] == "try-error"){
      fit1 = fit0
      fc1 <- data.frame(forecast(fit1, h=h))
    }else fc1 <- data.frame(forecast(fit1, h=h, xreg=xtest ))
    ###############################################
    #  Step 3: Forecast
    #if(hh==1){
    #fc1 <- data.frame(forecast(fit1, h=hh, xreg=t(xtest) ))
    #}else 
    
    
    if(plotReg){
      plot(y1.Obs, ylim = range(na.omit(y1.Obs), fc1[,1]), type = "o", pch=20)
      lines(fc1[,1], col="red", type = "o")
    }
    ResidualsCheck <- Residuals(fit1$residuals) 
    
    
    #model.1[hh,1:(ncol(model.1) -12+hh) ] <- c(Codigo = datos[ii,1],
    #                                           fechaI = fecha.ini, fechaF = fecha.fin,
    #                                           orden, ordenS, 
    #                                           round(ResidualsCheck,4),
    #                                           Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
    #                                           t(fc1[,1]), t(y1.Obs)) 
    model.1 <- rbind(model.1, c(Codigo = datos[ii,1],
                                fechaI = fecha.ini, fechaF = fecha.fin,
                                orden, ordenS, 
                                round(ResidualsCheck,4),
                                Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                t(fc1[,1]), y1.Obs) )
    
  }
  
  model.1 <- data.frame(model.1)
  #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
  
  colnames(model.1) <- c("Codigo","fechaI","fechaf",
                         "p","d","q","P","D","Q", "S", 
                         "NT_SWilk","LBlag5","LBlag10","LBlag15",
                         "ARCH", "ResiFail", paste("h",1:h,sep=""), paste("Mod",1:h,sep="")) 
  
  for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])  
  
  #return(list(FullModels = model.Full, BestModel = model.3))
  return(list( BestModel = model.1))
}

#### AR(p) directos


T_ModResiForestARp_direct <- function(datos, fecha.ini, fecha.fin, 
                                      h, plotReg=FALSE){
  
  # datosX <- YT1; datosX2 <- XT1
  model.3 <- model.Full <- NULL
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  
  for(ii in 1:nrow(datos)){
    #ii <- 1
    y1.0 <- t(datos[ii, FIni:FEnd])
    y1.0a <- y1.0
    y1.Obs <- t(datos[ii, (FEnd+1):(FEnd+h)])
    
    
    datSame <- embed(y1.0,13)
    complFirVals <- matrix(rep(datSame[1,],each=12),byrow=F,
                           nrow = 12,ncol = 13)
    datSame <- rbind(complFirVals, datSame)
    x <- data.matrix(datSame[,-1])
    colnames(x) <- paste("V", 1:ncol(x),sep = "")
    
    ## Fitting the best ARIMA model
    
    #fit0 <- auto.arima(y1.0)
    #fc0 <- data.frame(forecast(fit0, h=h))
    #orden0 <- OrdenArima(fit0)
    #orden <- orden0$orden
    #ordenS <- c(orden0$ordenS, orden0$Period)
    orden <- c(12,0,0)
    ordenS <- c(0,0,0,1)
    #ResidualsCheck <- Residuals(fit0$residuals) 
    #model.0 <- c(Codigo = datos[ii,1],
    #             fechaI = fecha.ini, fechaF = fecha.fin,
    #             orden, ordenS, 
    #             round(ResidualsCheck,4),
    #             Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
    #             t(fc0[,1]), y1.Obs)
    
    
    #define response variable
    model.1 <- matrix(NA, nrow=h, ncol = 40)
    
    for(hh in 1:h){
      #hh=1
      
      # Organizing the data to the direct forecast
      y <- y1.0[(hh+1):length(y1.0)]
      xtrain <- data.matrix(x[1:(nrow(x)-hh),])
      xtest <- data.matrix(x[(nrow(x)-hh+1):nrow(x),])
      
      ############################
      # Step 2: Fit the ARIMAX model
      
      #model <- glmnet(xtrain, y, alpha = KindRegre)
      fit1 <- forecast::Arima(y, order = c(0,0,0), 
                    xreg=xtrain, method = "ML")
      #if(class(fit1)[1] == "try-error") fit1 = fit0
      
      ###############################################
      #  Step 3: Forecast
      if(hh==1){
        fc1 <- data.frame(forecast::forecast(fit1, h=hh, xreg=t(xtest) ))
      }else fc1 <- data.frame(forecast::forecast(fit1, h=hh, xreg=xtest ))
      
      ResidualsCheck <- Residuals(fit1$residuals) 
      
      if(plotReg){
        y_Forest <- fc1[,1]
        if(hh == 1){
          plot(y1.Obs, ylim = range(na.omit(y1.Obs), y_Forest,0,6), type = "o", pch=20)
          lines(y_Forest, col=hh, type = "o")
        }else lines(y_Forest, col=hh, type = "o")
      }
      
      model.1[hh,1:(ncol(model.1) -12+hh) ] <- c(Codigo = datos[ii,1],
                                                 fechaI = fecha.ini, fechaF = fecha.fin,
                                                 orden, ordenS, 
                                                 round(ResidualsCheck,4),
                                                 Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                                 t(fc1[,1]), t(y1.Obs)) 
    }
    
    model.1 <- data.frame(model.1)
    #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
    
    colnames(model.1) <- c("Codigo","fechaI","fechaf",
                           "p","d","q","P","D","Q", "S", 
                           "NT_SWilk","LBlag5","LBlag10","LBlag15",
                           "ARCH", "ResiFail", paste("h",1:h,sep=""), paste("Mod",1:h,sep="")) 
    
    for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
    
    #model.1 <- rbind(model.1, model.1[1,])
    
    tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
    tvar1 <- "^Mod"; name1 <- c(grep(tvar1, names(model.1)))
    
    model.2 <- model.1[h,]
    for(ij in 1:h) model.2[1,name0[ij]] <- model.1[ij,name0[ij]]
    
    
    for(ij in 1:nrow(model.1)) model.1[ij,name1] <- ij
    
    #model.Full <- rbind(model.Full, model.1)
    model.3 <- rbind(model.3, model.2)      
    
  }
  
  #return(list(FullModels = model.Full, BestModel = model.3))
  return(list( BestModel = model.3))
}

#################################################################### 

# Esta función calcula pronosticos usando regresion Ridge (KindRegre = 0) 
#  y Lasso (KindRegre = 1)

T_ModResiForest_RL <- function(datos, fecha.ini, fecha.fin, 
                               h, KindRegre = 0, datosX=NULL, datosX2=NULL,
                               datosX3=NULL, 
                               IPPindv = FALSE, plotReg=FALSE){
  
  # KindRegre = 0 #Ridge; KindRegre = 1 #Lasso
  model.3 <- model.Full <- NULL
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  
  FIniX <- which(names(datosX) == fecha.ini)
  FEndX <- which(names(datosX) == fecha.fin)
  
  FIniX2 <- which(names(datosX2) == fecha.ini)
  FEndX2 <- which(names(datosX2) == fecha.fin)
  
  if(length(datosX) > 0 | length(datosX2) > 0){
    x <- data.matrix( t(datosX[, FIniX:FEndX]), t(datosX2[, FIniX2:FEndX2]) )
  }else x <- NULL
  
  for(ii in 1:nrow(datos)){
    #ii <- 1
    y1.0 <- t(datos[ii, FIni:FEnd])
    y1.0a <- y1.0
    y1.Obs <- t(datos[ii, (FEnd+1):(FEnd+h)])
  
    ## Fitting the best ARIMA model
    
    fit0 <- auto.arima(y1.0)
    #fc0 <- data.frame(forecast(fit0, h=h))
    orden0 <- OrdenArima(fit0)
    orden <- orden0$orden
    ordenS <- c(orden0$ordenS, orden0$Period)
    
    lagd <- orden[2]
    if(lagd > 0){
      y1.0 <- diff(y1.0, differences = lagd)
      y1.0 <- c(y1.0[1:lagd],y1.0)
    }
    
    # lagAR <- orden[1]+1
    # if(lage >= 2){
    #   yAR <- embed(y1.0,lagAR)
    #   yAR <- rbind( matrix(rep(y1.0[1:(lagAR-1)], each=ncol(yAR)), 
    #                         nrow= lagAR-1, ncol=ncol(yAR), byrow = T), yAR)
    #   yAR <- data.frame(yAR)
    # }
    
    # at <- rnorm(length(y1.0))
    # lagMA <- orden[3]+1
    # if(lage >= 2){
    #   yMA <- embed(at,lagMA)
    #   yMA <- data.frame(yMA)
    # }
    
    ############################
    # Step 1: Load the Data
    
    # To perform ridge regression, we’ll use functions from the glmnet package. 
    # This package requires the response variable to be a vector and the set 
    #  of predictor variables to be of the class data.matrix.
    
    #define response variable
    model.1 <- matrix(NA, nrow=h, ncol = 40)
    
    for(hh in 1:h){
      #hh=1
      
      y <- y1.0[(hh+1):length(y1.0)]
      xtrain <- data.matrix(x[1:(nrow(x)-hh),])
      xtest <- x[(nrow(x)-hh+1):nrow(x),]
      
      ############################
      # Step 2: Fit the Ridge Regression Model
      
      # Next, we’ll use the glmnet() function to fit the ridge regression model 
      #   and specify alpha=0.
      # Note that setting alpha equal to 1 is equivalent to using Lasso Regression 
      #   and setting alpha to some value between 0 and 1 is equivalent to using an elastic net.
      
      ## Also note that ridge regression requires the data to be standardized 
      #   such that each predictor variable has a mean of 0 and a standard deviation of 1.
      
      
      #fit ridge regression model
      model <- glmnet(xtrain, y, alpha = KindRegre)
      
      #view summary of model
      # summary(model)
      # plot(model, label = TRUE)
      #print(model)
      ###############################################
      #  Step 3: Choose an Optimal Value for Lambda
      
      #perform k-fold cross-validation to find optimal lambda value
      cv_model <- cv.glmnet(xtrain, y, alpha = 0)
      
      #find optimal lambda value that minimizes test MSE
      best_lambda <- cv_model$lambda.min
      # best_lambda
      
      #produce plot of test MSE by lambda value
      # plot(cv_model) 
      
      ###############################################
      #    Step 4: Analyze Final Model
      
      #find coefficients of best model
      
      best_model <- glmnet(xtrain, y, alpha = 0, lambda = best_lambda)
      # coef(best_model)
      
      #produce Ridge trace plot
      # plot(model, xvar = "lambda")
      
      y_predicted <- predict(model, s = best_lambda, newx = xtest)
      #y_Forest <- Integra(y1.Obs, y_predicted, y1.0a[length(y1.0a)])
      y_Forest <- Integra(y_predicted, y1.0a[length(y1.0a)])
      
      
      if(plotReg){
        if(hh == 1){
          plot(y1.Obs, ylim = range(na.omit(y1.Obs), y_Forest,0,6), type = "o", pch=20)
          lines(y_Forest, col=hh, type = "o")
        }else lines(y_Forest, col=hh, type = "o")
      }
      
      y_predictedF <- predict(model, s = best_lambda, newx = xtrain)
      Resi <- y-y_predictedF
      
      
      ResidualsCheck <- Residuals(Resi) 
      model.1[hh,1:(ncol(model.1) -12+hh) ] <- c(Codigo = datos[ii,1],
                                                 fechaI = fecha.ini, fechaF = fecha.fin,
                                                 orden, ordenS, 
                                                 round(ResidualsCheck,4),
                                                 Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                                 t(y_Forest), t(y1.Obs)) 
    }
    
    model.1 <- data.frame(model.1)
    #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
    
    colnames(model.1) <- c("Codigo","fechaI","fechaf",
                           "p","d","q","P","D","Q", "S", 
                           "NT_SWilk","LBlag5","LBlag10","LBlag15",
                           "ARCH", "ResiFail", paste("h",1:h,sep=""), paste("Mod",1:h,sep="")) 
    
    for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
    
    #model.1 <- rbind(model.1, model.1[1,])
    
    tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
    tvar1 <- "^Mod"; name1 <- c(grep(tvar1, names(model.1)))
    
    model.2 <- model.1[h,]
    for(ij in 1:h) model.2[1,name0[ij]] <- model.1[ij,name0[ij]]
    
    
    for(ij in 1:nrow(model.1)) model.1[ij,name1] <- ij
    
    #model.Full <- rbind(model.Full, model.1)
    model.3 <- rbind(model.3, model.2)      
    
  }
  
  #return(list(FullModels = model.Full, BestModel = model.3))
  return(list( BestModel = model.3))
}

### Version 2
# Esta función calcula pronosticos usando regresion Ridge (KindRegre = 0) 
#  y Lasso (KindRegre = 1)

T_ModResiForest_RL2 <- function(datos, fecha.ini, fecha.fin, 
                               h, KindRegre = 0, datosX=NULL, datosX2=NULL,
                               datosX3=NULL, datosCla=NULL, 
                               IPPindv = FALSE, plotReg=FALSE){
  
  # KindRegre = 0 #Ridge; KindRegre = 1 #Lasso
  model.3 <- model.Full <- NULL
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  
  FIniX <- which(names(datosX) == fecha.ini)
  FEndX <- which(names(datosX) == fecha.fin)
  
  FIniX2 <- which(names(datosX2) == fecha.ini)
  FEndX2 <- which(names(datosX2) == fecha.fin)
  
  for(ii in 1:nrow(datos)){
    #ii <- 17
    y1.0 <- t(datos[ii, FIni:FEnd])
    y1.0a <- y1.0
    y1.Obs <- t(datos[ii, (FEnd+1):(FEnd+h)])
    
    # This part consider the set of variables in a different
    # way for the aggregates and the desagregate information
    
    # This condition is for the aggregates
    if(datos[ii,"Codigo"] < 10){
      
      if(length(datosX) > 0 | length(datosX2) > 0){
        
        # This part classified the Xs by services, goods and others
        filtraSG <- datos[ii,c("Codigo","Tipo")] 
        if(filtraSG$Tipo == "Bienes"){
          filtraX <- datosCla %>% dplyr::filter(BIENES == 1)
          datosX2a <- datosX2 %>%
            dplyr::filter(rownames(datosX2) %in% filtraX$CodVariable)
        }else if(filtraSG$Tipo == "Servicios"){
          filtraX <- datosCla %>% dplyr::filter(SERVICIOS == 1)
          datosX2a <- datosX2 %>%
            dplyr::filter(rownames(datosX2) %in% filtraX$CodVariable)
        }else datosX2a <- datosX2
        
        x <- data.matrix( cbind( t(datosX[, FIniX:FEndX]), 
                                 t(datosX2a[, FIniX2:FEndX2]) ) )
      }else x <- NULL
      colnames(x) <- paste("V", 1:ncol(x),sep = "")
    
    # This condition is for the disagregates
    }else{
      if(length(datosX) > 0 | length(datosX2) > 0){
        
        Cod <- which(datosX[,"Codigo"] == datos[ii,"Codigo"])
        datSame <- embed(t(datosX[Cod, FIniX:FEndX]),13)
        complFirVals <- matrix(rep(datSame[1,],each=12),byrow=F,
                               nrow = 12,ncol = 13)
        datSame <- rbind(complFirVals, datSame)
        
        #Here is the filter for the Xs
      
        # This part classified the Xs by services, goods and others
        filtraSG <- datos[ii,c("Codigo","Tipo")] 
        if(filtraSG$Tipo == "b"){
          filtraX <- datosCla %>% dplyr::filter(BIENES == 1)
          datosX2a <- datosX2 %>%
            dplyr::filter(rownames(datosX2) %in% filtraX$CodVariable)
        }else if(filtraSG$Tipo == "s"){
          filtraX <- datosCla %>% dplyr::filter(SERVICIOS == 1)
          datosX2a <- datosX2 %>%
            dplyr::filter(rownames(datosX2) %in% filtraX$CodVariable)
        }else datosX2a <- datosX2
        
        x <- data.matrix(cbind(datSame[,-1] , 
                               t(datosX2a[, FIniX2:FEndX2]) ))
      }else x <- NULL
    
      colnames(x) <- paste("V", 1:ncol(x),sep = "")
    }    
    
    ## Fitting the best ARIMA model
    
    fit0 <- auto.arima(y1.0)
    #fc0 <- data.frame(forecast(fit0, h=h))
    orden0 <- OrdenArima(fit0)
    orden <- orden0$orden
    ordenS <- c(orden0$ordenS, orden0$Period)
    
    lagd <- orden[2]
    if(lagd > 0){
      y1.0 <- diff(y1.0, differences = lagd)
      y1.0 <- c(y1.0[1:lagd],y1.0)
    }
    
    # lagAR <- orden[1]+1
    # if(lage >= 2){
    #   yAR <- embed(y1.0,lagAR)
    #   yAR <- rbind( matrix(rep(y1.0[1:(lagAR-1)], each=ncol(yAR)), 
    #                         nrow= lagAR-1, ncol=ncol(yAR), byrow = T), yAR)
    #   yAR <- data.frame(yAR)
    # }
    
    # at <- rnorm(length(y1.0))
    # lagMA <- orden[3]+1
    # if(lage >= 2){
    #   yMA <- embed(at,lagMA)
    #   yMA <- data.frame(yMA)
    # }
    
    ############################
    # Step 1: Load the Data
    
    # To perform ridge regression, we’ll use functions from the glmnet package. 
    # This package requires the response variable to be a vector and the set 
    #  of predictor variables to be of the class data.matrix.
    
    #define response variable
    model.1 <- matrix(NA, nrow=h, ncol = 40)
    
    for(hh in 1:h){
      #hh=1
      
      y <- y1.0[(hh+1):length(y1.0)]
      xtrain <- data.matrix(x[1:(nrow(x)-hh),])
      xtest <- x[(nrow(x)-hh+1):nrow(x),]
      
      ############################
      # Step 2: Fit the Ridge Regression Model
      
      # Next, we’ll use the glmnet() function to fit the ridge regression model 
      #   and specify alpha=0.
      # Note that setting alpha equal to 1 is equivalent to using Lasso Regression 
      #   and setting alpha to some value between 0 and 1 is equivalent to using an elastic net.
      
      ## Also note that ridge regression requires the data to be standardized 
      #   such that each predictor variable has a mean of 0 and a standard deviation of 1.
      
      
      #fit ridge regression model
      model <- glmnet(xtrain, y, alpha = KindRegre)
      
      #view summary of model
      # summary(model)
      # plot(model, label = TRUE)
      #print(model)
      ###############################################
      #  Step 3: Choose an Optimal Value for Lambda
      
      #perform k-fold cross-validation to find optimal lambda value
      cv_model <- cv.glmnet(xtrain, y, alpha = 0)
      
      #find optimal lambda value that minimizes test MSE
      best_lambda <- cv_model$lambda.min
      # best_lambda
      
      #produce plot of test MSE by lambda value
      # plot(cv_model) 
      
      ###############################################
      #    Step 4: Analyze Final Model
      
      #find coefficients of best model
      
      best_model <- glmnet(xtrain, y, alpha = 0, lambda = best_lambda)
      # coef(best_model)
      
      #produce Ridge trace plot
      # plot(model, xvar = "lambda")
      
      y_predicted <- predict(model, s = best_lambda, newx = xtest)
      
      if(lagd > 0){
        y_Forest <- Integra(y_predicted, y1.0a[length(y1.0a)])
      }else y_Forest <- y_predicted
      
      if(plotReg){
        if(hh == 1){
          plot(y1.Obs, ylim = range(na.omit(y1.Obs), y_Forest,0,6), type = "o", pch=20)
          lines(y_Forest, col=hh, type = "o")
        }else lines(y_Forest, col=hh, type = "o")
      }
      
      y_predictedF <- predict(model, s = best_lambda, newx = xtrain)
      Resi <- y-y_predictedF
      
      
      ResidualsCheck <- Residuals(Resi) 
      model.1[hh,1:(ncol(model.1) -12+hh) ] <- c(Codigo = datos[ii,1],
                                                 fechaI = fecha.ini, fechaF = fecha.fin,
                                                 orden, ordenS, 
                                                 round(ResidualsCheck,4),
                                                 Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                                 t(y_Forest), t(y1.Obs)) 
    }
    
    model.1 <- data.frame(model.1)
    #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
    
    colnames(model.1) <- c("Codigo","fechaI","fechaf",
                           "p","d","q","P","D","Q", "S", 
                           "NT_SWilk","LBlag5","LBlag10","LBlag15",
                           "ARCH", "ResiFail", paste("h",1:h,sep=""), paste("Mod",1:h,sep="")) 
    
    for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
    
    #model.1 <- rbind(model.1, model.1[1,])
    
    tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
    tvar1 <- "^Mod"; name1 <- c(grep(tvar1, names(model.1)))
    
    model.2 <- model.1[h,]
    for(ij in 1:h) model.2[1,name0[ij]] <- model.1[ij,name0[ij]]
    
    
    for(ij in 1:nrow(model.1)) model.1[ij,name1] <- ij
    
    #model.Full <- rbind(model.Full, model.1)
    model.3 <- rbind(model.3, model.2)      
    
  }
  
  #return(list(FullModels = model.Full, BestModel = model.3))
  return(list( BestModel = model.3))
}

#################################################################### 

# Esta función calcula pronosticos usando random forest (KindRegre = 0) 
#  y Lasso (KindRegre = 1)

T_ModResiForest_RF <- function(datos, fecha.ini, fecha.fin, 
                               h, KindRegre = 0, datosX=NULL, datosX2=NULL,
                               datosX3=NULL, 
                               IPPindv = FALSE, plotReg=FALSE){
  
  # KindRegre = 0 #Ridge; KindRegre = 1 #Lasso
  model.3 <- model.Full <- NULL
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  
  FIniX <- which(names(datosX) == fecha.ini)
  FEndX <- which(names(datosX) == fecha.fin)
  
  FIniX2 <- which(names(datosX2) == fecha.ini)
  FEndX2 <- which(names(datosX2) == fecha.fin)
  
  if(length(datosX) > 0 | length(datosX2) > 0){
    x <- data.matrix( t(datosX[, FIniX:FEndX]), t(datosX2[, FIniX2:FEndX2]) )
  }else x <- NULL
  colnames(x) <- paste("V", colnames(x),sep = "")
  
  for(ii in 1:nrow(datos)){
    #ii <- 1
    y1.0 <- t(datos[ii, FIni:FEnd])
    y1.0a <- y1.0
    y1.Obs <- t(datos[ii, (FEnd+1):(FEnd+h)])
    
    # x <- data.matrix( t(datosX[, FIniX:FEndX]) )
    # colnames(x) <- paste("V", colnames(x),sep = "")
    ## Fitting the best ARIMA model
    
    fit0 <- auto.arima(y1.0)
    #fc0 <- data.frame(forecast(fit0, h=h))
    orden0 <- OrdenArima(fit0)
    orden <- orden0$orden
    ordenS <- c(orden0$ordenS, orden0$Period)
    
    lagd <- orden[2]
    if(lagd > 0){
      y1.0 <- diff(y1.0, differences = lagd)
      y1.0 <- c(y1.0[1:lagd],y1.0)
    }
    

    #define response variable
    model.1 <- matrix(NA, nrow=h, ncol = 40)
    
    for(hh in 1:h){
      #hh=1
      
      y <- y1.0[(hh+1):length(y1.0)]
      xtrain <- data.matrix(x[1:(nrow(x)-hh),])
      xtest <-  data.matrix(x[(nrow(x)-hh+1):nrow(x),])
      
      ############################
      # Step 2: 
      #fit the Random Forest regression model
      
      #model <- glmnet(xtrain, y, alpha = KindRegre)
      model <- randomForest( y ~ ., data = xtrain,
                               importance = TRUE, mtry=hh,
                               na.action = na.omit, ntree = 500)
      
      ############################
      # Step 3: Prediction 
      if(hh==1) xtest <- t(xtest)
      y_predicted <- predict(model, xtest)
      #y_Forest <- Integra(y1.Obs, y_predicted, y1.0a[length(y1.0a)])
      y_Forest <- Integra(y_predicted, y1.0a[length(y1.0a)])
      
      if(plotReg){
        if(hh == 1){
          plot(y1.Obs, ylim = range(na.omit(y1.Obs), y_Forest,0,6), type = "o", pch=20)
          lines(y_Forest, col=hh, type = "o")
        }else lines(y_Forest, col=hh, type = "o")
      }
      
      y_predictedF <- predict(model, xtrain)
      Resi <- y-y_predictedF
      
      
      ResidualsCheck <- Residuals(Resi) 
      model.1[hh,1:(ncol(model.1) -12+hh) ] <- c(Codigo = datos[ii,1],
                                                 fechaI = fecha.ini, fechaF = fecha.fin,
                                                 orden, ordenS, 
                                                 round(ResidualsCheck,4),
                                                 Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                                 t(y_Forest), t(y1.Obs)) 
    }
    
    model.1 <- data.frame(model.1)
    #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
    
    colnames(model.1) <- c("Codigo","fechaI","fechaf",
                           "p","d","q","P","D","Q", "S", 
                           "NT_SWilk","LBlag5","LBlag10","LBlag15",
                           "ARCH", "ResiFail", paste("h",1:h,sep=""), paste("Mod",1:h,sep="")) 
    
    for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
    
    #model.1 <- rbind(model.1, model.1[1,])
    
    tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
    tvar1 <- "^Mod"; name1 <- c(grep(tvar1, names(model.1)))
    
    model.2 <- model.1[h,]
    for(ij in 1:h) model.2[1,name0[ij]] <- model.1[ij,name0[ij]]
    
    
    for(ij in 1:nrow(model.1)) model.1[ij,name1] <- ij
    
    #model.Full <- rbind(model.Full, model.1)
    model.3 <- rbind(model.3, model.2)      
    
  }
  
  #return(list(FullModels = model.Full, BestModel = model.3))
  return(list( BestModel = model.3))
}

#################################################################### 

# Esta función calcula pronosticos usando random forest (KindRegre = 0) 
#  y Lasso (KindRegre = 1)

T_ModResiForest_RF2 <- function(datos, fecha.ini, fecha.fin, 
                               h, KindRegre = 0, datosX=NULL, datosX2=NULL,
                               datosX3=NULL, datosCla=NULL,
                               IPPindv = FALSE, plotReg=FALSE){
  
  # KindRegre = 0 #Ridge; KindRegre = 1 #Lasso
  model.3 <- model.Full <- NULL
  FIni <- which(names(datos) == fecha.ini)
  FEnd <- which(names(datos) == fecha.fin)
  
  FIniX <- which(names(datosX) == fecha.ini)
  FEndX <- which(names(datosX) == fecha.fin)
  
  FIniX2 <- which(names(datosX2) == fecha.ini)
  FEndX2 <- which(names(datosX2) == fecha.fin)
  
  for(ii in 1:nrow(datos)){
    #ii <- 1
    y1.0 <- t(datos[ii, FIni:FEnd])
    y1.0a <- y1.0
    y1.Obs <- t(datos[ii, (FEnd+1):(FEnd+h)])
    
    # This part consider the set of variables in a different
    # way for the aggregates and the desagregate information
    
    # This condition is for the aggregates
    if(datos[ii,"Codigo"] < 10){
      
      if(length(datosX) > 0 | length(datosX2) > 0){
        
        # This part classified the Xs by services, goods and others
        filtraSG <- datos[ii,c("Codigo","Tipo")] 
        if(filtraSG$Tipo == "Bienes"){
          filtraX <- datosCla %>% dplyr::filter(BIENES == 1)
          datosX2a <- datosX2 %>%
            dplyr::filter(rownames(datosX2) %in% filtraX$CodVariable)
        }else if(filtraSG$Tipo == "Servicios"){
          filtraX <- datosCla %>% dplyr::filter(SERVICIOS == 1)
          datosX2a <- datosX2 %>%
            dplyr::filter(rownames(datosX2) %in% filtraX$CodVariable)
        }else datosX2a <- datosX2
        
        x <- data.matrix( cbind( t(datosX[, FIniX:FEndX]), 
                                 t(datosX2a[, FIniX2:FEndX2]) ) )
      }else x <- NULL
      colnames(x) <- paste("V", 1:ncol(x),sep = "")
      
      # This condition is for the disagregates
    }else{
      if(length(datosX) > 0 | length(datosX2) > 0){
        
        Cod <- which(datosX[,"Codigo"] == datos[ii,"Codigo"])
        datSame <- embed(t(datosX[Cod, FIniX:FEndX]),13)
        complFirVals <- matrix(rep(datSame[1,],each=12),byrow=F,
                               nrow = 12,ncol = 13)
        datSame <- rbind(complFirVals, datSame)
        
        #Here is the filter for the Xs
        
        # This part classified the Xs by services, goods and others
        filtraSG <- datos[ii,c("Codigo","Tipo")] 
        if(filtraSG$Tipo == "b"){
          filtraX <- datosCla %>% dplyr::filter(BIENES == 1)
          datosX2a <- datosX2 %>%
            dplyr::filter(rownames(datosX2) %in% filtraX$CodVariable)
        }else if(filtraSG$Tipo == "s"){
          filtraX <- datosCla %>% dplyr::filter(SERVICIOS == 1)
          datosX2a <- datosX2 %>%
            dplyr::filter(rownames(datosX2) %in% filtraX$CodVariable)
        }else datosX2a <- datosX2
        
        x <- data.matrix(cbind(datSame[,-1] , 
                               t(datosX2a[, FIniX2:FEndX2]) ))
      }else x <- NULL
      
      colnames(x) <- paste("V", 1:ncol(x),sep = "")
    }    
    
    ## Fitting the best ARIMA model
    
    fit0 <- forecast::auto.arima(y1.0)
    #fc0 <- data.frame(forecast(fit0, h=h))
    orden0 <- OrdenArima(fit0)
    orden <- orden0$orden
    ordenS <- c(orden0$ordenS, orden0$Period)
    
    lagd <- orden[2]
    if(lagd > 0){
      y1.0 <- diff(y1.0, differences = lagd)
      y1.0 <- c(y1.0[1:lagd],y1.0)
    }
    
    
    #define response variable
    model.1 <- matrix(NA, nrow=h, ncol = 40)
    
    for(hh in 1:h){
      #hh=1
      
      y <- y1.0[(hh+1):length(y1.0)]
      xtrain <- data.matrix(x[1:(nrow(x)-hh),])
      xtest <-  data.matrix(x[(nrow(x)-hh+1):nrow(x),])
      
      ############################
      # Step 2: 
      #fit the Random Forest regression model
      
      #model <- glmnet(xtrain, y, alpha = KindRegre)
      model <- randomForest( y ~ ., data = xtrain,
                             importance = TRUE, 
                             na.action = na.omit, ntree = 500)
      
      ############################
      # Step 3: Prediction 
      if(hh==1) xtest <- t(xtest)
      y_predicted <- predict(model, xtest)
      #y_Forest <- Integra(y1.Obs, y_predicted, y1.0a[length(y1.0a)])
      if(lagd > 0){
        y_Forest <- Integra(y_predicted, y1.0a[length(y1.0a)])
      }else y_Forest <- y_predicted
        
      if(plotReg){
        if(hh == 1){
          plot(y1.Obs, ylim = range(na.omit(y1.Obs), y_Forest,0,14), type = "o", pch=20)
          lines(y_Forest, col=hh, type = "o")
        }else lines(y_Forest, col=hh, type = "o")
      }
      
      y_predictedF <- predict(model, xtrain)
      Resi <- y-y_predictedF
      
      
      ResidualsCheck <- Residuals(Resi) 
      model.1[hh,1:(ncol(model.1) -12+hh) ] <- c(Codigo = datos[ii,1],
                                                 fechaI = fecha.ini, fechaF = fecha.fin,
                                                 orden, ordenS, 
                                                 round(ResidualsCheck,4),
                                                 Fail = ifelse(any(ResidualsCheck < 0.05),1,0),
                                                 t(y_Forest), t(y1.Obs)) 
    }
    
    model.1 <- data.frame(model.1)
    #model.1$Fail <- ifelse(any(model.1[,11:15] < 0.05),1,0)
    
    colnames(model.1) <- c("Codigo","fechaI","fechaf",
                           "p","d","q","P","D","Q", "S", 
                           "NT_SWilk","LBlag5","LBlag10","LBlag15",
                           "ARCH", "ResiFail", paste("h",1:h,sep=""), paste("Mod",1:h,sep="")) 
    
    for(w in 4:ncol(model.1)) model.1[,w] <- as.numeric(model.1[,w])
    
    #model.1 <- rbind(model.1, model.1[1,])
    
    tvar0 <- "^h"; name0 <- c(grep(tvar0, names(model.1)))
    tvar1 <- "^Mod"; name1 <- c(grep(tvar1, names(model.1)))
    
    model.2 <- model.1[h,]
    for(ij in 1:h) model.2[1,name0[ij]] <- model.1[ij,name0[ij]]
    
    
    for(ij in 1:nrow(model.1)) model.1[ij,name1] <- ij
    
    #model.Full <- rbind(model.Full, model.1)
    model.3 <- rbind(model.3, model.2)      
    
  }
  
  #return(list(FullModels = model.Full, BestModel = model.3))
  return(list( BestModel = model.3))
}
#################################################################
###  Checking the forecast residuals
#################################################################

CheckResiFores <- function(outputForestE, h=12){

  tvar0 <- "^h"; name0 <- c(grep(tvar0, names(outputForestE)))
  tvar1 <- "^Obs"; name1 <- c(grep(tvar1, names(outputForestE)))
  
  ResiF <- outputForestE[,c(1:4, name0, name1, name1)]
  
  tvar2 <- "Obs1.1"; name2 <- c(grep(tvar2, names(ResiF)))
  ResiF[,name2:(name2+h-1)] <- outputForestE[,name1] - outputForestE[,name0]
  colnames(ResiF) <- c(colnames(ResiF)[1:(name2-1)], paste("ResiA",1:h,sep=""))
  
  # Filtering only the results from the Agregate Total = 0
  ResiA <- ResiF %>% dplyr::filter(Codigo == 0)
  
  # Filtering only the results from the desagregate and agregate the resuduals
  # from the disagregates
  
  #Pesos <- datos$Pesos2018[-1]/100
  ResiI <- ResiF %>% dplyr::filter(Codigo != 0) %>%
    dplyr::group_by(fechaf) %>% 
    summarise(across(starts_with("ResiA"), list(sum = ~sum(.x*Pesos2018,na.rm=TRUE))))
  # The previous line aplies over al columns that start with the name ResiA
  colnames(ResiI) <- c("fechaf", paste("ResiI",1:h,sep=""))
  
  # Joint the Residuals from aggregate forecast and the ones aggregating from
  #  dessagregates
  ResiAll <- merge(ResiA, ResiI, by.x="fechaf", by.y="fechaf")
  
  ## summary RMSE calculation
  tvar0 <- "^ResiA"; name0 <- c(grep(tvar0, names(ResiAll)))
  tvar1 <- "^ResiI"; name1 <- c(grep(tvar1, names(ResiAll)))
  
  RMSE <- function(x) sqrt(mean(x^2,na.rm = TRUE))
  RMSE.output <- apply(ResiAll[,c(name0, name1)],2, RMSE)
  
  # Diebold and Mariano test
  # For alternative="greater", the alternative hypothesis is that method 2 
  #   is more accurate than method 1. ECM1 > ECM2
  PvalDM <- c()
  for(fw in 1:h) PvalDM <- c(PvalDM, dm.test(na.omit(ResiAll[,c(name0[fw])]), 
                                             na.omit(ResiAll[,c(name1[fw])]), 
                                             h=fw, alternative ="greater")$p.value)
  
  Performance <- data.frame(h = 1:h, RMSEA = RMSE.output[1:h], 
                            RMSED = RMSE.output[(h+1):(2*h)])
  Performance$AccurateD <- ifelse(Performance$RMSEA > Performance$RMSED, 1,0)
  Performance$PvalDM <- PvalDM
  
  # Calcula tamaño fuera de muestra para cada horizonte
  tvar0 <- "^ResiA"; name0 <- c(grep(tvar0, names(ResiA)))
  SsFM <- c()
  for(kk in 1:length(name0)) SsFM <- c(SsFM, length(na.omit(ResiA[,name0[kk]])) )
  Performance$SsFM <- SsFM
  
  return(Performance)
}

###############################################################
#  This function calculates the forecast from the dissgregates
#  Keep the obs, the agre forecast (h) and the forecast from the dissgregates
#  Also, it includes the performance

CheckResiFores2 <- function(outputForestE, h=12){
  
  tvar0 <- "^h"; name0 <- c(grep(tvar0, names(outputForestE)))
  tvar1 <- "^Obs"; name1 <- c(grep(tvar1, names(outputForestE)))
  
  ####################################################
  
     #  Note1: El codigo del agregado debe ser 0
     #  Estas líneas agregan pronósticos por fecha 
  Forest_AfD0 <- outputForestE %>% dplyr::filter(Codigo != 0) %>%
    dplyr::group_by(fechaf) %>% 
    dplyr::summarise(across(starts_with("h"), list(sum = ~sum(.x*(Pesos2018/sum(Pesos2018)),na.rm=TRUE))))
  # The previous line aplies over al columns that start with the name ResiA
  colnames(Forest_AfD0) <- c("fechaf", paste("Forest_AfD",1:h,sep=""))
  
    # Filtra el total
  Forest_AfD1 <- outputForestE %>% dplyr::filter(Codigo == 0)
  
    # Une total y agregado desde desagregados (AD)
  Forest_AfD1 <- merge(Forest_AfD1, Forest_AfD0, by.x="fechaf", by.y="fechaf")
    # Organiza por fecha
  Forest_AfD1 <- Forest_AfD1[order(Forest_AfD1$ordenf),]
  OutForest <- Forest_AfD1
    # Define posiciones de variables cuyo nombre comienza por h y Forest...
  tvar0AD <- "^h"; name0AD <- c(grep(tvar0AD, names(Forest_AfD1)))
  tvar1AD <- "^Forest_AfD"; name1AD <- c(grep(tvar1AD, names(Forest_AfD1)))
  
    # Redefine Codigo == 1, asi los pronosticos agregados tiene la misma estructure
     ##    del tottal Codigo =0
  Forest_AfD1$Codigo <- 1 
    # Reemplaza los valores de h con los de los pronosticos AfD (aggre from Desa)
  Forest_AfD1[,name0AD] <- Forest_AfD1[,name1AD]
    # Elimina columnas "^Forest_AfD"
  Forest_AfD1 <- Forest_AfD1[,-name1AD]
  ####################################################
    # Une el agregado a los datos totales
  outputForestE <- rbind(Forest_AfD1, outputForestE)
  
    # Mantiene las siguientes columnas y duplica name1 para almacenar residuos
  ResiF <- outputForestE[,c(1:4, name0, name1, name1)]
    
  tvar2 <- "Obs1.1"; name2 <- c(grep(tvar2, names(ResiF)))
    # Calcula residuos
  ResiF[,name2:(name2+h-1)] <- outputForestE[,name1] - outputForestE[,name0]
  colnames(ResiF) <- c(colnames(ResiF)[1:(name2-1)], paste("ResiA",1:h,sep=""))
  
  # Filtering only the results from the Total = 0
  ResiA <- ResiF %>% dplyr::filter(Codigo <= 1)
  
  # Filtering only the results from the desagregate, and agregate the residuals
  # from the disagregates
  
  # Agrega residuos
  #Pesos <- datos$Pesos2018[-1]/100
  ResiI <- ResiF %>% dplyr::filter(Codigo > 1) %>%
    dplyr::group_by(fechaf) %>% 
    dplyr::summarise(across(starts_with("ResiA"), list(sum = ~sum(.x*(Pesos2018/sum(Pesos2018)),na.rm=TRUE))))
  # The previous line aplies over al columns that start with the name ResiA
  colnames(ResiI) <- c("fechaf", paste("ResiI",1:h,sep=""))
  
  # Joint the Residuals from aggregate forecast and the ones aggregating from
  #  dessagregates
  
  ResiAll <- merge(ResiA[ResiA$Codigo == 0, ], ResiI, by.x="fechaf", by.y="fechaf")
  ## summary RMSE calculation
  tvar0 <- "^ResiA"; name0 <- c(grep(tvar0, names(ResiAll)))
  tvar1 <- "^ResiI"; name1 <- c(grep(tvar1, names(ResiAll)))
  
  RMSE <- function(x) sqrt(mean(x^2,na.rm = TRUE))
  RMSE.output <- apply(ResiAll[,c(name0, name1)],2, RMSE)
  
  # Diebold and Mariano test
  # For alternative="greater", the alternative hypothesis is that method 2 
  #   is more accurate than method 1. ECM1 > ECM2
  PvalDM <- c()
  for(fw in 1:h) PvalDM <- c(PvalDM, dm.test(na.omit(ResiAll[,c(name0[fw])]), 
                                             na.omit(ResiAll[,c(name1[fw])]), 
                                             h=fw, alternative ="two.sided",
                                             varestimator = c("bartlett"))$p.value)
  
  Performance <- data.frame(h = 1:h, RMSEA = RMSE.output[1:h], 
                            RMSED = RMSE.output[(h+1):(2*h)])
  Performance$AccurateD <- ifelse(Performance$RMSEA > Performance$RMSED, 1,0)
  Performance$PvalDM <- PvalDM
  
  # Calcula tamaño fuera de muestra para cada horizonte
  tvar0 <- "^ResiA"; name0 <- c(grep(tvar0, names(ResiA)))
  SsFM <- c()
  for(kk in 1:length(name0)) SsFM <- c(SsFM, 
                                length(na.omit(ResiA[ResiA$Codigo == 0,name0[kk]])) )
  Performance$SsFM <- SsFM
  Performance0 <- Performance
  
  ###
  ###########
  tvar0 <- "^ResiA"; name0 <- c(grep(tvar0, names(ResiA)))
      # Filtra pronosticos agregados y residuos desde lo agregado
  ResiI0 <-  ResiA[ResiA$Codigo == 1, c(1,name0)]
  colnames(ResiI0) <- c("fechaf", paste("ResiI",1:h,sep=""))
  ResiAll <- merge(ResiA[ResiA$Codigo == 0, ], ResiI0, by.x="fechaf", by.y="fechaf")
  ## summary RMSE calculation
  tvar0 <- "^ResiA"; name0 <- c(grep(tvar0, names(ResiAll)))
  tvar1 <- "^ResiI"; name1 <- c(grep(tvar1, names(ResiAll)))
  
  RMSE <- function(x) sqrt(mean(x^2,na.rm = TRUE))
  RMSE.output <- apply(ResiAll[,c(name0, name1)],2, RMSE)
  
  # Diebold and Mariano test
  # For alternative="greater", the alternative hypothesis is that method 2 
  #   is more accurate than method 1. ECM1 > ECM2
  PvalDM <- c()
  for(fw in 1:h) PvalDM <- c(PvalDM, dm.test(na.omit(ResiAll[,c(name0[fw])]), 
                                             na.omit(ResiAll[,c(name1[fw])]), 
                                             h=fw, alternative ="two.sided",
                                             varestimator = c("bartlett"))$p.value)
  
  Performance <- data.frame(h = 1:h, RMSEA = RMSE.output[1:h], 
                            RMSED = RMSE.output[(h+1):(2*h)])
  Performance$AccurateD <- ifelse(Performance$RMSEA > Performance$RMSED, 1,0)
  Performance$PvalDM <- PvalDM
  
  # Calcula tamaño fuera de muestra para cada horizonte
  tvar0 <- "^ResiA"; name0 <- c(grep(tvar0, names(ResiA)))
  SsFM <- c()
  for(kk in 1:length(name0)) SsFM <- c(SsFM, 
                            length(na.omit(ResiA[ResiA$Codigo == 1,name0[kk]])) )
  Performance$SsFM <- SsFM
  
  Performance1 <- Performance
  
  # Nota: Performance0 contiene evaluación de pronósticos agregando errores
  #       Performance1 contiene evaluación de pronósticos agregando pronosticos
  #                    y calculando errores desde los agregados
  return(list(PerformanceAD = Performance0, PerformanceAyAfD = Performance1,
              ForesAyAfD = OutForest))
}

# This function calculates the individual RMSE
# outputForestE <- outputForestE.IPP 
# h =12
CheckResiFores.All <- function(outputForestE, h=12){
  
  tvar0 <- "^h"; name0 <- c(grep(tvar0, names(outputForestE)))
  tvar1 <- "^Obs"; name1 <- c(grep(tvar1, names(outputForestE)))
  
  # Keeps id variables and forecast, obs x2
  posicod <- which(colnames(outputForestE) == "Codigo")
  ResiF <- outputForestE[,c(posicod, name0, name1, name1)]
  
  tvar2 <- "Obs1.1"; name2 <- c(grep(tvar2, names(ResiF)))
  # calculate residuals
  ResiF[,name2:(name2+h-1)] <- outputForestE[,name1] - outputForestE[,name0]
  colnames(ResiF) <- c(colnames(ResiF)[1:(name2-1)], paste("ResiA",1:h,sep=""))
  
  RMSE <- function(x) sqrt(mean(x^2,na.rm = TRUE))
  
  RMSE.output <- ResiF %>% 
    dplyr::group_by(Codigo) %>% 
    summarise(across(starts_with("ResiA"), list(RMSE = ~RMSE(.x))))
  
  return(RMSE.output)
}

CheckResiFores.All.standa <- function(outputForestE, h=12){
  
  tvar0 <- "^h"; name0 <- c(grep(tvar0, names(outputForestE)))
  tvar1 <- "^Obs"; name1 <- c(grep(tvar1, names(outputForestE)))
  
  # Keeps id variables and forecast, obs x2
  ResiF <- outputForestE[,c(1:4, name0, name1, name1)]
  
  tvar2 <- "Obs1.1"; name2 <- c(grep(tvar2, names(ResiF)))
  # calculate residuals
  ResiF[,name2:(name2+h-1)] <- outputForestE[,name1] - outputForestE[,name0]
  colnames(ResiF) <- c(colnames(ResiF)[1:(name2-1)], paste("ResiA",1:h,sep=""))
  
  standa <- function(x) return( (x-mean(x,na.rm=T))/sd(x,na.rm=T) )
  tvar2 <- "^ResiA"; name2 <- c(grep(tvar2, names(ResiF)))
  for(ij in 1:length(name2)) ResiF[,name2[ij]] <- standa(ResiF[,name2[ij]])
  
  RMSE <- function(x) sqrt(mean(x^2,na.rm = TRUE))
  
  RMSE.output <- ResiF %>% 
    dplyr::group_by(Codigo) %>% 
    summarise(across(starts_with("ResiA"), list(RMSE = ~RMSE(.x))))
  
  return(RMSE.output)
}
