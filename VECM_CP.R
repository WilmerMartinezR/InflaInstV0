

#alpha <- 0.05
############################################################################
###   Functions
############################################################################

Check_UR <- function(serie, Alpha = 0.01, year0=1, month0=1, frequ=12){
  serie.1 <- ts(serie, start = c(year0, month0), frequency = frequ)
  suppressWarnings({ 
    # Code that generates warning messages 
    UnitRootsKpss <- data.frame(rbind(cbind(Pval=tseries::kpss.test(serie.1)$p.value,Diff=0),
                                      cbind(tseries::kpss.test(diff(serie.1,differences = 1))$p.value,Diff=1),
                                      cbind(tseries::kpss.test(diff(serie.1,differences = 2))$p.value,Diff=2)))
  })
  
  suppressWarnings({ 
    # Code that generates warning messages 
    UnitRootsADF <- data.frame(rbind(cbind(Pval=tseries::adf.test(serie.1)$p.value,Diff=0),
                                     cbind(tseries::adf.test(diff(serie.1,differences = 1))$p.value,Diff=1),
                                     cbind(tseries::adf.test(diff(serie.1,differences = 2))$p.value,Diff=2)))
  })
  
  UnitRootsKpss$EstacKpss <- ifelse(UnitRootsKpss$Pval >= Alpha,1,0)
  UnitRootsADF$EstacADF <- ifelse(UnitRootsADF$Pval <= Alpha,1,0)
  
  UnitRoots <- data.frame(UnitRootsKpss[,-1],EstacADF=UnitRootsADF$EstacADF)
  UnitRoots$Rule <- ifelse(UnitRoots$EstacKpss == 1 & 
                             UnitRoots$EstacADF == 1,1, 
                           ifelse(UnitRoots$EstacKpss == 1 & 
                                    UnitRoots$EstacADF == 0,2,0))
  ## This part keep the minimum lag
  PosiD <- anyDuplicated(UnitRoots$Rule)
  if(PosiD > 0) UnitRoots$Rule[PosiD] <- 0
  PosiD2 <- anyDuplicated(UnitRoots$Rule)
  if(PosiD2 > 0) UnitRoots$Rule[PosiD2] <- 0
  
  return(UnitRoots)
}  

# Test unit root kpss and Differenciation
KP.Pval <- function(X1,X2){
  #Xs=x1 
  kp1 <- Check_UR(X1)
  kp2 <- Check_UR(X2)
  
  d1 <- which.max(kp1$Rule)[1]-1
  d2 <- which.max(kp2$Rule)[1]-1
  
  X1.d <- X1
  X2.d <- X2
  # if(d1 > 0) X1.d <- diff(X1, diff = d1) else X1.d <- X1
  # if(d2 > 0) X2.d <- diff(X2, diff = d2) else X2.d <- X2
  # 
  # if(d1 > d2) X1.d <- X1.d[-c(1:(d1-d2))]
  # if(d1 < d2) X2.d <- X2.d[-c(1:(d2-d1))]
  
  return(list(Dat = cbind(X1.d, X2.d), Dif=c(d1,d2)))
}




# Fit an VAR(Lags) model beteen a answer (x1) and a predictor (x2) and its
#  lags according to the the coincident profile analysis.
#  Lags <- abs(LagA)+1

# The output is a lits with two matrices (parameters and forecast) 

VECMCP.fit <- function(FechaIni, FechaFin, Nombre1, x1,x2,LagA,a, hf, 
                       GB_Obs, alpha=0.05){
  
  # FechaIni: Fecha inicial de la muestra 
  # FechaFin: Fecha final de la muestra 
  # Nombre1: Nombre de la variable de an?lisis 
  # x1: datos de la variable respuesta, vector
  # x2: datos de la variable respuesta, vector
  # LagA: Rezago encontrado segun perfil coincidente
  # a: Valor de la transformacion inflainst
  # hf: Valor del horizonte de pronostico
  
  
  #### Ajuste VAR
  Lags <- abs(LagA)+1
  
  #XY <- KP.Pval(x1,x2)
  XY <- cbind(x1,x2)

  ecdet <- "const"
  exogen <- NULL
  alpha_coi='5pct'
  pvar <- Lags
  test_coi<-test_cointeg(j = 1, series= XY, CI = pvar, ecdet = ecdet, 
                         exogen = exogen, alpha_coint = alpha_coi, 
                         summary = FALSE)[[1]] # Los objetos salida de la funcion test_cointeg se guardan en call
  # orden coint
  r=test_coi$r
  
  if(r >= 1){
  	if(r == dim(XY)[2]) r = r - 1
  
    # Para el modelo vec2_Var
    test_type=test_coi$test_type
    sjf.vecm1 <- urca::ca.jo(XY, type=test_type, ecdet=ecdet, K=pvar+1, 
                             spec="transitory",dumvar=exogen) # Asigna a sjf.vecm el modelo VEC estimado con el procedimiento de Johansen
    vecm1<-VECM(XY,r=r,lag=pvar,estim = "ML",include="none",
                LRinclude=ecdet,exogen=exogen) # Asigna a vecm1 el modelo VEC estimado con el procedimiento de Johansen
    #summary(vecm1)
    #coefB(vecm1)
    #coefA(vecm1)
    #coefPI(vecm1)
    
    vecm<-cajorls(sjf.vecm1,r=r) # Asigna a vecm el ajuste del modelo VEC. vecm1
            # y vecm son equivalentes. De vecm1 se extraer?n m?sdelante los 
            # coeficientes de cointegraci?n y de vecm se extraer?n los 
            # coeficientes de la variable ex?gena si es el caso
    ParamV <- as.matrix.data.frame(coeftest(vecm$rlm))[1:(2*pvar+1),c(1,4)]
    ParamV <- data.frame(ParamV)
    colnames(ParamV) <- c("coef","Pval")
    rownames(ParamV) <- c("const",
                          paste(rep(c("x1","x2"),pvar), 
                                rep(1:pvar,each=2),sep="_")  )
    
    ParamV$Sig <- ifelse(ParamV$Pval < alpha,1, 0)
    
    ParamB <- data.frame(round(c(vecm1$model.specific$beta[1,1],
                  -vecm1$model.specific$beta[-1,1]),6), 
                  c(rep(NA,3)), c(rep(NA,3)))
    ParamA <- data.frame( round(vecm1$coefficients[,1:r],6), 
                         c(rep(NA,2)), c(rep(NA,2)))
    colnames(ParamB) <- c("coef","Pval","Sig")
    rownames(ParamB) <- paste(c("x1","x2","const"),"B",sep="")
    colnames(ParamA) <- c("coef","Pval","Sig")
    rownames(ParamA) <- paste(c("x1","x2"),"A",sep="")
  
    Params <- rbind(ParamV, ParamB, ParamA)
    Params$Npar <- 1:nrow(Params)
    
    resi_modVEC<-vecm$rlm$residuals[,1] #residuales de la variable dependiente
    ResidualsCheck <- Residuals(resi_modVEC)
    
    ########################################################################
    #### Pronos
    ########################################################################
    fFin <- which(substr(GB_Obs$Fecha,1,7) == FechaFin)
    y1.Obs <- t(GB_Obs[(fFin+1):(fFin+hf),Nombre1])
    colnames(y1.Obs) <- paste("Obs",1:hf,sep = "")
    
    vec2_var_mod<-vec2var(sjf.vecm1,r=r)  #representacion VAR del modelo VEC
    fc01 = predict(vec2_var_mod, n.ahead = hf)[[1]][[1]][,1]
    
    ForesC <- data.frame(t(fc01))
    colnames(ForesC) <- paste("h",1:hf,sep = "")
    
    Out0 <- data.frame(FechaIni = FechaIni,  
                       FechaFin = FechaFin, Agregado = Nombre1, 
                       MaxLT = LagA,
                       MaxL = Lags-1,
                       a=a, 
                       t(ResidualsCheck), y1.Obs, ForesC)
    
    Output <- list(ParametrosE = Params, Pronosticos = Out0)
  }else{
    Output <- list(ParametrosE = NA, Pronosticos = NA)
  }
  return(Output)
}

##########################################################################

# This function eval if a VAR or VECM are appropriate. The order of the model
# is taken according to CP result

VECMCP.fit2 <- function(FechaIni, FechaFin, Nombre1, x1,x2,LagA,a, hf, 
                       GB_Obs, alpha=0.05){
  
  # FechaIni: Fecha inicial de la muestra 
  # FechaFin: Fecha final de la muestra 
  # Nombre1: Nombre de la variable de an?lisis 
  # x1: datos de la variable respuesta, vector
  # x2: datos de la variable respuesta, vector
  # LagA: Rezago encontrado segun perfil coincidente
  # a: Valor de la transformacion inflainst
  # hf: Valor del horizonte de pronostico
  
  XY <- cbind(x1,x2)

  #### Ajuste VAR
  Lags <- abs(LagA)+1
  ecdet <- "const"
  exogen = newdata_exo = NULL
  #x1 <- X1; x2 <- X2
  OrdenP <- vars::VARselect(XY)
  pvar <- Lags
  EvalUR <- KP.Pval(x1,x2)
  
  # If there are unit roots
  if(EvalUR$Dif[1] > 0 || EvalUR$Dif[2] > 0){
    
    alpha_coi='5pct'
    test_coi<-test_cointeg(j = 1, series= XY, CI = pvar, ecdet = ecdet, 
                           exogen = exogen, alpha_coint = alpha_coi, 
                           summary = FALSE)[[1]] # Los objetos salida de la funcion test_cointeg se guardan en call
    # orden coint
    r=test_coi$r
    if(r >= 1){
      if(r == dim(XY)[2]) r = r - 1
      
      # Para el modelo vec2_Var
      test_type=test_coi$test_type
      sjf.vecm1 <- urca::ca.jo(XY, type=test_type, ecdet=ecdet, K=pvar+1, 
                               spec="transitory",dumvar=exogen) # Asigna a sjf.vecm el modelo VEC estimado con el procedimiento de Johansen
      vecm1<-VECM(XY,r=r,lag=pvar,estim = "ML",include="none",
                  LRinclude=ecdet,exogen=exogen) # Asigna a vecm1 el modelo VEC estimado con el procedimiento de Johansen
      #summary(vecm1)
      #coefB(vecm1)
      #coefA(vecm1)
      #coefPI(vecm1)
      
      vecm<-cajorls(sjf.vecm1,r=r) # Asigna a vecm el ajuste del modelo VEC. vecm1
      # y vecm son equivalentes. De vecm1 se extraer?n m?sdelante los 
      # coeficientes de cointegraci?n y de vecm se extraer?n los 
      # coeficientes de la variable ex?gena si es el caso
      ParamV <- as.matrix.data.frame(coeftest(vecm$rlm))[1:(2*pvar+1),c(1,4)]
      ParamV <- data.frame(ParamV)
      colnames(ParamV) <- c("coef","Pval")
      rownames(ParamV) <- c("const",
                            paste(rep(c("x1","x2"),pvar), 
                                  rep(1:pvar,each=2),sep="_")  )
      
      ParamV$Sig <- ifelse(ParamV$Pval < alpha,1, 0)
      
      ParamB <- data.frame(round(c(vecm1$model.specific$beta[1,1],
                                   -vecm1$model.specific$beta[-1,1]),6), 
                           c(rep(NA,3)), c(rep(NA,3)))
      ParamA <- data.frame( round(vecm1$coefficients[,1:r],6), 
                            c(rep(NA,2)), c(rep(NA,2)))
      colnames(ParamB) <- c("coef","Pval","Sig")
      rownames(ParamB) <- paste(c("x1","x2","const"),"B",sep="")
      colnames(ParamA) <- c("coef","Pval","Sig")
      rownames(ParamA) <- paste(c("x1","x2"),"A",sep="")
      
      Params <- rbind(ParamV, ParamB, ParamA)
      Params$Npar <- 1:nrow(Params)
      
      resi_modVEC<-vecm$rlm$residuals[,1] #residuales de la variable dependiente
      ResidualsCheck <- Residuals(resi_modVEC)
      
      ########################################################################
      #### Pronos
      ########################################################################
      fFin <- which(substr(GB_Obs$Fecha,1,7) == FechaFin)
      y1.Obs <- t(GB_Obs[(fFin+1):(fFin+hf),Nombre1])
      colnames(y1.Obs) <- paste("Obs",1:hf,sep = "")
      
      vec2_var_mod<-vec2var(sjf.vecm1,r=r)  #representacion VAR del modelo VEC
      fc01 = predict(vec2_var_mod, n.ahead = hf)[[1]][[1]][,1]
      
      ForesC <- data.frame(t(fc01))
      colnames(ForesC) <- paste("h",1:hf,sep = "")
      
      Out0 <- data.frame(FechaIni = FechaIni,  
                         FechaFin = FechaFin, Agregado = Nombre1, 
                         MaxLT = LagA,
                         MaxL = Lags-1,
                         a=a, 
                         t(ResidualsCheck), y1.Obs, ForesC)
      
      UR = EvalUR$Dif
      names(UR) <- c("d1","d2")
      Output <- list(ParametrosE = Params, Pronosticos = Out0, 
                     Orden = cbind(t(OrdenP$selection-1), UR = t(UR), Coint_r = r))
      #return(Output)
    }else{
      UR = EvalUR$Dif
      names(UR) <- c("d1","d2")
      Output <- list(ParametrosE = NA, Pronosticos = NA, 
                     Orden = cbind(t(OrdenP$selection-1), UR = t(UR), Coint_r = r ))
    }
  }else if(EvalUR$Dif[1] == 0 && EvalUR$Dif[2] == 0){
    
    modelo_VAR=vars::VAR(XY,p=pvar,type=ecdet) # Estima el modelo VAR y guarda la lista de salida en modelo_VAR. 
    
    ParamV <- as.matrix.data.frame(coeftest(modelo_VAR))[1:(2*pvar+1),c(1,4)]
    ParamV <- data.frame(ParamV)
    colnames(ParamV) <- c("coef","Pval")
    rownames(ParamV) <- c("const",
                          paste(rep(c("x1","x2"),pvar), 
                                rep(1:pvar,each=2),sep="_")  )
    
    ParamV$Sig <- ifelse(ParamV$Pval < alpha,1, 0)
    Params <- ParamV
    Params$Npar <- 1:nrow(Params)

    res_modVAR=modelo_VAR$varresult[[1]]$residuals # Guarda los residuales de la variable Dependiente
    ResidualsCheck <- Residuals(res_modVAR)
    
    ########################################################################
    #### Pronos
    ########################################################################
    fFin <- which(substr(GB_Obs$Fecha,1,7) == FechaFin)
    y1.Obs <- t(GB_Obs[(fFin+1):(fFin+hf),Nombre1])
    colnames(y1.Obs) <- paste("Obs",1:hf,sep = "")
    
    fc01 = predict(modelo_VAR, n.ahead = hf)[[1]][[1]][,1]
    
    ForesC <- data.frame(t(fc01))
    colnames(ForesC) <- paste("h",1:hf,sep = "")
    
    Out0 <- data.frame(FechaIni = FechaIni,  
                       FechaFin = FechaFin, Agregado = Nombre1, 
                       MaxLT = LagA,
                       MaxL = Lags-1,
                       a=a, 
                       t(ResidualsCheck), y1.Obs, ForesC)
    
    UR = EvalUR$Dif
    names(UR) <- c("d1","d2")
    Output <- list(ParametrosE = Params, Pronosticos = Out0, 
                   Orden = cbind(t(OrdenP$selection-1), UR = t(UR), Coint_r = 0))
    
  }
  return(Output)
}

##########################################################################

# This function eval if a VAR or VECM are appropriate. The order of the model
# is taken according to the varSelect function

VECMCP.fit3 <- function(FechaIni, FechaFin, Nombre1, x1,x2,LagA,a, hf, 
                        GB_Obs, alpha=0.05, fFin){
  
  # FechaIni: Fecha inicial de la muestra 
  # FechaFin: Fecha final de la muestra 
  # Nombre1: Nombre de la variable de an?lisis 
  # x1: datos de la variable respuesta, vector
  # x2: datos de la variable respuesta, vector
  # LagA: Rezago encontrado segun perfil coincidente
  # a: Valor de la transformacion inflainst
  # hf: Valor del horizonte de pronostico
  
  XY <- cbind(x1,x2)
  
  #### Ajuste VAR
  Lags <- abs(LagA)+1
  ecdet <- "const"
  exogen = newdata_exo = NULL
  #x1 <- X1; x2 <- X2
  OrdenP <- vars::VARselect(XY)
  #pvar <- Lags
  EvalUR <- KP.Pval(x1,x2)
  
  # If there are unit roots
  if(EvalUR$Dif[1] > 0 || EvalUR$Dif[2] > 0){
    
    pvar <- as.numeric(OrdenP$selection[1])-1
    alpha_coi='5pct'
    test_coi<-test_cointeg(j = 1, series= XY, CI = pvar, ecdet = ecdet, 
                           exogen = exogen, alpha_coint = alpha_coi, 
                           summary = FALSE)[[1]] # Los objetos salida de la funcion test_cointeg se guardan en call
    # orden coint
    r=test_coi$r
    if(r >= 1){
      if(r == dim(XY)[2]) r = r - 1
      
      # Para el modelo vec2_Var
      test_type=test_coi$test_type
      sjf.vecm1 <- urca::ca.jo(XY, type=test_type, ecdet=ecdet, K=pvar+1, 
                               spec="transitory",dumvar=exogen) # Asigna a sjf.vecm el modelo VEC estimado con el procedimiento de Johansen
      vecm1<-VECM(XY,r=r,lag=pvar,estim = "ML",include="none",
                  LRinclude=ecdet,exogen=exogen) # Asigna a vecm1 el modelo VEC estimado con el procedimiento de Johansen
      #summary(vecm1)
      #coefB(vecm1)
      #coefA(vecm1)
      #coefPI(vecm1)
      
      vecm<-cajorls(sjf.vecm1,r=r) # Asigna a vecm el ajuste del modelo VEC. vecm1
      # y vecm son equivalentes. De vecm1 se extraer?n m?sdelante los 
      # coeficientes de cointegraci?n y de vecm se extraer?n los 
      # coeficientes de la variable ex?gena si es el caso
      ParamV <- as.matrix.data.frame(coeftest(vecm$rlm))[1:(2*pvar+1),c(1,4)]
      ParamV <- data.frame(ParamV)
      colnames(ParamV) <- c("coef","Pval")
      rownames(ParamV) <- c("const",
                            paste(rep(c("x1","x2"),pvar), 
                                  rep(1:pvar,each=2),sep="_")  )
      
      ParamV$Sig <- ifelse(ParamV$Pval < alpha,1, 0)
      
      ParamB <- data.frame(round(c(vecm1$model.specific$beta[1,1],
                                   -vecm1$model.specific$beta[-1,1]),6), 
                           c(rep(NA,3)), c(rep(NA,3)))
      ParamA <- data.frame( round(vecm1$coefficients[,1:r],6), 
                            c(rep(NA,2)), c(rep(NA,2)))
      colnames(ParamB) <- c("coef","Pval","Sig")
      rownames(ParamB) <- paste(c("x1","x2","const"),"B",sep="")
      colnames(ParamA) <- c("coef","Pval","Sig")
      rownames(ParamA) <- paste(c("x1","x2"),"A",sep="")
      
      Params <- rbind(ParamV, ParamB, ParamA)
      Params$Npar <- 1:nrow(Params)
      
      resi_modVEC<-vecm$rlm$residuals[,1] #residuales de la variable dependiente
      ResidualsCheck <- Residuals(resi_modVEC)
      
      ########################################################################
      #### Pronos
      ########################################################################
      fFin <- which(substr(GB_Obs$Fecha,1,7) == FechaFin)
      y1.Obs <- t(GB_Obs[(fFin+1):(fFin+hf),Nombre1])
      colnames(y1.Obs) <- paste("Obs",1:hf,sep = "")
      
      vec2_var_mod<-vec2var(sjf.vecm1,r=r)  #representacion VAR del modelo VEC
      fc01 = predict(vec2_var_mod, n.ahead = hf)[[1]][[1]][,1]
      
      ForesC <- data.frame(t(fc01))
      colnames(ForesC) <- paste("h",1:hf,sep = "")
      
      Out0 <- data.frame(FechaIni = FechaIni,  
                         FechaFin = FechaFin, Agregado = Nombre1, 
                         MaxLT = LagA,
                         MaxL = Lags-1,
                         a=a, 
                         t(ResidualsCheck), y1.Obs, ForesC)
      
      UR = EvalUR$Dif
      names(UR) <- c("d1","d2")
      Output <- list(ParametrosE = Params, Pronosticos = Out0, 
                     Orden = cbind(t(OrdenP$selection-1), UR = t(UR), Coint_r = r))
      #return(Output)
    }else{
      UR = EvalUR$Dif
      names(UR) <- c("d1","d2")
      Output <- list(ParametrosE = NA, Pronosticos = NA, 
                     Orden = cbind(t(OrdenP$selection-1), UR = t(UR), Coint_r = r ))
    }
  }else if(EvalUR$Dif[1] == 0 && EvalUR$Dif[2] == 0){
    
    pvar <- as.numeric(OrdenP$selection[1])
    modelo_VAR=vars::VAR(XY,p=pvar,type=ecdet) # Estima el modelo VAR y guarda la lista de salida en modelo_VAR. 
    
    ParamV <- as.matrix.data.frame(coeftest(modelo_VAR))[1:(2*pvar+1),c(1,4)]
    ParamV <- data.frame(ParamV)
    colnames(ParamV) <- c("coef","Pval")
    rownames(ParamV) <- c("const",
                          paste(rep(c("x1","x2"),pvar), 
                                rep(1:pvar,each=2),sep="_")  )
    
    ParamV$Sig <- ifelse(ParamV$Pval < alpha,1, 0)
    Params <- ParamV
    Params$Npar <- 1:nrow(Params)
    
    res_modVAR=modelo_VAR$varresult[[1]]$residuals # Guarda los residuales de la variable Dependiente
    ResidualsCheck <- Residuals(res_modVAR)
    
    ########################################################################
    #### Pronos
    ########################################################################
    fFin <- which(substr(GB_Obs$Fecha,1,7) == FechaFin)
    y1.Obs <- t(GB_Obs[(fFin+1):(fFin+hf),Nombre1])
    colnames(y1.Obs) <- paste("Obs",1:hf,sep = "")
    
    fc01 = predict(modelo_VAR, n.ahead = hf)[[1]][[1]][,1]
    
    ForesC <- data.frame(t(fc01))
    colnames(ForesC) <- paste("h",1:hf,sep = "")
    
    Out0 <- data.frame(FechaIni = FechaIni,  
                       FechaFin = FechaFin, Agregado = Nombre1, 
                       MaxLT = LagA,
                       MaxL = Lags-1,
                       a=a, 
                       t(ResidualsCheck), y1.Obs, ForesC)
    
    UR = EvalUR$Dif
    names(UR) <- c("d1","d2")
    Output <- list(ParametrosE = Params, Pronosticos = Out0, 
                   Orden = cbind(t(OrdenP$selection-1), UR = t(UR), Coint_r = 0))
    
  }
  return(Output)
}

VARMod <- FALSE
if(VARMod){
  OP1 <- MTS::VARorder(XY)
  OP2 <- vars::VARselect(XY)
  
  M1 <- MTS::VAR(XY, p=Lags, output = F)
  M2 <- vars::VAR(XY, p=Lags)
  #M2 <- vars:: (XY, p=Lags)
  
  OCoint <- urca::ca.jo(XY,K=Lags)
  summary(OCoint)
  
  coefs<- data.frame(summary(M2)[[2]][[1]][['coefficients']])
  colnames(coefs) <- c("Estimate","se","tval","pval")
  coefs$Sig<- ifelse(coefs$pval < alpha, 1,0)
  
  coefs<-t(data.frame(summary(M2)[[2]][[1]][['coefficients']][,c('Estimate','Pr(>|t|)')])) # Guarda el valor del coeficiente y el p_valor del mismo.
  rownames(coefs)<-c('coef','pval')
  
  
  M1$coef
  M1$secoef
  
  vars::causality(XY)
  
}
