

# Fit an ARIMAX model beteen a answer (x1) and a predictor (x2) and its
#  lags according to the the coincident profile analysis.

# The output is a lits with two matrices (parameters and forecast) 

BestARIMA.fit <- function(FechaIni, FechaFin, Nombre1, x1, hf, GB_Obs, a){
  
  # FechaIni: Fecha inicial de la muestra 
  # FechaFin: Fecha final de la muestra 
  # Nombre1: Nombre de la variable de an?lisis 
  # x1: datos de la variable respuesta, vector
  # hf: Valor del horizonte de pronostico
  
    YY <- x1
  
  fit01 <- auto.arima(YY) 
  #fit01 <- auto.arima(YY, xreg = XY)
  orden0 <- OrdenArima(fit01)
  orden <- orden0$orden
  if(orden[1]==0 && orden[3]==0){
    coefs <- 0
    SeCoefs <- 0
  }else{
    coefs <- fit01$coef
    SeCoefs <- sqrt(diag(fit01$var.coef))
  }
  ordenS <- c(orden0$ordenS, orden0$Period)
  
  ResidualsCheck <- Residuals(fit01$residuals)
  
  Params01 <- data.frame(FechaIni, FechaFin, Nombre1, coef = coefs, se = SeCoefs)
  colnames(Params01) <- c("FechaI", "FechaF","Agregado", "coef","se")
  Params01$Sig <- ifelse( abs(Params01$coef) > 2*Params01$se,1, 0)
  Params01$Npar <- 1:nrow(Params01)
  
  ########################################################################
  #### Pronos
  ########################################################################
  fFin <- which(substr(GB_Obs$Fecha,1,7) == FechaFin)
  y1.Obs <- t(GB_Obs[(fFin+1):(fFin+hf),Nombre1])
  colnames(y1.Obs) <- paste("Obs",1:hf,sep = "")
  
  fc01 <- data.frame(forecast::forecast(fit01, h=hf))
  ForesC <- data.frame(t(fc01[,1]))
  colnames(ForesC) <- paste("h",1:hf,sep = "")
  
  FitOrd <- data.frame(t(orden), t(ordenS))
  colnames(FitOrd) <- c("p","d","q","P","D","Q", "S")
  
  Out0 <- data.frame(FechaIni = FechaIni,  
                     FechaFin = FechaFin, Agregado = Nombre1, 
                     MaxLT = NA,
                     MaxL = NA,
                     a=a, FitOrd,
                     t(ResidualsCheck), y1.Obs, ForesC)
  
  Output <- list(ParametrosE = Params01, Pronosticos = Out0)
  return(Output)
}