

# Fit an ARIMAX model beteen a answer (x1) and a predictor (x2) and its
#  lags according to the the coincident profile analysis.

# The output is a lits with two matrices (parameters and forecast) 

ARIMAX.fit <- function(FechaIni, FechaFin, Nombre1, x1,x2,LagA,a, hf,GB_Obs){

  # FechaIni: Fecha inicial de la muestra 
  # FechaFin: Fecha final de la muestra 
  # Nombre1: Nombre de la variable de an?lisis 
  # x1: datos de la variable respuesta, vector
  # x2: datos de la variable respuesta, vector
  # LagA: Rezago encontrado segun perfil coincidente
  # hf: Valor del horizonte de pronostico
  
  
  #### Ajuste ARIMAX
  Lags <- abs(LagA)+1
  
  XY <- matrix(x2,ncol = 1)
  XY <- embed(XY, Lags)
  YY <- x1[-c(1:(Lags-1))]
  #if(Lags > 1) YY <- x1[-c(1:(Lags-1))] else YY <- x1[-c(1:(Lags))]
  
  #fit0 <- auto.arima(YY) 
  
  fit01 <- auto.arima(YY, xreg = XY)
  orden0 <- OrdenArima(fit01)
  orden <- orden0$orden
  ordenS <- c(orden0$ordenS, orden0$Period)
  
  ResidualsCheck <- Residuals(fit01$residuals)
  
  Params01 <- data.frame(FechaIni, FechaFin, Nombre1, coef = fit01$coef, se = sqrt(diag(fit01$var.coef)))
  colnames(Params01) <- c("FechaI", "FechaF","Agregado", "coef","se")
  Params01$Sig <- ifelse( abs(Params01$coef) > 2*Params01$se,1, 0)
  Params01$Npar <- 1:nrow(Params01)

  ########################################################################
  #### Pronos
  ########################################################################
  fFin <- which(substr(GB_Obs$Fecha,1,7) == FechaFin)
  y1.Obs <- t(GB_Obs[(fFin+1):(fFin+hf),Nombre1])
  colnames(y1.Obs) <- paste("Obs",1:hf,sep = "")
  
  BaseF0 <- XY[(nrow(XY)-abs(LagA)):nrow(XY),]
  BaseF1 <- matrix(rep(XY[nrow(XY),],each=hf-Lags), ncol = Lags)
  fitX01 <- as.matrix(rbind(BaseF0, BaseF1))
  
  fc01 <- data.frame(forecast::forecast(fit01, h=hf, xreg=fitX01))
  ForesC <- data.frame(t(fc01[,1]))
  colnames(ForesC) <- paste("h",1:hf,sep = "")
  
  FitOrd <- data.frame(t(orden), t(ordenS))
  colnames(FitOrd) <- c("p","d","q","P","D","Q", "S")
    
  Out0 <- data.frame(FechaIni = FechaIni,  
                FechaFin = FechaFin, Agregado = Nombre1, 
                MaxLT = LagA,
                MaxL = Lags-1,
                a=a, FitOrd,
                t(ResidualsCheck), y1.Obs, ForesC)
  
  Output <- list(ParametrosE = Params01, Pronosticos = Out0)
  return(Output)
}