# ####################################################################### #
#                           CARGA DE LIBRERIAS                      ####
# ###################################################################### ##

# Lista de paquetes que se desean cargar
load.lib <- c("readxl", "multDM", "gtools", "plyr", "mvnTest", "rlist",
              "tsDyn", "MTS", "vars", "forecast", "portes", "urca",
              "imputeTS", "dplyr", "bvarsv", "tidyverse", "aTSA", 
              "data.table", "tseries", "mFilter", "DescTools", "svars",
              "ggplot2", "moments", "mlVAR", "jcolors", "tsDyn", "sarima",
              "sysid","lawstat","FinTS", "xtable", "glmnet","randomForest",
              "astsa","Coinprofile", "seasonal", "highcharter", "writexl", "lubridate", "xts")

# Identificar los paquetes no instalados
install.lib <- load.lib[!load.lib %in% installed.packages()]

# Imprimir los nombres de los paquetes no instalados o confirmar que todos están instalados
if (length(install.lib) > 0) {
  cat("Los siguientes paquetes no estan instalados:\n")
  cat(install.lib, sep = "\n")
  
  # Instalar los paquetes no instalados
  install.packages(install.lib)
} else {
  cat("Todos los paquetes estan instalados.\n")
}

# Cargar las librerias requeridas
sapply(load.lib, require, character = TRUE)


###################################################################### #
# Funciones que adaptan la salida de los datos para plataforma   ####
###################################################################### #
GuardaOutExcel <- function(FechaFinAll, NombreOut, Grupo, BaseG, path.out){
  IniF <- FechaFinAll[1]
  FinF <- FechaFinAll[length(FechaFinAll)]
  file.out <- paste(NombreOut, "_", Grupo, "_",IniF,"_",FinF, ".xlsx",sep = "")
  write_xlsx(BaseG, paste(path.out, file.out, sep = ""))
}

f_pronos <- function(nombreDF, DF, hf){ 
  selectcol <- c("FechaFin",paste("h", 1:hf, sep =""),"Obs1","a")
  DF2 <- DF %>% dplyr::select(selectcol)
  #DF2 <- DF[,selectcol] 
  for (pronos in 1:hf){
    DF2[,pronos+1] <- lag(DF2[,pronos+1], n=pronos)
    if (pronos==hf) {DF2[,pronos+2] <- lag(DF2[,pronos+2], n=1)}}
  colnames(DF2) <- c("Fecha", paste0("pron_", 1:hf), "inf_obs", "a")
  DF2$Fecha <- as.Date(paste(DF2$Fecha, "-1", sep = ""))
  DF2$Fecha <- format(DF2$Fecha,format = "%d/%m/%Y")
  if(nombreDF=="ARIMAX"){DF2$model <- rep("ARIMAX", nrow(DF2))} 
  else if (nombreDF=="BARIMA"){DF2$model <- rep("BARIMA", nrow(DF2))}
  else if (nombreDF=="VECM_CP"){DF2$model <- rep("VECM_CP", nrow(DF2))}
  else if (nombreDF=="VECM_CI"){DF2$model <- rep("VECM_CI", nrow(DF2))}
  # Eliminar N'as y organizar salida
  DF2 <- DF2[-c(1:hf),c(1,ncol(DF2), 2:(ncol(DF2)-1))] #-c(1:12)
  #DF2 <- DF2 %>% filter(complete.cases(.))
  return(DF2)
} 


###################################################################### #
# Funcion que ejecuta los modelos segun especificaciones  ####
###################################################################### #

f_inflains <- function(agregado, GB_I0a ,AA, FechaIni, AO_INI_f, AO_FIN_f, faltanM, hf, modelos){
  
  # Variables necesarias para el proceso no definidas en argumentos
  Names0 <- names(GB_I0a)
  #agregado <- Names0[2]
  
  # Ver progreso del vector a corrido en el for
  #ncat <- floor(length(AA)/10)
  #progreso <- 0
  
  ###################################################################### #
  # Ajuste por estacionalidad ####
  ###################################################################### #
  #if(FitAS){
    GB_saj <- GB_I0a
    year <- as.numeric(substr(GB_I0a$Fecha[1],1,4))
    month0 <- as.numeric(substr(GB_I0a$Fecha[1],6,7))
    frequ <- 12
    j=2
    Vtemp <- ts(GB_I0a[,j], start = c(year,month0), 
                frequency = frequ)
    
    m <- seas(Vtemp, transform.function = "auto")
    
    # serie ajustada por estacionalidad
    GB_saj[,j] <- final(m)
    GB_saj <- data.frame(GB_saj)
    #GB_saj$Fecha <- as.character(GB_saj$Fecha)
    GB_VM <- GB_saj[-1,]
  #}else GB_VM <- GB_I0a[-1,]
  
  ###################################################################### #
  # Variaciones mensuales ####
  ###################################################################### #
  
  for(j in 2:ncol(GB_saj)) GB_VM[,j] <- GB_saj[2:nrow(GB_saj),j] / GB_saj[1:(nrow(GB_saj)-1),j]
  
  ###################################################################### #
  # Inflacion observada  ####
  ###################################################################### #
  
  GB_Obs <- GB_VM[-c(1:11),]
  for(j in 2:ncol(GB_Obs)){
    for(i in 1:nrow(GB_Obs)){
      GB_Obs[i,j] <- InflaIns(GB_VM[c(i:(i+11)),j],rep(1,12))
    }}
  
  ########################################################################## #
  # Creacion de objetos que guarden las salidas de los modelos por cada a ####
  ########################################################################## #
  # BEST ARIMA
  if ("BAR"  %in% modelos)    OutModBARIMApara <- OutModBARIMApronos <- PlatfModBARIMApronos <- NULL
  # ARIMAX
  if ("ARIMAX"  %in% modelos) OutModARXpara    <- OutModARXpronos    <- PlatfModARXpronos    <- NULL
  # VECM_CP
  if ("VECM_CP" %in% modelos) OutModVECMparaB  <- OutModVECMpronosB  <- OutModVECMordenB     <- PlatfModVECMpronosB <- NULL
  # VECM_CI
  if ("VECM_CI" %in% modelos) OutModVECMparaC  <- OutModVECMpronosC  <- OutModVECMordenC     <- PlatfModVECMpronosC <- NULL
  #OutModVECMparaA <- OutModVECMpronosA <- NULL
  
  # Rezago optimo
  Out2AA  <- NULL
  
  # Creacion para resultados funcion y plataforma 
  resultados <- list()
  
  # Almacenar las inflaciones instantaneas
  AA_out <- NULL
  
  for(a in AA){
    ###################################################################### #
    # Inflacion instantanea  ####
    ###################################################################### #
    kapa <- KernelA(a,12)
    GB_I <-  GB_VM[-c(1:11),]
    for(j in 2:ncol(GB_I)){
      for(i in 1:nrow(GB_I)){
        GB_I[i,j] <- InflaIns(GB_VM[c(i:(i+11)),j],sort(kapa))
      }}
    colnames(GB_I) <- paste(Names0,"_Ins",sep="")
    GB_F <- cbind(GB_Obs, GB_I[,-1])
    colnames(GB_F) <- c(colnames(GB_Obs), colnames(GB_I)[-1])
    
    ###################################################################### #
    # Variables auxiliares para visualizacion ####
    ###################################################################### #
    tvar <-paste("^", colnames(GB_F)[2],sep=""); name1 <- c(grep(tvar, colnames(GB_F)))
    frequ <- 12
    #a_out  <- xts(x = GB_F[,name1[2]], order.by = GB_F$Fecha)
    #AA_out <- cbind(AA_out, a_out)
    AA_out <- cbind(AA_out, GB_F[,name1[2]])
    
    
    ############################################################################## #
    # Creacion de objetos que guarden las salidas de los modelos por cada fecha ####
    ############################################################################## #
    # ARIMAX
    if ("ARIMAX"  %in% modelos) ModARIMAX.Para  <- ModARIMAX.Pronos <- RezOptimo        <- NULL 
    # VECM_CP
    if ("VECM_CP" %in% modelos) ModVECM.ParaB   <- ModVECM.PronosB  <- ModVECM.OrdenB   <- NULL
    # VECM_CI
    if ("VECM_CI" %in% modelos) ModVECM.ParaC   <- ModVECM.PronosC  <- ModVECM.OrdenC   <- NULL
    #ModVECM.ParaA <- ModVECM.PronosA <- NULL
    
    ###################################################################### #
    # Perfil coincidente ####
    ###################################################################### #
    
    FechaFinAll <- c(paste(c(rep(AO_INI_f:AO_FIN_f,each=12)),c("01","02","03",
                                                               "04","05","06",
                                                               "07","08","09",
                                                               "10","11","12"),sep="-"))
    FechaFinAll <- c(paste(AO_INI_f-1,"12",sep="-"),FechaFinAll)
    FechaFinAll <- FechaFinAll[ c(1:(length(FechaFinAll)-(12-faltanM)) )]
    
    for(FechaFin in FechaFinAll){
      fIni <- which(substr(GB_Obs$Fecha,1,7) == FechaIni)
      fFin <- which(substr(GB_Obs$Fecha,1,7) == FechaFin)
      
      GB_F <- cbind(GB_Obs[fIni:fFin,], GB_I[fIni:fFin,-1])
      colnames(GB_F) <- c(colnames(GB_Obs), colnames(GB_I)[-1])
      name1 <- c(grep(tvar, colnames(GB_F)))
      
      yearIni <- as.numeric(substr(GB_F$Fecha[1],1,4))
      monthIni <- as.numeric(substr(GB_F$Fecha[1],6,7))
      yearFin <- as.numeric(substr(GB_F$Fecha[nrow(GB_F)],1,4))
      monthFin <- as.numeric(substr(GB_F$Fecha[nrow(GB_F)],6,7))
      frequ <- 12
      
      x1 <- GB_F[,name1[1]]
      x2 <- GB_F[,name1[2]]
      
      Out.profile <- try(coincident_profile(x2,x1, frequ,12,nvar1 = "Ins",nvar2 = "Obs",
                                            print.graf = F,iyear = yearIni, lyear = yearFin,
                                            imonth = 12-monthIni+1, lmonth = monthFin,
                                            tit1 = paste(colnames(GB_F)[name1[1]]," Ins",sep=""),
                                            tit2 = paste(colnames(GB_F)[name1[1]]," Obs",sep=""),
                                            tit3 = "Ins vs. Obs"),TRUE)
      
      if(class(Out.profile)[1] != "try-error"){
        LagA_Orig <- LagA <- as.numeric(Out.profile$MainLag[1])
      }else LagA_Orig <- LagA <- 0 
      #LagA_Orig <- LagA <- -6
      #if(LagA == 0) LagA = -1 else LagA = LagA 
      if(LagA >= 0){
        message("There is no leading index")
        LagA = -1
      }
      
      ###################################################################### #
      # Salidas de los modelos por cada fecha ####
      ###################################################################### #
      
      # BEST ARIMA
      if("BAR" %in% modelos && a == AA[1] ){ #Unica vez
        ModBARIMA.0        <- BestARIMA.fit(FechaIni, FechaFin, agregado, x1, hf, GB_Obs, a)
        OutModBARIMApara   <- rbind(OutModBARIMApara, ModBARIMA.0$ParametrosE)
        OutModBARIMApronos <- rbind(OutModBARIMApronos, ModBARIMA.0$Pronosticos)
      }
      
      # ARIMAX
      if ("ARIMAX" %in% modelos){
        ModARIMAX.0      <- ARIMAX.fit(FechaIni, FechaFin, agregado, x1,x2,LagA,a, hf, GB_Obs)
        ModARIMAX.Para   <- rbind(ModARIMAX.Para, ModARIMAX.0$ParametrosE)
        ModARIMAX.Pronos <- rbind(ModARIMAX.Pronos, ModARIMAX.0$Pronosticos)
      }
      
      # VECM_CP
      if ("VECM_CP" %in% modelos){
        ModVECM.0b <- VECMCP.fit2(FechaIni, FechaFin, agregado, x1,x2,LagA,a, hf, GB_Obs,alpha = 0.05)
        if(!any(is.na(ModVECM.0b))){
          ModVECM.ParaB   <- rbind(ModVECM.ParaB, ModVECM.0b$ParametrosE)
          ModVECM.PronosB <- rbind(ModVECM.PronosB, ModVECM.0b$Pronosticos)
        }
        ModVECM.OrdenB <- rbind(ModVECM.OrdenB, ModVECM.0b$Orden)
      }
      
      # VECM_CI
      if ("VECM_CI" %in% modelos){
        ModVECM.0c <- VECMCP.fit3(FechaIni, FechaFin, agregado, x1,x2,LagA,a, hf, GB_Obs,alpha = 0.05)
        if(!any(is.na(ModVECM.0c))){
          ModVECM.ParaC   <- rbind(ModVECM.ParaC, ModVECM.0c$ParametrosE)
          ModVECM.PronosC <- rbind(ModVECM.PronosC, ModVECM.0c$Pronosticos)
        }
        ModVECM.OrdenC    <- rbind(ModVECM.OrdenC, ModVECM.0c$Orden)
      }
      
      # Rezago Optimo
      RezOptimo <- rbind(RezOptimo, data.frame(FechaIni,FechaFin,agregado, a, LagA, LagA_Orig))
      
    } #for Fecha
    
    ############################################################################## #
    # Salidas de los modelos por cada a y ajuste de pronosticos para plataforma ####
    ############################################################################## #
    # BEST ARIMA
    if ("BAR" %in% modelos && a == AA[1] ){
      aux_PlatfModBARIMApronos <- f_pronos("BARIMA", OutModBARIMApronos, hf)
      PlatfModBARIMApronos <- rbind(PlatfModBARIMApronos, aux_PlatfModBARIMApronos)
    }
    
    # ARIMAX
    if ("ARIMAX" %in% modelos){
      OutModARXpara         <- rbind(OutModARXpara, ModARIMAX.Para)
      OutModARXpronos       <- rbind(OutModARXpronos, ModARIMAX.Pronos)
      aux_PlatfModARXpronos <- f_pronos("ARIMAX", ModARIMAX.Pronos, hf)
      PlatfModARXpronos     <- rbind(PlatfModARXpronos, aux_PlatfModARXpronos)
    }
    
    # VECM_CP
    if ("VECM_CP" %in% modelos){
      OutModVECMparaB         <- rbind(OutModVECMparaB, ModVECM.ParaB)
      OutModVECMpronosB       <- rbind(OutModVECMpronosB, ModVECM.PronosB)
      OutModVECMordenB        <- rbind(OutModVECMordenB, ModVECM.OrdenB)
      aux_PlatfModVECMpronosB <- f_pronos("VECM_CP", ModVECM.PronosB, hf)
      PlatfModVECMpronosB     <- rbind(PlatfModVECMpronosB, aux_PlatfModVECMpronosB)
    }
    
    # VECM_CI
    if ("VECM_CI" %in% modelos){
      OutModVECMparaC <- rbind(OutModVECMparaC, ModVECM.ParaC)
      OutModVECMpronosC <- rbind(OutModVECMpronosC, ModVECM.PronosC)
      OutModVECMordenC <- rbind(OutModVECMordenC, ModVECM.OrdenC)
      aux_PlatfModVECMpronosC <- f_pronos("VECM_CI", ModVECM.PronosC, hf)
      PlatfModVECMpronosC     <- rbind(PlatfModVECMpronosC, aux_PlatfModVECMpronosC)
    }
    
    # Rezago optimo
    Out2AA <- rbind(Out2AA, RezOptimo)
    
    ###################################################################### #
    # Progreso del for  ####
    ###################################################################### #
    #progreso = progreso + 1
    #if (progreso%%ncat == 0) 
    #  cat(100*round(progreso/length(AA),1), "% del vector AA corrido ... \n", sep = "")
  } #for AA
  
  
  ###################################################################### #
  # Visualizar GB_Obs y la inflacion instantanea (Variando a) ####
  ###################################################################### #
  colnames(AA_out) <- paste0("a: ", AA)
  graf_Inf_IO <- hchart(xts(x = GB_Obs[,agregado], order.by = GB_Obs$Fecha), name = agregado, color = "black")
  
  for (i in 1:length(AA)){
    graf_Inf_IO <- graf_Inf_IO %>% 
      #hc_add_series(AA_out[,i], type = "line", name = names(AA_out)[i])
      hc_add_series(xts(x = AA_out[,i], order.by = GB_Obs$Fecha), type = "line", name = colnames(AA_out)[i])
  }
  graf_Inf_IO <- graf_Inf_IO %>% hc_add_theme(hc_theme_google()) %>% 
    hc_exporting(enabled = TRUE, filename = "FileName" ) %>% 
    hc_title(text = paste("IPC", agregado, "ajustado estacionalmente e inflaciones instantáneas", sep = " ")) %>% 
    hc_tooltip(valueDecimals = 2, borderWidth = 2 , 
               crosshairs = TRUE, table = TRUE) %>%
    hc_legend(enabled = TRUE) %>%
    hc_colors(RColorBrewer::brewer.pal(n = 9, "Set1"))
  
  
  ###################################################################### #
  # Guardar salidas de cada modelo ####
  ###################################################################### #
  # BEST ARIMA
  if ("BAR" %in% modelos){
    resultados$OutModBARIMApara      <- OutModBARIMApara
    resultados$OutModBARIMApronos    <- OutModBARIMApronos
    resultados$PlatfModBARIMApronos  <- PlatfModBARIMApronos
  }
  
  # ARIMAX
  if ("ARIMAX" %in% modelos){
    resultados$OutModARXpara     <- OutModARXpara
    resultados$OutModARXpronos   <- OutModARXpronos
    resultados$PlatfModARXpronos <- PlatfModARXpronos
  }
  
  # VECM_CP
  if ("VECM_CP" %in% modelos){
    resultados$OutModVECMparaB     <- OutModVECMparaB
    resultados$OutModVECMpronosB   <- OutModVECMpronosB
    resultados$OutModVECMordenB    <- OutModVECMordenB
    resultados$PlatfModVECMpronosB <- PlatfModVECMpronosB 
  }
  
  # VECM_CI
  if ("VECM_CI" %in% modelos){
    resultados$OutModVECMparaC     <- OutModVECMparaC
    resultados$OutModVECMpronosC   <- OutModVECMpronosC
    resultados$OutModVECMordenC    <- OutModVECMordenC
    resultados$PlatfModVECMpronosC <- PlatfModVECMpronosC
  }
  
  # Rezago optimo
  resultados$Behave_a <- Out2AA
  
  # DF Inflacion instantanea
  infIns_AA <- cbind(GB_Obs, AA_out)
  colnames(infIns_AA) <- c(colnames(GB_Obs), colnames(AA_out))
  resultados$infIns_AA <- infIns_AA
  resultados$Grafico_InsObs <- graf_Inf_IO
  resultados$FechaFinAll <- FechaFinAll#[length(FechaFinAll)]
  return(resultados)
  
} # Funcion f_inflains
