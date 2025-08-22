test_cointeg<-function(j,series,CI,ecdet,exogen,alpha_coint,summary=c(TRUE,FALSE)){
  #j: Numero de  modelo
  #CI: Rezagos a incluir
  #ecdet: inclusión de componentes deterministicos 
  #exogen: vars dummies 
  #alpha_coint:percentil (1,5,10 pct) para evaluar estadistica
  sjf.vecm_traza <- ca.jo(series, type="trace", ecdet=ecdet, K=CI, spec="transitory",dumvar=exogen) # Se realiza el procedimiento de Johansen bajo el test de la traza y una especificaci?n del modelo VEC transitoria. Los an?lisis que ?ste proporciona se guardan en sjf.vecm_traza
  O_traza=summary(sjf.vecm_traza) # El resumen del an?lisis del proceso de Johansen se guarda en O_traza
  sjf.vecm_eigen <- ca.jo(series, type="eigen", ecdet=ecdet, K=CI, spec="transitory",dumvar=exogen) # Se realiza el procedimiento de Johansen bajo el test delm?ximo valorpropio y una especificaci?n del modelo VEC transitoria. Los an?lisis que ?ste proporciona se guardan en sjf.vecm_eigen
  O_eigen=summary(sjf.vecm_eigen) # El resumen del an?lisis del proceso de Johansen se guarda en O_eigen
  
  # Selecci?n de n?mero de relaciones de cointegraci?n test de la traza
  fila_acepta_t<-max(which(O_traza@teststat<O_traza@cval[,alpha_coint]))  # Guarda en fila_acepta_t la posici?n en la que no se rechaza que hay cointegraci?n seg?n el test de la traza
  if(fila_acepta_t==-Inf){ # Si nunca deja de rechazarse (no hay relacion de cointegracion) entonces
    rt<-ncol(series) # se establece que el n?mero de relaciones de cointegraci?n es igual al n?mero de variables
  } else{ # Si deja de rechazar en alg?n momento entonces
    x<-rownames(O_traza@cval)[fila_acepta_t] # Se guarda en x el n?mero de relaciones de cointegraci?n resultante
    rt<-as.numeric(sapply(strsplit(as.character(x), " "),unlist)[3]) # Se convierte en formato num?rico el valor de x encontrado anteriormente asociado aln?mero de relaciones de cointegraci?n y lo guarda en rt
  }
  
  # Slecci?n de n?mero de relaciones de cointegraci?n test eigen
  fila_acepta_e<-max(which(O_eigen@teststat<O_eigen@cval[,alpha_coint])) # Guarda en fila_acepta_e la posici?n en la que no se rechaza que hay cointegraci?n seg?n el test eigen
  if(fila_acepta_e==-Inf){ # Si nunca deja de rechazarse entonces
    re<-ncol(series) # se establece que el n?mero de relaciones de cointegraci?n es igual al n?mero de variables
  } else{ # Si deja de rechazar en alg?n momento entonces
    x<-rownames(O_eigen@cval)[fila_acepta_e] # Se guarda en x el n?mero de relaciones de cointegraci?n resultante
    re<-as.numeric(sapply(strsplit(as.character(x), " "),unlist)[3]) # Se convierte en formato num?rico el valor de x encontrado anteriormente asociado aln?mero de relaciones de cointegraci?n y lo guarda en re
  }
  
  if(rt!=0 & re!=0) r=min(rt,re) else r=max(rt,re) # Si hay cointegraci?n en los dos test, se establece el n?mero de relaciones menor que se encontr? en los dos test y se guarda en r
  if(r==rt) test_type="trace" else if(r==re) test_type="eigen" # Se selecciona como tipo de test "test_type" el test que arroj? menor n?mero de relaciones de cointegraci?n con base en lo concluido en la l?nea anterior
  
  outp_tci <- list("modelo"=paste0("modelo",j),"lag"=CI,"r"=as.numeric(r), "test_type" = test_type) # Guarda en outp_tci el nombre del modelo, el n?mero de rezagos de su representaci?n VAR, eln?mero de relaciones determinadas y el test 
  
  if(summary==TRUE){ # Si se desea hacer un resumen general delas relaciones de cointegraci?n seg?n varios tipos de modelo entonces
    
    sfj.vecm_none<-ca.jo(series, type =test_type, K =CI, ecdet="none", spec="transitory", dumvar=exogen); O_none=summary(sfj.vecm_none) # Se realiza el procedimiento de Johansen bajo un tipo de modelo sin constante el test determinado anteriormente (el que cumple tener menor n?mero de relaciones de cintegraci?n) y una especificaci?n del modelo VEC transitoria. Los an?lisis que ?ste proporciona se guardan en sfj.vecm_none
    sfj.vecm_const<-ca.jo(series, type =test_type, K =CI, ecdet="const", spec="transitory", dumvar=exogen); O_const=summary(sfj.vecm_const)  # Se realiza el procedimiento de Johansen bajo un tipo de modelo con constante el test determinado anteriormente (el que cumple tener menor n?mero de relaciones de cintegraci?n) y una especificaci?n del modelo VEC transitoria. Los an?lisis que ?ste proporciona se guardan en sfj.vecm_const
    sfj.vecm_trend<-ca.jo(series, type =test_type, K =CI, ecdet="trend", spec="transitory", dumvar=exogen); O_trend=summary(sfj.vecm_trend) # Se realiza el procedimiento de Johansen bajo un tipo de modelo con tendencia el test determinado anteriormente (el que cumple tener menor n?mero de relaciones de cintegraci?n) y una especificaci?n del modelo VEC transitoria. Los an?lisis que ?ste proporciona se guardan en sfj.vecm_trend
    
    summ_teststat<-cbind(O_none@teststat,O_const@teststat,O_trend@teststat) # Guarda en summ_teststat los estad?sticos estimados para los tres test
    summ_cval<-cbind(O_none@cval[,alpha_coint],O_const@cval[,alpha_coint],O_trend@cval[,alpha_coint]) # Guarda en summ_cval los valores cr?ticos para los tres tipos de test
    
    sfila_acepta=NULL;sx=NULL;sr=NULL # Crea objetos nulos para guardar el n?mero de relaciones de cointegraci?n por cada test
    for(i in 1:3){ #Loop que recorre las columnas de las matrices que contienen los estad?sticos y los valores cr?ticos 
      sfila_acepta[i]<-max(which(summ_teststat[,i]<summ_cval[,i])) # Selecciona la fila en la que deja de rechazarse el i-?simo test
      sx[i]<-rownames(summ_cval)[sfila_acepta[i]]  # Se guarda en sx el n?mero de relaciones de cointegraci?n resultante
      sr[i]<-as.numeric(sapply(strsplit(as.character(sx[i]), " "),unlist)[3]) # Se convierte en formato num?rico el valor de sx encontrado anteriormente asociado aln?mero de relaciones de cointegraci?n y lo guarda como elemento i-?simo de sr
    }
    df.sr <- as.data.frame(as.numeric(sr)); rownames(df.sr) <- c("none","const","trend"); colnames(df.sr)<-"r" # Se reunen en una matriz eln?mero de relaciones de cointegraci?n concluido por los tres tes (un test por cada tipo de modelo, son 3 tipos de modelo)
    outp_tci_s <- list("modelo"=paste0("modelo",j),"lag"=CI,"r"=t(df.sr)) # Se guarda en outp_tci_s una lista que informa elnombre delmodelo, el n?mero de rezagos de la representaci?n VAR y la matriz de relaciones de cointegraci?n creada en la l?nea anterior
    
  }
  
  if(summary==FALSE) return(list(outp_tci,list())) else if(summary==TRUE) return(list(outp_tci,outp_tci_s)) # La salida de esta funci?n es la lista outp_tci el nombre del modelo, el n?mero de rezagos de su representaci?n VAR, eln?mero de relaciones determinadas y el test 
  # ? si se pidi? un resumen general se proporciona adicionalmente, la lista outp_tci_s una lista que informa elnombre delmodelo, el n?mero de rezagos de la representaci?n VAR y la matriz de relaciones de cointegraci?n creada en la l?nea anterior
}
