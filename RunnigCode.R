#*********************************************************************#
#               Reading CPI data                                   ####
#.      Aggregates: headline, core, goods, etc...
#*********************************************************************#

#*********************************************************************#
# 1. Path definition ####
#*********************************************************************#

# Defina la ruta de la carpeta donde estan sus datos, codigos, etc
path_permaNent <- "C:/..."

# Define the folder where you stored the data
path_Datos <- paste(path_permaNent,"Data/",sep = "")
# Define the folder where you stored the codes
path.Code <- paste(path_permaNent,"Code/",sep = "")

# Name of directory to be created
Name_Newdir <- "Outputs/"
# To create folders to store results
dir.create(paste(path_permaNent,Name_Newdir,sep = ""))

# If you need to store data and plots create the respective folders
Name_Newdir2 <- "DataO/"
Name_Newdir3 <- "Plots/"
dir.create(paste(path_permaNent,Name_Newdir,Name_Newdir2,sep = ""))
dir.create(paste(path_permaNent,Name_Newdir,Name_Newdir3,sep = ""))

path.out <- paste(path_permaNent,Name_Newdir,sep = "")
path.out.datos <- paste(path.out,"DataO/",sep = "")
path.out.plots <- paste(path.out,"Plots/",sep = "")

#*********************************************************************#
# 2. Uploading libraries and source codes                          ####
#*********************************************************************#

# Libraries
source(paste(path.Code, "InsumoPaquetes.R", sep=""))

# Source functions
source(paste(path.Code, "Model1.R", sep=""))
source(paste(path.Code, "test_cointeg.R", sep=""))
source(paste(path.Code, "ARIMAX_CP.R", sep=""))
source(paste(path.Code, "BestARIMA_CP.R", sep=""))
source(paste(path.Code, "VECM_CP.R", sep=""))
source(paste(path.Code, "Funciones_Base0.R", sep=""))

#*********************************************************************#
# 3. Reading data                                                  ####
#*********************************************************************#

# Note: the index to upload should be aggregates
# Name of the file including the extension .xlsx
archivo_datos <- "DatosCO_ex.xlsx"  
# Name of the sheet
hoja_datos    <- "Hoja1"
# Define the reading range
rango_lectura <- "A1:C503"


# Reading the data
GB_I0 <- read_excel(paste(path_Datos,archivo_datos,sep=""),
                    sheet = hoja_datos,
                    range = rango_lectura)

# Define the name of the variables. The first should be dates
Names0 <- c("Fecha","Total","Core")
colnames(GB_I0) <- Names0

# Fitting the date format  %Y-%m-%d
GB_I0$Fecha <- parse_date_time(GB_I0$Fecha, orders = c("ymd", "dmy", "myd"))

# Define the aggregate to analyze and the column dates
agregado <- "Total"
Names0 <- c("Fecha", agregado)
GB_I0a <- GB_I0[,Names0]

# Check the uploaded data
head(GB_I0)

#*********************************************************************#
# 4. Selection of a values                                         ####
#*********************************************************************#

AA <- c(0.25,1,2,3,100)
names(AA) <- paste0("AA: ", AA)

#*********************************************************************#
# 5. Define the sample period and the dates to forecast out-of-sample #
#*********************************************************************#

# Initial date of the sample period
FechaIni <- "2010-01" #YYYY-MM

# Define the initial year to forecast
AO_INI_f <- 2018

# Define the final year to forecast
AO_FIN_f <- 2023

# Specify the number of months on the last forecast year
faltanM <- 12  

#*********************************************************************#
# 6. Define the forecast horizon                                   ####
#*********************************************************************#

hf <- 12

#*********************************************************************#
# 7. Define the models to fit                                      ####
#*********************************************************************#

# Available models: BAR, VECM_CP, VECM_CI, ARIMAX
modelos <- c("BAR", "ARIMAX", "VECM_CP", "VECM_CI")

#*********************************************************************#
# 8. Function to run the procedure according to the previous 
#     specifications  ####
#*********************************************************************#

start_time <- Sys.time()
GuardaResults <- f_inflains(agregado, GB_I0a, AA, FechaIni, AO_INI_f, 
                            AO_FIN_f, faltanM, hf, modelos)
end_time <- Sys.time()
end_time - start_time

#*********************************************************************#
# 9. Storing the results                                            ###
#*********************************************************************#

### Final date to forecast ####
FechaFinAll <- GuardaResults$FechaFinAll
FechaFinAll <- FechaFinAll[length(FechaFinAll)]

### Grafico IPC (agregado) ajustado estacionalmente e inflaciones instantaneas ####
GuardaResults$Grafico_InsObs

### Tabla inflaciones instantaneas ####
infIns_AA <- GuardaResults$infIns_AA
GuardaOut(FechaFinAll, agregado, "infIns_AA_CP_CO", infIns_AA, path.out.datos)
#data.frame(Fecha = index(GuardaResults$infIns_AA), coredata(GuardaResults$infIns_AA))

### Best ARIMA ####
if ("BAR" %in% modelos){

  # Forecasting best ARIMA (.txt)
  OutModBARIMApronos <- GuardaResults$OutModBARIMApronos
  GuardaOut(FechaFinAll, agregado, "BAR_pronos_CP_CO", OutModBARIMApronos, path.out.datos)
  
}

### ARIMAX ####
if ("ARIMAX" %in% modelos){
  
  # Forecasting ARIMAX model (.txt)
  OutModARXpronos <- GuardaResults$OutModARXpronos
  GuardaOut(FechaFinAll, agregado, "ARX_pronos_CP_CO", OutModARXpronos, path.out.datos)
  
}


### VECM_CP ####
if ("VECM_CP" %in% modelos){

  # Forecasting VECM_CP model (.txt)
  OutModVECMpronosB <- GuardaResults$OutModVECMpronosB
  GuardaOut(FechaFinAll, agregado, "VECM_pronos_CP_CO_B", OutModVECMpronosB, path.out.datos)
  
  # Order VECM_CP model (.txt)
  OutModVECMordenB <- GuardaResults$OutModVECMordenB
  GuardaOut(FechaFinAll, agregado, "VECM_orden_CP_CO_B", OutModVECMordenB, path.out.datos)
  
}

### VECM_CI ####
if ("VECM_CI" %in% modelos){

  # Forecasting VECM_CI model (.txt)
  OutModVECMpronosC <- GuardaResults$OutModVECMpronosC
  GuardaOut(FechaFinAll, agregado, "VECM_pronos_CI_CO_C", OutModVECMpronosC, path.out.datos)
  
  # Order VECM_CI model (.txt)
  OutModVECMordenC <- GuardaResults$OutModVECMordenC
  GuardaOut(FechaFinAll, agregado, "VECM_orden_CI_CO_C", OutModVECMordenC, path.out.datos)
  
} 


### Behavior of the lags according to the Coincident Profile analysis ####
Out2AA <- GuardaResults$Behave_a
GuardaOut(FechaFinAll, agregado, "Behave_a_CP_CO", Out2AA, path.out.datos)


