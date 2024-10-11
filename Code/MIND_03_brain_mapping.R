################################################################################
# Script to plot from .csv files the regional brain maps of: (1) MIND networks 
# and (2) effect sizes of MIND degree after stratifying by cognition and symptoms 
################################################################################

#  Copyright (C) 2024 University of Seville
# 
#  Written by Natalia García San Martín (ngarcia1@us.es)
# 
#  This file is part of Neurobiology MIND Psychosis toolkit.
# 
#  Neurobiology MIND Mapping toolkit is free software: 
#  you can redistribute it and/or modify it under the terms of the 
#  GNU General Public License as published by the Free Software Foundation, 
#  either version 3 of the License, or (at your option) any later version.
# 
#  Neurobiology MIND Mapping toolkit is distributed in the hope that 
#  it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Neurobiology MIND Mapping toolkit. If not, see 
#  <https://www.gnu.org/licenses/>.


rm(list=ls()) # Previous data cleaning in memory

location <- 'C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/2. MIND/'

source(paste0(location,'regional_brainmap_representation_borders.R'))
source(paste0(location,'regional_brainmap_representation.R'))


folders <-c('degree')

degrees <-c('degree_68_CN','degree_68_FEP')

# Select effsizes_morphometries as desired:
effsizes_morphometries <-c('degree_68')
# effsizes_morphometries <-c('degree_68_high_cognition','degree_68_low_cognition',
#                            'degree_68_high_BPRS','degree_68_low_BPRS',
#                            'degree_68_high_SAPS','degree_68_low_SAPS',
#                            'degree_68_high_SANS','degree_68_low_SANS')
effsizes_morphometries <-c('degree_68_low_high_88_cognition','degree_68_low_high_88_BPRS',
                           'degree_68_low_high_88_SAPS','degree_68_low_high_88_SANS',
                           'degree_68_low_high__cognition','degree_68_low_high__BPRS',
                           'degree_68_low_high__SAPS','degree_68_low_high__SANS')

effsizes_morphometries <- paste0('effsizes_', effsizes_morphometries)

saving_location = paste0(location,'Data/degree/plots/')


for (folder in folders){
  # Change file path as desired
  location = paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/',folder,'/')
  location = paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/2. MIND/Data/',folder,'/')
  
  setwd(location)
  
  
  # Libraries
  library(dplyr)
  library(ggseg)
  library(ggplot2)
  
  ############### MIND NETWORKS ###############
  
  morphometries <- c(degrees,effsizes_morphometries)
  
  for (i_tot in 1:length(morphometries)){
    morphometry <- morphometries[i_tot]
    print(morphometry)
  
  
    # We load MIND networks
    NEWDATA <- read.csv(paste0(location,morphometry,".csv"), row.names = 1)
    NEWDATA <- colMeans(NEWDATA,na.rm = TRUE)
    
    if (morphometry=='effsizes_degree_68'){
      Y_pred <-read.csv('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/CCAs/CCA_degrees/Y_pred_effsizes_degrees.csv', row.names = 1)
      map<-regional_brainmap_representation(t(Y_pred),paste0(morphometry,' predicted'),max(rbind(Y_pred,NEWDATA)),min(rbind(Y_pred,NEWDATA)),0)
      print(map)
      ggsave(paste0(saving_location,paste0(morphometry,'_predicted.png')), width = 8, height = 5, dpi = 300)
      
    }else if (morphometry=='degree_68_CN'){
      Y_pred <-read.csv('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/CCAs/CCA_degrees/Y_pred_CN_degrees.csv', row.names = 1)
      map<-regional_brainmap_representation(t(Y_pred),paste0(morphometry,' predicted'),max(rbind(Y_pred,NEWDATA)),min(rbind(Y_pred,NEWDATA)),0.17)
      print(map)
      ggsave(paste0(saving_location,paste0(morphometry,'_predicted.png')), width = 8, height = 5, dpi = 300)
    }

    # Configuration of the desired limits of representation (depends on whether it is for degree or effsizes)
    if (grepl('effsizes',morphometry) & grepl('_68',morphometry)){
      effsizes_degree_68 = read.csv(paste0(location,"effsizes_degree_68.csv"), row.names = 1)

      effsizes_degree_68_low_cognition = read.csv(paste0(location,"effsizes_degree_68_low_cognition.csv"), row.names = 1)
      effsizes_degree_68_high_cognition = read.csv(paste0(location,"effsizes_degree_68_high_cognition.csv"), row.names = 1)
      effsizes_degree_68_low_BPRS = read.csv(paste0(location,"effsizes_degree_68_low_BPRS.csv"), row.names = 1)
      effsizes_degree_68_high_BPRS = read.csv(paste0(location,"effsizes_degree_68_high_BPRS.csv"), row.names = 1)
      effsizes_degree_68_low_SAPS = read.csv(paste0(location,"effsizes_degree_68_low_SAPS.csv"), row.names = 1)
      effsizes_degree_68_high_SAPS = read.csv(paste0(location,"effsizes_degree_68_high_SAPS.csv"), row.names = 1)
      effsizes_degree_68_low_SANS = read.csv(paste0(location,"effsizes_degree_68_low_SANS.csv"), row.names = 1)
      effsizes_degree_68_high_SANS = read.csv(paste0(location,"effsizes_degree_68_high_SANS.csv"), row.names = 1)
     
      effsizes_degree_68_low_high_cognition = read.csv(paste0(location,"effsizes_degree_68_low_high_88_cognition.csv"), row.names = 1)
      effsizes_degree_68_low_high_BPRS = read.csv(paste0(location,"effsizes_degree_68_low_high_88_BPRS.csv"), row.names = 1)
      effsizes_degree_68_low_high_SAPS = read.csv(paste0(location,"effsizes_degree_68_low_high_88_SAPS.csv"), row.names = 1)
      effsizes_degree_68_low_high_SANS = read.csv(paste0(location,"effsizes_degree_68_low_high_88_SANS.csv"), row.names = 1)

      total_effsizes_68 = cbind(t(effsizes_degree_68),t(effsizes_degree_68_low_cognition),t(effsizes_degree_68_high_cognition),t(effsizes_degree_68_low_BPRS),
                                t(effsizes_degree_68_high_BPRS),t(effsizes_degree_68_low_SAPS),t(effsizes_degree_68_high_SAPS),
                                t(effsizes_degree_68_low_SANS),t(effsizes_degree_68_high_SANS))
      total_effsizes_68 = cbind(t(effsizes_degree_68_low_high_cognition),t(effsizes_degree_68_low_high_BPRS),t(effsizes_degree_68_low_high_SAPS),
                                t(effsizes_degree_68_low_high_SANS))
      limits <- c(min(total_effsizes_68), max(total_effsizes_68))
      
    } else{
      degree_CN = read.csv(paste0(location,"degree_68_CN.csv"), row.names = 1)
      degree_FEP = read.csv(paste0(location,"degree_68_FEP.csv"), row.names = 1)
      total_morphometries = rbind(degree_FEP, degree_CN)
      limits <- c(min(total_morphometries), max(total_morphometries))
    }
      
    
    # Map the data
    if (grepl('effsizes',morphometry)){
      if (grepl('low_high',morphometry)){
        significant_NEWDATA <- read.csv(paste0(location,morphometry,"_ttest_significant.csv"), row.names = 1)
      }else{
        significant_NEWDATA <- read.csv(paste0(location,morphometry,"_significant.csv"), row.names = 1)
      }
      data<-data.frame(values=NEWDATA,sig=t(!is.na(significant_NEWDATA)))
      names(data)[names(data) == "effsizes"] <- "sig"
      if (morphometry=='effsizes_degree_68'){
        limits <- c(min(rbind(Y_pred,NEWDATA)),max(rbind(Y_pred,NEWDATA)))
      }
      map<- regional_brainmap_representation_borders(data,morphometry,limits[2],limits[1],0,'vertical')
      print(map)
      ggsave(paste0(saving_location,paste0(morphometry,'_significant.png')), width = 8, height = 5, dpi = 300)
    
    } else {
      
      if (morphometry=='degree_68_CN'){
        limits <- c(min(rbind(Y_pred,NEWDATA)),max(rbind(Y_pred,NEWDATA)))

      }
      map<-regional_brainmap_representation(NEWDATA,morphometry,limits[2],limits[1],0.17)
      print(map)
      ggsave(paste0(saving_location,paste0(morphometry,'.png')), width = 8, height = 5, dpi = 300)
    }
    browser()
  }
  
}


# FEP centiles

location <- 'C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/2. MIND/Data/centiles/'

centiles <- read.csv(paste0(location,"centiles_CorticalMeasuresENIGMA_GrayAvg.csv"), row.names = 1)
centiles<-centiles[centiles$session==1 & centiles$dx=='PS',18:85]
significant_centiles<-read.csv(paste0(location,"significant_values_68.csv"),row.names=1)

data<-data.frame(values=colMeans(centiles),sig=!is.na(significant_centiles[,2]))
map<-regional_brainmap_representation_borders(data,'FEP centiles',max(data$values),min(data$values),0.5,'vertical')
print(map)

'Done!'

