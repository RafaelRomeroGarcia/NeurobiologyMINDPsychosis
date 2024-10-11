################################################################################
# Script to plot from .csv files the regional brain maps of connectivity to the 
# disease epicenters after stratifying by cognition and symptoms 
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

rm(list=ls())

location <- 'C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/2. MIND/'
source(paste0(location,'regional_brainmap_representation_borders.R'))

epicenter<-c('centiles_MIND_CN',
            'centiles_MIND_FEP_FEP',

'centiles_MIND_FEP_FEP_low_cognition',
'centiles_MIND_FEP_FEP_high_cognition',
'centiles_MIND_FEP_FEP_low_BPRS',
'centiles_MIND_FEP_FEP_high_BPRS',
'centiles_MIND_FEP_FEP_low_SAPS',
'centiles_MIND_FEP_FEP_high_SAPS',
'centiles_MIND_FEP_FEP_low_SANS',
'centiles_MIND_FEP_FEP_high_SANS')

for (i in epicenter){

  epicenters<-read.csv(paste0(location,'Data/connectivity_to_epicenters/epicenters_',i,'.csv'))
  significant_epicenters<-read.csv(paste0(location,'Data/connectivity_to_epicenters/significant_epicenters_',i,'.csv'))
  epicenters<-data.frame(values=t(epicenters),sig = t(!is.na(significant_epicenters)))
  map<-regional_brainmap_representation_borders(epicenters,i,0.42,-0.52,0,'vertical')
  
  print(map)
  
  ggsave(paste0(location,'Data/connectivity_to_epicenters/plots/',i,'.png'), width = 8, height = 5, dpi = 300)
}
