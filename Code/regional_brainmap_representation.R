################################################################################
# Function to plot regional brain maps from .csv files
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


regional_brainmap_representation <- function(v,text_Grap,sup_lim,inf_lim,midd_p){
  zones <- c( 'bankssts',
              'caudal anterior cingulate',
              'caudal middle frontal',
              'cuneus',
              'entorhinal',
              'frontal pole',
              'fusiform',
              'inferior parietal',
              'inferior temporal',
              'insula',
              'isthmus cingulate',
              'lateral occipital',
              'lateral orbitofrontal',
              'lingual',
              'medial orbitofrontal',
              'middle temporal',
              'paracentral',
              'parahippocampal',
              'pars opercularis',
              'pars orbitalis',
              'pars triangularis',
              'pericalcarine',
              'postcentral',
              'posterior cingulate',
              'precentral',
              'precuneus',
              'rostral anterior cingulate',
              'rostral middle frontal',
              'superior frontal',
              'superior parietal',
              'superior temporal',
              'supramarginal',
              'temporal pole',
              'transverse temporal')
  
  if(length(v) == 34){
    someData <- tibble(
      region = zones, 
      values = v,
      groups = c(rep(text_Grap, length(v)))
      )
      
  }else if (length(v) == 68){
    someData <- tibble(
      label = c(paste0('lh_',gsub(' ','',zones)),paste0('rh_',gsub(' ','',zones))),
      values = v,
      groups = c(rep(text_Grap, length(v)))
    )
  }
  
  # Brain map configuration
  color_palette <- colorRampPalette(c("white", "pink","red")) # scale_fill_gradientn
  someData %>%
    group_by(groups) %>%
    ggplot() + scale_fill_gradientn(
      
      colors = color_palette(100), # scale_fill_gradientn
      
      # low = "blue",
      # mid = "white",
      # high = "red",
      
      # low = "white",
      # mid= 'red',
      # high = "black",
      
      # low = "white",
      # mid = "turquoise4",
      # high = "#005654",

      # low = "white",
      # mid = "#FF97CB",
      # high = "#E60073",
      

      # midpoint = midd_p,
      limits=c(inf_lim,sup_lim),
      space = "Lab",
      na.value = "grey50",
      # guide = guide_colorbar(title = "age (years)"),
      guide = guide_colorbar(title = "",
                             direction = 'horizontal',
                             barwidth = 10,
                             label.theme = element_text(color = "black"),
                             title.theme = element_text(color = "black")),
      # guide = "colourbar",
      aesthetics = "fill"
    ) +
    geom_brain( 
      atlas = dk, 
     position = position_brain(hemi ~ side),show.legend = TRUE,
     aes(fill = values)) +
    facet_wrap(~groups)+theme_brain2(plot.background = "white")
  
}
