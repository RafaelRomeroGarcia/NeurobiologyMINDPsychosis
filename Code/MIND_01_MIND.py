#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 10:21:18 2023

@author: natalia

Script to compute MIND networks

Copyright (C) 2024 University of Seville

Written by Natalia García San Martín (ngarcia1@us.es)

This file is part of Neurobiology MIND Psychosis toolkit.

Neurobiology MIND Mapping toolkit is free software: 
you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation, 
either version 3 of the License, or (at your option) any later version.

Neurobiology MIND Mapping toolkit is distributed in the hope that 
it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Neurobiology MIND Mapping toolkit. If not, see 
<https://www.gnu.org/licenses/>.

"""
import os
import sys
import pandas as pd
sys.path.insert(1, '/home/rafael/Downloads/MIND-master/')
from MIND import compute_MIND

main_path = '/home/rafael/data/PAFIP/SURFER/' # freesurfer directory
saving_path = '/home/rafael/Downloads/MIND_networks_PAFIP/' # saving path

folders = [subject for subject in os.listdir(main_path) if os.path.isdir(os.path.join(main_path, subject))]
folders.remove('fsaverage_HCP')
# folders = folders[142:]

MIND = pd.DataFrame()
errors = []
for subject in folders:
    ## Specify path to surfer directory. This is the path to the Freesurfer 
    ## folder containing all standard output directories (e.g. surf, mri, label)
    path_to_surf_dir =  main_path + subject
    
    ## Specify features to include in MIND calculation. The following 
    ## abbreviations specifies the ?h.thickness, ?h.curv, ?h.volume, ?h.sulc, 
    ## and ?h.area Freesurfer surface features.
    features = ['CT','MC','Vol','SD','SA'] 
    
    ## Select which parcellation to use. This has been tested on 
    ## Desikan Killiany (DK), HCP-Glasser, DK-308 and DK-318 parcellations.
    parcellation = 'aparc' 
    
    try: 
        ## Returns a dataframe of regions X regions containing the final MIND network.
        MIND = compute_MIND(path_to_surf_dir, features, parcellation) 
        
        ## Save MIND network   
        MIND.to_csv(saving_path+'aparc/'+subject+'.csv')
    except Exception:
        errors.append(subject)
        continue
    
