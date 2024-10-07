'''
Analysis modules that are shared functions between the 
nue selections, this includes the categorization, 
conversion distance, and counting particles

'''

import numpy as np
import pandas as pd
import yaml, os, sys, re
import glob
import math
from scipy.spatial import distance


class nue_analysis:
    def true_category(self,interaction,topology):
        category = ''
        #topology = interaction.topology
        if interaction.nu_id >= 0:
            particles = self.count_particles(topology)
            if particles[1] == 1 and particles[2] == 0:
                if particles[3] == 0 and particles[0] == 0 and particles[4] == 1 and interaction.is_contained == True:
                    category = '1e1p'
                elif particles[3] == 0 and particles[0] == 0 and particles[4] == 1 and interaction.is_contained == False:
                    category = 'uncontained 1e1p'
                elif particles[3] == 0 and particles[0] == 0 and particles[4] == 0 and interaction.is_contained == True:
                    category = '1e'
                elif particles[3] == 0 and particles[0] == 0 and particles[4] == 0 and interaction.is_contained == False:
                    category = 'uncontained 1e'
                elif particles[3] == 0 and particles[0] == 0 and particles[4] >1 and interaction.is_contained == True:
                    category = '1eNp'
                elif particles[3] == 0 and particles[0] == 0 and particles[4] >1 and interaction.is_contained == False:
                    category = 'uncontained 1eNp'
                elif particles[3] == 1 and particles[0] == 0 and particles[4] ==1 and interaction.is_contained == True:
                    category = '1e1pi1p' 
                elif particles[3] == 1 and particles[0] == 0 and particles[4] ==1 and interaction.is_contained == False:
                    category = 'uncontained 1e1pi1p'
                elif interaction.current_type == 0 and interaction.is_contained == True:
                    category = 'Nue Other'
                elif interaction.current_type == 0 and interaction.is_contained == False:
                    category = 'uncontained Nue Other'
                elif interaction.current_type == 1 and interaction.is_contained == True:
                    category = 'NC'
            elif particles[2] != 0 and interaction.current_type == 0:
                category = 'Numu'
            elif interaction.current_type == 1 and interaction.is_contained == True:
                category = 'NC'
            elif interaction.is_contained == True and abs(interaction.pdg_code) == 12: 
                category = 'Nue Other'
            elif interaction.is_contained == False and abs(interaction.pdg_code) == 12: 
                category = 'uncontained Nue Other'
            elif interaction.is_contained == True and abs(interaction.pdg_code) == 14: 
                category = 'Numu'
            else:
                print(topology)
    
        else:
            category= 'cosmic'
        return category

    def reco_category(self,interaction,topology,event_status):
        #topology = interaction.topology
        catergory = ''
        particles = self.count_particles(topology)
        if topology == None:
            category = None
        elif (particles[1] == 1 and particles[2] == 0 and event_status != True):
            if interaction.is_contained == False:
                category = 'uncontained'
            elif particles[1] == 1 and particles[0] == 0 and particles[3] == 0 and particles[4] == 1 and interaction.is_contained == True:
                category = '1e1p'
            elif particles[1] == 1 and particles[0] == 0 and particles[3] == 0 and particles[4] ==0 and interaction.is_contained == True:
                category = '1e'
            elif particles[1] == 1 and particles[0] == 0 and  particles[3] == 0 and particles[4] >1 and interaction.is_contained == True:
                category = '1eNp'
            elif particles[1] == 1 and particles[0] == 0 and particles[3] == 1 and particles[4] ==1 and interaction.is_contained == True:
                category = '1e1pi1p'
            else: 
                category = 'Nue Other'
        elif interaction.flash_time >9.6 or interaction.flash_time <0:
            category = 'Flash fail'
        elif interaction.is_contained == False:
            category = 'Contain fail'
        elif particles[1] >1:
            category = 'Too many Electron fail'
        elif particles[1] <1 and particles[0] != 0:
            category = 'E to G fail'
        elif particles[1] <1:
            category = 'No Electron fail'
        else:
            category = 'Other fail'

        return category
                
    def count_particles(self,topology):
    #num_e, num_p, num_g, num_pi, num_m = 0,0,0,0,0
        count = [0]*5
        for i in range(len(topology)):
            if topology[i] == 'g':
                count[0] +=int(topology[i-1])
            elif topology[i] == 'e':
                count[1] +=int(topology[i-1])
            elif topology[i] == 'm':
                count[2] +=int(topology[i-1])
            elif topology[i:i+2] == 'pi':
                count[3] +=int(topology[i-1])
            elif topology[i] == 'p' :
                count[4] +=int(topology[i-1])
        return count
        
    def NuMI_angle(self,vertex):
        x = (31512.0380) - vertex[0]
        y = (3364.4912) - vertex[1]
        z = (73363.2532) - vertex[2]
        r = np.sqrt(x**2+y**2+z**2)
        x = x/r
        y = y/r
        z = z/r
        return math.acos(x *vertex[0] + y *vertex[1]+z *vertex[2]);

    def conversion_dist(self,particle,vertex):
        dist = distance.cdist(np.array([vertex]), particle.points)
        #print(dist[0])
        closest_point  = np.min(dist[0])
        return closest_point

    def min_vertex_dist(self,particle,vertex):
        dist = distance.cdist(np.array([vertex]), particle.points)
        #print(dist[0])
        closest_point  = np.min(dist[0])
        return closest_point