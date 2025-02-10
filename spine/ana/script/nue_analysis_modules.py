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
from scipy.spatial.distance import cdist
from scipy.spatial import distance


class nue_analysis:
    def true_category(self,interaction,topology,containment):
        category = ''
        #topology = interaction.topology
        if interaction.nu_id >= 0:
            particles = self.count_particles(topology)
            if particles[1] == 1 and particles[2] == 0:
                if particles[3] == 0 and particles[0] == 0 and particles[4] == 1 and containment == True:
                    category = '1e1p'
                elif particles[3] == 0 and particles[0] == 0 and particles[4] == 1 and containment == False:
                    category = 'uncontained 1e1p'
                elif particles[3] == 0 and particles[0] == 0 and particles[4] == 0 and containment == True:
                    category = '1e'
                elif particles[3] == 0 and particles[0] == 0 and particles[4] == 0 and containment == False:
                    category = 'uncontained 1e'
                elif particles[3] == 0 and particles[0] == 0 and particles[4] >1 and containment == True:
                    category = '1eNp'
                elif particles[3] == 0 and particles[0] == 0 and particles[4] >1 and containment == False:
                    category = 'uncontained 1eNp'
                elif particles[3] == 1 and particles[0] == 0 and particles[4] ==1 and containment == True:
                    category = '1e1pi1p' 
                elif particles[3] == 1 and particles[0] == 0 and particles[4] ==1 and containment == False:
                    category = 'uncontained 1e1pi1p'
                elif interaction.current_type == 0 and containment == True:
                    category = 'Nue Other'
                elif interaction.current_type == 0 and containment == False:
                    category = 'uncontained Nue Other'
                elif interaction.current_type == 1 and containment == True:
                    category = 'NC'
            elif particles[2] != 0 and interaction.current_type == 0:
                category = 'Numu'
            elif interaction.current_type == 1 and containment == True:
                category = 'NC'
            elif containment == True and abs(interaction.pdg_code) == 12: 
                category = 'Nue Other'
            elif containment == False and abs(interaction.pdg_code) == 12: 
                category = 'uncontained Nue Other'
            elif containment == True and abs(interaction.pdg_code) == 14: 
                category = 'Numu'
            #else:
            #    print(topology)
    
        else:
            category= 'cosmic'
        return category

    def reco_category(self,interaction,topology,event_status,containment):
        #topology = interaction.topology
        catergory = ''
        particles = self.count_particles(topology)
        if topology == None:
            category = None
        elif (particles[1] == 1 and particles[2] == 0 and event_status != True):
            if containment == False:
                category = 'uncontained'
            elif particles[1] == 1 and particles[0] == 0 and particles[3] == 0 and particles[4] == 1 and containment == True:
                category = '1e1p'
            elif particles[1] == 1 and particles[0] == 0 and particles[3] == 0 and particles[4] ==0 and containment == True:
                category = '1e'
            elif particles[1] == 1 and particles[0] == 0 and  particles[3] == 0 and particles[4] >1 and containment == True:
                category = '1eNp'
            #elif particles[1] == 1 and particles[0] == 0 and particles[3] >= 1 and particles[4] ==0 and containment == True:
            #    category = '1eNpi'
            elif particles[1] == 1 and particles[0] == 0 and particles[3] == 1 and particles[4] >0 and containment == True:
                category = '1e1piNp'
            #elif particles[1] == 1 and particles[0] > 0 and particles[3] >= 1 and particles[4] >0 and containment == True:
            #    category = 'Ng1eNpiNp'
            else: 
                category = 'Nue Other'
        elif interaction.flash_time >9.6 or interaction.flash_time <0:
            category = 'Flash fail'
        elif containment == False:
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
        x = (31512.0380) + vertex[0]
        y = (3364.4912) + vertex[1]
        z = (73363.2532) + vertex[2]
        r = np.sqrt(x**2+y**2+z**2)
        x = x/r
        y = y/r
        z = z/r
        return math.acos(x *vertex[0] + y *vertex[1]+z *vertex[2]);
    def opening_angle(self,leading_e, leading_p):
        return math.acos(leading_e.start_dir[0]*leading_p.start_dir[0]+leading_e.start_dir[1]*leading_p.start_dir[1]+leading_e.start_dir[2]*leading_p.start_dir[2]);

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
    def leading_particle_index(self,interaction,pid,true_particle):
        leading_ke = 0
        index =0
        for i in range(len(interaction.particles)):
            p = interaction.particles[i];
            energy = p.calo_ke
            if true_particle:
                energy = p.energy_init
            if(p.pid == pid and energy > leading_ke and p.is_primary):
            
                leading_ke = energy;
                index = i;
                
            
            return index;
    def electron_transverse_momentum(self,interaction,true_int):

        dir_vector = [0,0,0]
        dir_vector[0] = interaction.vertex[0] + (31512.0380)
        dir_vector[1] = interaction.vertex[1] + (3364.4912)
        dir_vector[2] = interaction.vertex[2] + (73363.2532)
        r = np.sqrt(pow(dir_vector[0], 2)+pow(dir_vector[1], 2)+pow(dir_vector[2], 2))
        dir_vector[0] = dir_vector[0]/r
        dir_vector[1] = dir_vector[1]/r
        dir_vector[2] = dir_vector[2]/r                                                                                                                                                                                             
        beamdir = [dir_vector[0], dir_vector[1], dir_vector[2]] #NuMI  
        
        i = self.leading_particle_index(interaction, 1, true_int)
        particle = interaction.particles[i];

        # pT = p - pL                               
        #    = p-(p dot beamdir) * beamdir     
        
        p = [0,0,0]
        pL = [0,0,0]
        pT = [0,0,0]

        p[0] = np.sqrt(particle.calo_ke**2 + 2*particle.calo_ke*0.511998)*particle.start_dir[0]
        p[1] = np.sqrt(particle.calo_ke**2 + 2*particle.calo_ke*0.511998)*particle.start_dir[1]
        p[2] = np.sqrt(particle.calo_ke**2 + 2*particle.calo_ke*0.511998)*particle.start_dir[2]

        pL = np.dot(p,beamdir) * np.array(beamdir)
        pT = p - pL
        
        return pT

    def proton_transverse_momentum(self,interaction):
        dir_vector = [0,0,0]
        dir_vector[0] = interaction.vertex[0] + (31512.0380)
        dir_vector[1] = interaction.vertex[1] + (3364.4912)
        dir_vector[2] = interaction.vertex[2] + (73363.2532)
        r = np.sqrt(pow(dir_vector[0], 2)+pow(dir_vector[1], 2)+pow(dir_vector[2], 2))
        dir_vector[0] = dir_vector[0]/r
        dir_vector[1] = dir_vector[1]/r
        dir_vector[2] = dir_vector[2]/r                                                                                                                                                                                             
        #beamdir = [dir_vector[0], dir_vector[1], dir_vector[2]] #NuMI                                                                                                                                                      0.39431672, 0.04210058, 0.91800973                                           
        beamdir = [0.39431672, 0.04210058, 0.91800973 ]
        
        pT0,pT1,pT2 = 0,0,0                    
        
        for particle in interaction.particles:
            if particle.pid == 4 and particle.is_primary:
        # pT = p - pL 
        #    = p-(p dot beamdir) * beamdir  
                
                p = [0,0,0]
                pL = [0,0,0]
                pT = [0,0,0]
                
                p[0] = particle.momentum[0]
                p[1] = particle.momentum[1]
                p[2] = particle.momentum[2]
        
                pL = np.dot(p,beamdir) * np.array(beamdir)
                pT = p - pL
                pT0 += pT[0]
                pT1 += pT[1]
                pT2 += pT[2]
        ppT = [pT0,pT1,pT2]
        
        return ppT
        
    def delta_pT(self,interaction,true_int):
        pLT = self.electron_transverse_momentum(interaction, true_int)
        ppT = self.proton_transverse_momentum(interaction)
        delta_p = np.array([pLT[0]+ppT[0],pLT[1]+ppT[1],pLT[2]+ppT[2]])
        mag = np.linalg.norm(delta_p)
        return mag
    def delta_alphaT(self,interaction,true_int):
        pLT = np.array(self.electron_transverse_momentum(interaction, true_int))
        ppT = np.array(self.proton_transverse_momentum(interaction))
        delta_p = np.array([pLT[0]+ppT[0],pLT[1]+ppT[1],pLT[2]+ppT[2]])
        print(np.dot(-pLT,delta_p)/(np.linalg.norm(pLT)*np.linalg.norm(delta_p)))
        if -1 > np.dot(-pLT,delta_p)/(np.linalg.norm(pLT)*np.linalg.norm(delta_p)):
            return math.acos(-1)
        elif np.dot(-pLT,delta_p)/(np.linalg.norm(pLT)*np.linalg.norm(delta_p)) > 1:
            return math.acos(1)
        
        alphaT = math.acos(np.dot(-pLT,delta_p)/(np.linalg.norm(pLT)*np.linalg.norm(delta_p)))
        return alphaT
    def delta_phiT(self,interaction,true_int):
        pLT = np.array(self.electron_transverse_momentum(interaction, true_int))
        ppT = np.array(self.proton_transverse_momentum(interaction))
        if -1 > np.dot(-pLT,ppT)/(np.linalg.norm(pLT)*np.linalg.norm(ppT)):
            return math.acos(-1)
        elif np.dot(-pLT,ppT)/(np.linalg.norm(pLT)*np.linalg.norm(ppT)) > 1:
            return math.acos(1)
        phiT = math.acos(np.dot(-pLT,ppT)/(np.linalg.norm(pLT)*np.linalg.norm(ppT)))
        return phiT
    
    def polar_angle(self,particle):
        return math.acos(particle.start_dir[2])
        
    def azimuthal_angle(self,particle):
        if particle.start_dir[1] > 0:
            return math.acos(particle.start_dir[0] / np.sqrt(pow(particle.start_dir[0], 2) + pow(particle.start_dir[1], 2)))
        else:
            return -math.acos(particle.start_dir[0] / np.sqrt(pow(particle.start_dir[0], 2) + pow(particle.start_dir[1], 2)))
    def dedx(self, particle,vertex,start_dist = 10):
        vertex = np.array(vertex)
        points = np.array(particle.points)
        dist = cdist(vertex,points)
        mask = dist < start_dist
        dedx = sum(particle.depositions[mask])/start_dist
        