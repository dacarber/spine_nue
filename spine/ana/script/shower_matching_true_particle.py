"""Analysis module template.

Use this template as a basis to build your own analysis script. An analysis
script takes the output of the reconstruction and the post-processors and
performs basic selection cuts and store the output to a CSV file.
"""

# Add the imports specific to this module here
import numpy as np
import pandas as pd
import yaml, os, sys, re
from spine.ana.script.nue_analysis_modules import nue_analysis
# Must import the analysis script base class
from spine.ana.base import AnaBase

# Must list the post-processor(s) here to be found by the factory.
# You must also add it to the list of imported modules in the
# `spine.ana.factories`!
__all__ = ['showert2rAna']


class showert2rAna(AnaBase):
    """Script to make numu CC pi0 selection"""
    name = 'shower_t2r'

    def __init__(self, **kwargs):
        """Initialize the analysis script.
        """
        # Initialize the parent class
        super().__init__('interaction', 'both', **kwargs)
        cfg_path = '/sdf/home/d/dcarber/spine_nue/spine/ana/nue_anaconfig.cfg'
        cfg = yaml.safe_load(open(cfg_path, 'r'))
        # Initialize the CSV writer(s) you want
        job_id = os.environ['SLURM_JOB_ID']
        #job_id = 1
        #local_id =  os.environ['SLURM_LOCALID']
        #print(os.environ['SLURM_JOB_ID']," ", os.environ['SLURM_JOB_GID'], " ",  os.environ['SLURM_LOCALID'])
        #job_id = 1
        # Initialize the CSV writer(s) you want
        
        self.initialize_writer(f'log_{job_id}')
        
        #self.keys['particle_matches_t2r']= True
        self.keys['interaction_matches_t2r'] = True

    def process(self, data):
        """Pass data products corresponding to one entry through the analysis.

        Parameters
        ----------
        data : dict
            Dictionary of data products
        """
        


        # Loop over matched interactions (r2t)
        #r2t_matched_particles = data['particle_matches_r2t']
        #for match in r2t_matched_particles:
        selection = nue_analysis()
        reco_nue_dict = {}
        t2r_matched_interactions = data['interaction_matches_t2r']
        for match in t2r_matched_interactions:
            
            # Get match components
            reco_inter = match[1]
            true_inter = match[0]
            if reco_inter == None:
                continue
            #reco_par = reco_inter.particles
            #true_par = true_inter.particles
            #reco_par = match[0]
            #true_par = match[1]
            pi0_1 = -1
            pi0_2 = -1
            for true_par in true_inter.particles:
                if true_par.pdg_code == 22 and true_par.is_primary and true_par.creation_process == 'Decay' and true_par.ancestor_pdg_code == 111:
                    pi0_tag = True
                    if pi0_1 == -1:
                        pi0_1 = true_par.energy_init
                    elif pi0_2 == -1:
                        pi0_2 = true_par.energy_init
            for true_par in true_inter.particles:
                matched_par = None
                for reco_par in reco_inter.particles:
                    if  len(true_par.match_ids)>0:
                        if reco_par.id == true_par.match_ids[0]:
                            matched_par = reco_par
                            break
                if matched_par == None:
                    continue
                if true_par == None:
                    continue
                #if reco_par == None:
                #    continue
                # Containment cut
                #if not true_par.is_contained : continue
    
                if not true_par.is_primary: continue
                
                #Photon Shower
                if true_par.pid > 1 and true_par.pid<4:
                    continue
            
                #if reco_par.pid == 1 and matched_par.pid == 2:
                #    print("True interaction",true_inter)
                #    print("Reco Interaction",reco_inter)
                #    print("Data file",data['file_index'])
                reco_conversion_dist = selection.conversion_dist(reco_par, reco_inter.vertex)
  
            # Fiducial cut
            #if not reco_inter.is_fiducial : continue
          
            # This is our pi0 event
            #reco_pi0_photons = reco_photons[:2]
            
            # reco dict corresponding to a CSV row
            
                reco_nue_dict['reco_id'] = matched_par.id
                reco_nue_dict['reco_contain'] = matched_par.is_contained
                reco_nue_dict['reco_calo_ke'] = matched_par.calo_ke
                reco_nue_dict['reco_match'] = matched_par.match_overlaps[0]
                reco_nue_dict['reco_start_x'] = matched_par.start_point[0]
                reco_nue_dict['reco_start_y'] = matched_par.start_point[1]
                reco_nue_dict['reco_start_z'] = matched_par.start_point[2]
                reco_nue_dict['reco_primary_score'] = matched_par.primary_scores[0]
                reco_nue_dict['reco_pid_score'] = matched_par.pid_scores[0]
                reco_nue_dict['reco_pid'] = matched_par.pid
                reco_nue_dict['reco_start_dir_x'] = matched_par.start_dir[0]
                reco_nue_dict['reco_start_dir_y'] = matched_par.start_dir[1]
                reco_nue_dict['reco_start_dir_z'] = matched_par.start_dir[2]
                reco_nue_dict['reco_size'] = matched_par.size
                reco_nue_dict['reco_shape'] = matched_par.shape
                reco_nue_dict['reco_conversion_dist'] = reco_conversion_dist
    
    
                reco_nue_dict['true_id'] = true_par.id
                reco_nue_dict['true_contain'] = true_par.is_contained
                reco_nue_dict['true_energy_init'] = true_par.energy_init
                reco_nue_dict['true_energy_deposit'] = true_par.energy_deposit
                reco_nue_dict['true_match'] = true_par.match_overlaps[0]
                reco_nue_dict['true_start_x'] = true_par.start_point[0]
                reco_nue_dict['true_start_y'] = true_par.start_point[1]
                reco_nue_dict['true_start_z'] = true_par.start_point[2]
                reco_nue_dict['true_primary'] = true_par.is_primary
                reco_nue_dict['true_pid'] = true_par.pid
                reco_nue_dict['true_start_dir_x'] = true_par.start_dir[0]
                reco_nue_dict['true_start_dir_y'] = true_par.start_dir[1]
                reco_nue_dict['true_start_dir_z'] = true_par.start_dir[2]
                reco_nue_dict['true_size'] = true_par.size
                reco_nue_dict['true_shape'] = true_par.shape
                reco_nue_dict['true_time'] = true_par.t
                reco_nue_dict['reco_vertex_x'] = reco_inter.vertex[0]
                reco_nue_dict['reco_vertex_y'] = reco_inter.vertex[1]
                reco_nue_dict['reco_vertex_z'] = reco_inter.vertex[2]
    
                reco_nue_dict['true_vertex_x'] = true_inter.vertex[0]
                reco_nue_dict['true_vertex_y'] = true_inter.vertex[1]
                reco_nue_dict['true_vertex_z'] = true_inter.vertex[2]
                reco_nue_dict['true_neutrino_id'] = true_inter.nu_id
                reco_nue_dict['pi0_energy_init'] = pi0_1 + pi0_2
           


            ### To-do ##########################################################
            # Add truth info, so we can tell if we are selecting pi0s properly
            # Angle wrt NuMI beam
            # Opening Angle
            # Azmuthal Angle
            # Polar Angle
            # Momentum
            # Change proton energy reconstruction to csda
            # True containment
            # 
            # Energy bias
            # Maybe primary id
            # 
            # 
            # Do everything for truth to reco matching
            #
            ####################################################################
                job_id = os.environ['SLURM_JOB_ID']
                #job_id = 1
                # Append row to CSV
                self.append(f'log_{job_id}', **reco_nue_dict)

    