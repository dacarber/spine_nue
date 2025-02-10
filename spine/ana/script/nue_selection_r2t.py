"""Analysis module template.

Use this template as a basis to build your own analysis script. An analysis
script takes the output of the reconstruction and the post-processors and
performs basic selection cuts and store the output to a CSV file.
"""

# Add the imports specific to this module here
import numpy as np
import pandas as pd
import yaml, os, sys, re
import glob
import math
from scipy.spatial import distance
#sys.path.append('/sdf/home/d/dcarber/spine_nue/')
# Must import the analysis script base class
from spine.ana.base import AnaBase
from spine.ana.script.nue_analysis_modules import nue_analysis

# Must list the post-processor(s) here to be found by the factory.
# You must also add it to the list of imported modules in the
# `spine.ana.factories`!
__all__ = ['nuer2tAna']


class nuer2tAna(AnaBase):
    """Script to make numu CC pi0 selection"""
    name = 'nue_r2t'

    def __init__(self, **kwargs):
        """Initialize the analysis script.
        """
        # Initialize the parent class
        super().__init__('interaction', 'both', **kwargs)
        job_id = os.environ['SLURM_JOB_ID']
        local_id =  os.environ['SLURM_LOCALID']
        print(os.environ['SLURM_JOB_ID']," ", os.environ['SLURM_JOB_GID'], " ",  os.environ['SLURM_LOCALID'])
        #job_id = 1
        # Initialize the CSV writer(s) you want
        
        self.initialize_writer(f'log_{job_id}_{local_id}')
        #self.append = True
        self.keys['interaction_matches_r2t'] = True
        #self.keys['interaction_matches_t2r'] = True

    def process(self, data):
        """Pass data products corresponding to one entry through the analysis.

        Parameters
        ----------
        data : dict
            Dictionary of data products
        """
        

        selection = nue_analysis()

        interaction_matches_r2t = data['interaction_matches_r2t']
        for match in interaction_matches_r2t:
            
            # Get match components
            reco_inter = match[0]
            true_inter = match[1]
            event_status = ''
            if true_inter == None:
                continue
            
            # Containment cut
            reco_containment = True
            for par in reco_inter.particles:
                if par.pid >1 and par.is_contained == False and par.is_primary:
                    event_status+='c'
                    reco_containment = False
                    continue
            true_containment = True
            for par in true_inter.particles:
                if par.pid >1 and par.is_contained == False and par.is_primary:
                    true_containment = False
            
            #if not reco_inter.is_contained : 
            #    event_status+='c'
            #    continue
            
            # Fiducial cut
            #if not reco_inter.is_fiducial : continue
            
            # Flash cut
            if reco_inter.flash_time < 0 or reco_inter.flash_time > 9.6 : 
                event_status+='f'
                continue
            # Grabs electrons and requires that there is only 1 primary electron shower
            reco_electron = [p for p in reco_inter.particles if (p.pid == 1) and (p.is_primary) and (p.calo_ke > 70)]
            print(reco_inter.particles[0].__dir__())
            if len(reco_electron) != 1 : 
                reco_electron = sorted([te for te in reco_electron], key=lambda te : te.calo_ke, reverse=True)
                event_status+='e'
                continue
            
            true_electrons = [t for t in true_inter.particles if (t.pid == 1) and (t.is_primary) and (t.energy_init > 70)]
            true_electrons = sorted([te for te in true_electrons], key=lambda te : te.energy_deposit, reverse=True)
            matched_electron = None
            #print(reco_electron)
            if len(reco_electron) > 0:
                for true_e in true_inter.particles:
                    if (len(reco_electron[0].match_ids) <1):
                        break
                    if true_e.id == reco_electron[0].match_ids[0]:
                        matched_electron = true_e
                        break

            reco_protons = [p for p in reco_inter.particles if (p.pid == 4) and (p.is_primary) and (p.csda_ke > 40)]
            true_protons = [t for t in true_inter.particles if (t.pid == 4) and (t.is_primary) and (t.energy_init > 40)]
            
            #Sort protons by highest energy
            reco_protons = sorted([rp for rp in reco_protons], key=lambda rp : rp.csda_ke, reverse=True) 
            true_protons = sorted([tp for tp in true_protons], key=lambda tp : tp.energy_deposit, reverse=True)

            reco_photons = [p for p in reco_inter.particles if (p.pid == 0) and (p.is_primary) and (p.calo_ke > 25)]
            reco_muon = [p for p in reco_inter.particles if (p.pid == 2) and (p.is_primary) and (p.csda_ke > 25)]
            reco_pion = [p for p in reco_inter.particles if (p.pid == 3) and (p.is_primary) and (p.csda_ke > 25)]
            topology = f'{len(reco_photons)}g{len(reco_electron)}e{len(reco_muon)}m{len(reco_pion)}pi{len(reco_protons)}p'


            #Grabs true photons and tags if they come from a Pi0
            pi0_tag = False
            true_photons = [p for p in true_inter.particles if (p.pid == 0) and (p.is_primary)and (p.energy_init > 25)]
            if len(true_photons) >0:
                for tph in true_photons:
                    if tph.pdg_code == 22 and tph.is_primary and tph.creation_process == 'Decay' and tph.ancestor_pdg_code == 111:
                        pi0_tag = True
            true_muon = [p for p in true_inter.particles if (p.pid == 2) and (p.is_primary) and (p.energy_init > 25)]
            true_pion = [p for p in true_inter.particles if (p.pid == 3) and (p.is_primary) and (p.energy_init> 25)]
            true_topology =  f'{len(true_photons)}g{len(true_electrons)}e{len(true_muon)}m{len(true_pion)}pi{len(true_protons)}p'
            print(true_topology, " : ", true_inter.topology)

            #Finds the category of interactions i.e. 1e1p, backgrounds
            reco_category = selection.reco_category(reco_inter,topology, event_status,reco_containment)
            true_category = selection.true_category(true_inter,true_topology,true_containment)
            
            # Reco energy calculation
            reco_total_energy = 0
            for particles in reco_inter.particles:
                if particles.pid < 2 and particles.is_primary:
                    reco_total_energy += particles.calo_ke/77.0777/0.81*81.3955
                elif particles.pid ==2 and particles.is_primary:
                    reco_total_energy += particles.csda_ke
                elif particles.pid ==3 and particles.is_primary:
                    reco_total_energy += particles.csda_ke
                elif particles.pid ==4 and particles.is_primary:
                    reco_total_energy += particles.csda_ke

            # Calculation of conversion dist
            if len(reco_electron) > 0:
                reco_conversion_dist = selection.conversion_dist(reco_electron[0], reco_inter.vertex)

            #Finds angles of particles
            reco_open_angle = -1
            if len(reco_electron) > 0:
                if len(reco_protons) > 0:
                    reco_open_angle = selection.opening_angle(reco_electron[0],reco_protons[0])
                reco_e_polar = selection.polar_angle(reco_electron[0])
                reco_e_azim = selection.azimuthal_angle(reco_electron[0])
            true_open_angle = -1
            if len(true_electrons) > 0:
                if len(true_protons) > 0:
                    true_open_angle = selection.opening_angle(true_electrons[0],true_protons[0])
                true_e_polar = selection.polar_angle(true_electrons[0])
                true_e_azim = selection.azimuthal_angle(true_electrons[0])
            # reco dict corresponding to a CSV row
            reco_nue_dict = {}
            reco_nue_dict['reco_event_status'] = event_status
            reco_nue_dict['reco_category'] = reco_category
            reco_nue_dict['reco_flash_time'] = reco_inter.flash_time
            reco_nue_dict['reco_containment'] = reco_inter.is_contained
            reco_nue_dict['reco_topology'] = topology
            reco_nue_dict['reco_number_protons'] = len(reco_protons)
            reco_nue_dict['reco_interaction_id'] = match[0].id
            reco_nue_dict['reco_vertex_x'] = reco_inter.vertex[0]
            reco_nue_dict['reco_vertex_y'] = reco_inter.vertex[1]
            reco_nue_dict['reco_vertex_z'] = reco_inter.vertex[2]
            reco_nue_dict['reco_nu_energy'] = reco_total_energy
            reco_nue_dict['reco_dpT'] = selection.delta_pT(reco_inter,true_int=False)
            reco_nue_dict['reco_dalphaT'] = selection.delta_alphaT(reco_inter,true_int=False)
            reco_nue_dict['reco_dphiT'] = selection.delta_phiT(reco_inter,true_int=False)
            reco_nue_dict['reco_e_pT'] = np.linalg.norm(np.array(selection.electron_transverse_momentum(reco_inter,true_int=False)))
            reco_nue_dict['reco_p_pT'] = np.linalg.norm(np.array(selection.proton_transverse_momentum(reco_inter)))
            reco_nue_dict['reco_opening_anglee'] = reco_open_angle
            


            if len(reco_protons)>0:
                reco_nue_dict['reco_leading_proton_start_point_x'] = reco_protons[0].start_point[0]
                reco_nue_dict['reco_leading_proton_start_point_y'] = reco_protons[0].start_point[1]
                reco_nue_dict['reco_leading_proton_start_point_z'] = reco_protons[0].start_point[2]
                reco_nue_dict['reco_leading_proton_energy'] = reco_protons[0].csda_ke
                reco_nue_dict['reco_leading_proton_pid_score'] = reco_protons[0].pid_scores[4]
                reco_nue_dict['reco_leading_proton_length'] = reco_protons[0].length
                
            else: 
                reco_nue_dict['reco_leading_proton_start_point_x'] = None
                reco_nue_dict['reco_leading_proton_start_point_y'] = None
                reco_nue_dict['reco_leading_proton_start_point_z'] = None
                reco_nue_dict['reco_leading_proton_energy'] = None
                reco_nue_dict['reco_leading_proton_pid_score'] = None
                reco_nue_dict['reco_leading_proton_length'] = None
            if len(reco_electron)>0:
                reco_nue_dict['reco_electron_start_point_x'] = reco_electron[0].start_dir[0]
                reco_nue_dict['reco_electron_start_point_y'] = reco_electron[0].start_dir[1]
                reco_nue_dict['reco_electron_start_point_z'] = reco_electron[0].start_dir[2]
                reco_nue_dict['reco_electron_energy'] = reco_electron[0].calo_ke
                reco_nue_dict['reco_electron_size'] = reco_electron[0].size
                reco_nue_dict['reco_electron_pid_score'] = reco_electron[0].pid_scores[0]
                reco_nue_dict['reco_electron_shape'] = reco_electron[0].shape
                reco_nue_dict['reco_electron_conversion_dist'] = reco_conversion_dist
                reco_nue_dict['reco_electron_polar'] = reco_e_polar
                reco_nue_dict['reco_electron_azimuthal'] = reco_e_azim
            else:
                reco_nue_dict['reco_electron_start_point_x'] = None
                reco_nue_dict['reco_electron_start_point_y'] = None
                reco_nue_dict['reco_electron_start_point_z'] = None
                reco_nue_dict['reco_electron_energy'] = None
                reco_nue_dict['reco_electron_size'] = None
                reco_nue_dict['reco_electron_pid_score'] = None
                reco_nue_dict['reco_electron_shape'] = None
                reco_nue_dict['reco_electron_conversion_dist'] = None
                reco_nue_dict['reco_electron_polar'] = None
                reco_nue_dict['reco_electron_azimuthal'] = None
            
            reco_nue_dict['true_topology'] = true_topology
            reco_nue_dict['true_category'] = true_category
            reco_nue_dict['true_flash_time'] = true_inter.particles[0].parent_t
            reco_nue_dict['true_containment'] = true_containment
            reco_nue_dict['true_fiducial'] = true_inter.is_fiducial
            reco_nue_dict['true_neutrino_id'] = true_inter.nu_id
            reco_nue_dict['true_nu_energy'] = true_inter.energy_init
            reco_nue_dict['true_dpT'] = selection.delta_pT(true_inter,true_int=True)
            reco_nue_dict['true_dalphaT'] = selection.delta_alphaT(true_inter,true_int=True)
            reco_nue_dict['true_dphiT'] = selection.delta_phiT(true_inter,true_int=True)
            reco_nue_dict['true_e_pT'] = np.linalg.norm(np.array(selection.electron_transverse_momentum(true_inter,true_int=True)))
            reco_nue_dict['true_p_pT'] = np.linalg.norm(np.array(selection.proton_transverse_momentum(true_inter)))
            reco_nue_dict['true_opening_angle'] = true_open_angle

            if len(true_protons)>0:
                reco_nue_dict['true_leading_proton_energy_deposit'] = true_protons[0].energy_deposit
                reco_nue_dict['true_leading_proton_energy_init'] = true_protons[0].energy_init
                reco_nue_dict['true_leading_proton_start_point_x'] = true_protons[0].start_point[0]
                reco_nue_dict['true_leading_proton_start_point_y'] = true_protons[0].start_point[1]
                reco_nue_dict['true_leading_proton_start_point_z'] = true_protons[0].start_point[2]
            else:
                reco_nue_dict['true_leading_proton_energy_deposit'] = None
                reco_nue_dict['true_leading_proton_energy_init'] = None
                reco_nue_dict['true_leading_proton_start_point_x'] = None
                reco_nue_dict['true_leading_proton_start_point_y'] = None
                reco_nue_dict['true_leading_proton_start_point_z'] = None

            if len(true_electrons)>0:
                reco_nue_dict['true_leading_electron_energy_deposit'] = true_electrons[0].energy_deposit
                reco_nue_dict['true_leading_electron_energy_init'] = true_electrons[0].energy_init
                reco_nue_dict['true_electron_start_point_x'] = true_electrons[0].start_point[0]
                reco_nue_dict['true_electron_start_point_y'] = true_electrons[0].start_point[1]
                reco_nue_dict['true_electron_start_point_z'] = true_electrons[0].start_point[2]
                reco_nue_dict['true_electron_size'] = true_electrons[0].size
                reco_nue_dict['true_electron_shape'] = true_electrons[0].shape
                reco_nue_dict['true_electron_pid'] = true_electrons[0].pid
                reco_nue_dict['true_electron_polar'] = true_e_polar
                reco_nue_dict['true_electron_azimuthal'] = true_e_azim

                
            else:
                reco_nue_dict['true_leading_electron_energy_deposit'] = None
                reco_nue_dict['true_leading_electron_energy_init'] = None
                reco_nue_dict['true_electron_start_point_x'] = None
                reco_nue_dict['true_electron_start_point_y'] = None
                reco_nue_dict['true_electron_start_point_z'] = None
                reco_nue_dict['true_electron_size'] = None
                reco_nue_dict['true_electron_shape'] = None
                reco_nue_dict['true_electron_pid'] = None
                reco_nue_dict['true_electron_polar'] = None
                reco_nue_dict['true_electron_azimuthal'] = None
            if matched_electron is not None:
                reco_nue_dict['true_match_electron_energy_deposit'] = matched_electron.energy_deposit
                reco_nue_dict['true_match_leading_electron_energy_init'] = matched_electron.energy_init
                reco_nue_dict['true_match_electron_start_point_x'] = matched_electron.start_point[0]
                reco_nue_dict['true_match_electron_start_point_y'] = matched_electron.start_point[1]
                reco_nue_dict['true_match_electron_start_point_z'] = matched_electron.start_point[2]
                reco_nue_dict['true_match_electron_size'] = matched_electron.size
                reco_nue_dict['true_match_electron_shape'] = matched_electron.shape
                reco_nue_dict['true_match_electron_pid'] = matched_electron.pid
            else:
                reco_nue_dict['true_match_electron_energy_deposit'] = None
                reco_nue_dict['true_match_leading_electron_energy_init'] = None
                reco_nue_dict['true_match_electron_start_point_x'] = None
                reco_nue_dict['true_match_electron_start_point_y'] = None
                reco_nue_dict['true_match_electron_start_point_z'] = None
                reco_nue_dict['true_match_electron_size'] = None
                reco_nue_dict['true_match_electron_shape'] = None
                reco_nue_dict['true_match_electron_pid'] = None
            reco_nue_dict['true_interaction_id'] = None
            reco_nue_dict['true_pi0'] = pi0_tag
            reco_nue_dict['true_current_type'] = true_inter.current_type
            reco_nue_dict['true_vertex_x'] = true_inter.vertex[0]
            reco_nue_dict['true_vertex_y'] = true_inter.vertex[1]
            reco_nue_dict['true_vertex_z'] = true_inter.vertex[2]


            ### To-do ##########################################################
            # Add truth info, so we can tell if we are selecting pi0s properly
            # Angle wrt NuMI beam
            # Opening Angle
            # Azmuthal Angle
            # Polar Angle
            # Momentum
            # 
            # 
            # 
            # Energy bias
            # Maybe primary id
            # 
            # 
            # Do everything for truth to reco matching
            #
            ####################################################################
            job_id = os.environ['SLURM_JOB_ID']
            local_id =  os.environ['SLURM_LOCALID']
            #job_id = 1
            # Append row to CSV
            self.append(f'log_{job_id}_{local_id}', **reco_nue_dict)
            
    
    
        
            
