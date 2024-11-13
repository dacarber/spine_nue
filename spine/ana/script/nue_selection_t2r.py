"""Analysis module template.
Nue analysis script with truth to reco matching.
This is used to find the efficiency of the nue script and
will save data products into a csv file
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
__all__ = ['nuet2rAna']


class nuet2rAna(AnaBase):
    name = 'nue_t2r'

    def __init__(self, **kwargs):
        """Initialize the analysis script.
        """
        # Initialize the parent class
        super().__init__('interaction', 'both', **kwargs)
        cfg_path = '/sdf/home/d/dcarber/spine_nue/spine/ana/nue_anaconfig.cfg'
        cfg = yaml.safe_load(open(cfg_path, 'r'))
        # Initialize the CSV writer(s) you want
        job_id = os.environ['SLURM_JOB_ID']
        #local_id =  os.environ['SLURM_LOCALID']
        print(os.environ['SLURM_JOB_ID']," ", os.environ['SLURM_JOB_GID'], " ",  os.environ['SLURM_LOCALID'])
        #job_id = 1
        # Initialize the CSV writer(s) you want
        
        self.initialize_writer(f'log_{job_id}')
        
        self.keys['interaction_matches_t2r'] = True

    def process(self, data):
        """Pass data products corresponding to one entry through the analysis.

        Parameters
        ----------
        data : dict
            Dictionary of data products
        """
        

        selection = nue_analysis()
        # Loop over matched interactions (r2t)
        interaction_matches_t2r = data['interaction_matches_t2r']
        for match in interaction_matches_t2r:
            
            # Get match components
            
            reco_inter = match[1]
            true_inter = match[0]
            event_status = None
            if reco_inter == None:
                continue
            if true_inter.nu_id < 0:
                continue
            energy = 0
            for par in true_inter.particles:
                if par.is_primary:
                    energy += par.energy_init
            #df_int_e.append(energy)
            # Containment cut
            #if not true_inter.is_contained : continue
            true_containment = True
            for par in true_inter.particles:
                if par.pid >1 and par.is_contained == False and par.is_primary:
                    #print(par)
                    true_containment = False
                    #continue
            reco_containment = True
            bad_interaction = False
            for par in reco_inter.particles:
                if par.pid >1 and par.is_contained == False and par.is_primary:
                    reco_containment = False
                    bad_interaction = True
                    
            # Primary electron cut
            true_electron = [p for p in true_inter.particles if (p.pid == 1) and (p.is_primary) and (p.energy_init > 70)]
            #if len(true_electron) != 1 : continue

            # Primary photons cut
            true_protons = [t for t in true_inter.particles if (t.pid == 4) and (t.is_primary) and (t.energy_init > 40)]
            
            #Sort protons by highest energy
            true_protons = sorted([tp for tp in true_protons], key=lambda tp : tp.energy_deposit, reverse=True)


            
             # Containment cut
            #if not reco_inter.is_contained : 
            #    bad_interaction = True
            
            # Fiducial cut
            #if not reco_inter.is_fiducial : continue
            
            # Flash cut
            if reco_inter.flash_time < 0 or reco_inter.flash_time > 9.6 : 
                bad_interaction = True
            reco_electrons = [t for t in reco_inter.particles if (t.pid == 1) and (t.is_primary) and (t.calo_ke > 70)]
            if len(reco_electrons) != 1 :
                bad_interaction = True
            reco_electrons = sorted([te for te in reco_electrons], key=lambda te : te.calo_ke, reverse=True)
            matched_electron = None
            
            for reco_e in reco_inter.particles:
                if len(true_electron) <= 0:
                    break
                if len(true_electron[0].match_ids) > 0:
                    if reco_e.id == true_electron[0].match_ids[0]:
                        matched_electron = reco_e
                        break
            
            reco_protons = [p for p in reco_inter.particles if (p.pid == 4) and (p.is_primary) and (p.csda_ke > 40)]
            reco_protons = sorted([rp for rp in reco_protons], key=lambda rp : rp.calo_ke, reverse=True) 
            

            reco_photons = [p for p in reco_inter.particles if (p.pid == 0) and (p.is_primary) and (p.calo_ke > 25)]
            #reco_photons.append([p for p in reco_inter.particles if (p.pid == 1) and (p.is_primary) and (p.calo_ke > 80) and (selection.conversion_dist(p, reco_inter.vertex) > 2)])
            reco_muon = [p for p in reco_inter.particles if (p.pid == 2) and (p.is_primary) and (p.csda_ke > 25)]
            reco_pion = [p for p in reco_inter.particles if (p.pid == 3) and (p.is_primary) and (p.csda_ke > 25)]
            topology = f'{len(reco_photons)}g{len(reco_electrons)}e{len(reco_muon)}m{len(reco_pion)}pi{len(reco_protons)}p'

            true_photons = [p for p in true_inter.particles if (p.pid == 0) and (p.is_primary) and (p.energy_init > 25)]
            true_muon = [p for p in true_inter.particles if (p.pid == 2) and (p.is_primary) and (p.energy_init > 25)]
            true_pion = [p for p in true_inter.particles if (p.pid == 3) and (p.is_primary) and (p.energy_init > 25)]
            true_topology =  f'{len(true_photons)}g{len(true_electron)}e{len(true_muon)}m{len(true_pion)}pi{len(true_protons)}p'

             #Finds the category of interactions i.e. 1e1p, backgrounds
            reco_category = selection.reco_category(reco_inter,topology, bad_interaction,reco_containment)
            true_category = selection.true_category(true_inter,true_topology,true_containment)
            
            # This is our pi0 event
            reco_total_energy = 0
            for particles in reco_inter.particles:
                if particles.pid < 2 and particles.is_primary:
                    reco_total_energy += particles.calo_ke
                elif particles.pid ==2 and particles.is_primary:
                    reco_total_energy += particles.csda_ke
                elif particles.pid ==3 and particles.is_primary:
                    reco_total_energy += particles.csda_ke
                elif particles.pid ==4 and particles.is_primary:
                    reco_total_energy += particles.csda_ke


            # Calculation of conversion dist
            if len(reco_electrons) > 0 :
                reco_conversion_dist = selection.conversion_dist(reco_electrons[0], reco_inter.vertex)


            
            # reco dict corresponding to a CSV row
            reco_nue_dict = {}
            reco_nue_dict['reco_bad_event'] = bad_interaction
            reco_nue_dict['reco_interaction_id'] = match[0].id
            reco_nue_dict['reco_flash_time'] = reco_inter.flash_time
            reco_nue_dict['reco_containment'] = reco_containment
            reco_nue_dict['reco_topology'] = topology
            reco_nue_dict['reco_number_protons'] = len(reco_protons)
            reco_nue_dict['reco_vertex_x'] = reco_inter.vertex[0]
            reco_nue_dict['reco_vertex_y'] = reco_inter.vertex[1]
            reco_nue_dict['reco_vertex_z'] = reco_inter.vertex[2]
            reco_nue_dict['reco_nu_energy'] = reco_total_energy
            reco_nue_dict['reco_event_status'] = event_status
            reco_nue_dict['reco_category'] = reco_category

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
            if len(reco_electrons)>0:
                reco_nue_dict['reco_leading_electron_start_point_x'] = reco_electrons[0].start_point[0]
                reco_nue_dict['reco_leading_electron_start_point_y'] = reco_electrons[0].start_point[1]
                reco_nue_dict['reco_leading_electron_start_point_z'] = reco_electrons[0].start_point[2]
                reco_nue_dict['reco_leading_electron_energy'] = reco_electrons[0].calo_ke
                reco_nue_dict['reco_leading_electron_size'] = reco_electrons[0].size
                reco_nue_dict['reco_leading_electron_pid_score'] = reco_electrons[0].pid_scores[0]
                reco_nue_dict['reco_leading_electron_shape'] = reco_electrons[0].shape
                reco_nue_dict['reco_leading_electron_conversion_dist'] = reco_conversion_dist
            else:
                reco_nue_dict['reco_leading_electron_start_point_x'] = None
                reco_nue_dict['reco_leading_electron_start_point_y'] = None
                reco_nue_dict['reco_leading_electron_start_point_z'] = None
                reco_nue_dict['reco_leading_electron_energy'] = None
                reco_nue_dict['reco_leading_electron_size'] = None
                reco_nue_dict['reco_leading_electron_pid_score'] = None
                reco_nue_dict['reco_leading_electron_shape'] = None
                reco_nue_dict['reco_leading_electron_conversion_dist'] = None
                
            reco_nue_dict['true_topology'] = true_topology
            reco_nue_dict['true_containment'] = true_inter.is_contained
            reco_nue_dict['true_neutrino_id'] = true_inter.nu_id
            reco_nue_dict['true_nu_energy'] = true_inter.energy_init
            reco_nue_dict['true_category'] = true_category
            reco_nue_dict['true_flash_time'] = true_inter.particles[0].parent_t
            reco_nue_dict['true_neutrino_energy'] = energy


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

            if len(true_electron)>0:
                reco_nue_dict['true_leading_electron_energy_deposit'] = true_electron[0].energy_deposit
                reco_nue_dict['true_leading_electron_energy_init'] = true_electron[0].energy_init
                reco_nue_dict['true_electron_start_point_x'] = true_electron[0].start_point[0]
                reco_nue_dict['true_electron_start_point_y'] = true_electron[0].start_point[1]
                reco_nue_dict['true_electron_start_point_z'] = true_electron[0].start_point[2]
                reco_nue_dict['true_electron_size'] = true_electron[0].size
                reco_nue_dict['true_electron_shape'] = true_electron[0].shape
                reco_nue_dict['true_electron_pid'] = true_electron[0].pid

            else:
                reco_nue_dict['true_leading_electron_energy_deposit'] = None
                reco_nue_dict['true_leading_electron_energy_init'] = None
                reco_nue_dict['true_electron_start_point_x'] = None
                reco_nue_dict['true_electron_start_point_y'] = None
                reco_nue_dict['true_electron_start_point_z'] = None
                reco_nue_dict['true_electron_size'] = None
                reco_nue_dict['true_electron_shape'] = None
                reco_nue_dict['true_electron_pid'] = None
            
            if matched_electron is not None:
                reco_nue_dict['reco_match_electron_energy'] = matched_electron.calo_ke
                reco_nue_dict['reco_match_electron_start_point_x'] = matched_electron.start_point[0]
                reco_nue_dict['reco_match_electron_start_point_y'] = matched_electron.start_point[1]
                reco_nue_dict['reco_match_electron_start_point_z'] = matched_electron.start_point[2]
                reco_nue_dict['reco_match_electron_size'] = matched_electron.size
                reco_nue_dict['reco_match_electron_shape'] = matched_electron.shape
                reco_nue_dict['reco_match_electron_pid'] = matched_electron.pid
            else:
                reco_nue_dict['reco_match_electron_energy'] = None
                reco_nue_dict['reco_match_electron_start_point_x'] = None
                reco_nue_dict['reco_match_electron_start_point_y'] = None
                reco_nue_dict['reco_match_electron_start_point_z'] = None
                reco_nue_dict['reco_match_electron_size'] = None
                reco_nue_dict['reco_match_electron_shape'] = None
                reco_nue_dict['reco_match_electron_pid'] = None
                
            reco_nue_dict['true_interaction_id'] = true_inter.id
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
            for keys in reco_nue_dict.keys():
                print(keys," : ",reco_nue_dict[keys])
            job_id = os.environ['SLURM_JOB_ID']
            # Append row to CSV
            self.append(f'log_{job_id}', **reco_nue_dict)

    
    
    def opening_angle(particle_1,particle_2):
        particle_1_dir = particle_1.start_dir / np.linalg.norm(particle_1.start_dir)
        particle_2_dir = particle_2.start_dir / np.linalg.norm(particle_2.start_dir)
        cos_opening_angle = np.dot(particle_1_dir, particle_2_dir)
        opening_angle = np.degrees(np.arccos(cos_opening_angle))
        return opening_angle

    
