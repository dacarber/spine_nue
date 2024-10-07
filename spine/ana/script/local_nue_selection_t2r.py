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
from spine.ana.script.nue_analysis_modules import nue_analysis
#sys.path.append('/sdf/home/d/dcarber/spine_nue/')
# Must import the analysis script base class
from spine.ana.base import AnaBase

# Must list the post-processor(s) here to be found by the factory.
# You must also add it to the list of imported modules in the
# `spine.ana.factories`!
__all__ = ['local_nuet2rAna']


class local_nuet2rAna(AnaBase):
    """Script to make numu CC pi0 selection"""
    name = 'local_nue_t2r'

    def __init__(self, **kwargs):
        """Initialize the analysis script.
        """
        # Initialize the parent class
        super().__init__('interaction', 'both', **kwargs)
        #job_id = os.environ['SLURM_JOB_ID']
        #local_id =  os.environ['SLURM_LOCALID']
        #print(os.environ['SLURM_JOB_ID']," ", os.environ['SLURM_JOB_GID'], " ",  os.environ['SLURM_LOCALID'])
        #job_id = 1
        # Initialize the CSV writer(s) you want
        
        self.initialize_writer(f'log_part_t2r_v2')
        #self.append = True
        self.keys['interaction_matches_t2r'] = True
        #self.keys['interaction_matches_t2r'] = True

    def process(self, data):
        """Pass data products corresponding to one entry through the analysis.

        Parameters
        ----------
        data : dict
            Dictionary of data products
        """
        


        #interaction_matches_t2r = data['interaction_matches_t2r']
        #for match in interaction_matches_r2t:
        #    
        #    # Get match components
        #    reco_inter = match[0]
        #    true_inter = match[1]

        # Loop over matched interactions (r2t)
        #csv_files = glob.glob('nue*.{}'.format('csv'))
        #print(csv_files)
        selection = nue_analysis()
        interaction_matches_t2r = data['interaction_matches_t2r']
        for match in interaction_matches_t2r:
            
            # Get match components
            reco_inter = match[1]
            true_inter = match[0]
            event_status = None
            if reco_inter == None:
                continue
            
            # Containment cut
            if not true_inter.is_contained : continue
            
            # Primary electron cut
            true_electron = [p for p in true_inter.particles if (p.pid == 1) and (p.is_primary) and (p.energy_deposit > 80)]
            if len(true_electron) != 1 : continue

            # Primary photons cut
            print('1')
            true_protons = [t for t in true_inter.particles if (t.pid == 4) and (t.is_primary) and (t.energy_deposit > 50)]
            
            #Sort protons by highest energy
            
            true_protons = sorted([tp for tp in true_protons], key=lambda tp : tp.energy_deposit, reverse=True)


            bad_interaction = False
             # Containment cut
            if not reco_inter.is_contained : 
                bad_interaction = True
            
            # Fiducial cut
            #if not reco_inter.is_fiducial : continue
            
            # Flash cut
            if reco_inter.flash_time < 0 or reco_inter.flash_time > 9.6 : 
                bad_interaction = True
            reco_electrons = [t for t in reco_inter.particles if (t.pid == 1) and (t.is_primary) and (t.calo_ke > 80)]
            if len(reco_electrons) != 1 :
                bad_interaction = True
            reco_electrons = sorted([te for te in reco_electrons], key=lambda te : te.calo_ke, reverse=True)
            matched_electron = None
            for reco_e in reco_inter.particles: 
                if len(true_electron[0].match_ids) > 0:
                    if reco_e.id == true_electron[0].match_ids[0]:
                        matched_electron = reco_e
                        break
            
            reco_protons = [p for p in reco_inter.particles if (p.pid == 4) and (p.is_primary) and (p.csda_ke > 50)]
            reco_protons = sorted([rp for rp in reco_protons], key=lambda rp : rp.calo_ke, reverse=True) 
            

            reco_photons = [p for p in reco_inter.particles if (p.pid == 0) and (p.is_primary) and (p.calo_ke > 25)]
            reco_muon = [p for p in reco_inter.particles if (p.pid == 2) and (p.is_primary) and (p.csda_ke > 25)]
            reco_pion = [p for p in reco_inter.particles if (p.pid == 3) and (p.is_primary)]
            topology = f'{len(reco_photons)}g{len(reco_electrons)}e{len(reco_muon)}m{len(reco_pion)}pi{len(reco_protons)}p'

            true_photons = [p for p in true_inter.particles if (p.pid == 0) and (p.is_primary) and (p.energy_deposit > 25)]
            true_muon = [p for p in true_inter.particles if (p.pid == 2) and (p.is_primary) and (p.energy_deposit > 25)]
            true_pion = [p for p in true_inter.particles if (p.pid == 3) and (p.is_primary)]
            true_topology =  f'{len(true_photons)}g{len(true_electron)}e{len(true_muon)}m{len(true_pion)}pi{len(true_protons)}p'

            reco_category = selection.reco_category(reco_inter,topology, event_status)
            true_category = selection.true_category(true_inter,true_topology)

            if (true_category == 'Nue Other' or true_category == '1e' or true_category == '1e1pi1p' or true_category == '1e1p' or true_category == '1eNp'):
                print("Nue Inclusive",reco_category,":",data['file_index']," Topology: ", topology)
            
            # This is our pi0 event
            reco_total_energy = 0
            for particles in reco_inter.particles:
                if particles.pid < 2 and particles.is_primary:
                    reco_total_energy += particles.calo_ke
                elif particles.pid >=2 and particles.is_primary:
                    reco_total_energy += particles.csda_ke
            #print(reco_total_energy)
            if len(reco_electrons) > 0:
                reco_conversion_dist = selection.conversion_dist(reco_electrons[0], reco_inter.vertex)
            #reco_pi0_photons = reco_photons[:2]
            
            # reco dict corresponding to a CSV row
            reco_nue_dict = {}
            reco_nue_dict['File_id'] = data['file_index']
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

            if len(reco_protons)>0:
                reco_nue_dict['reco_leading_proton_start_point_x'] = reco_protons[0].start_point[0]
                reco_nue_dict['reco_leading_proton_start_point_y'] = reco_protons[0].start_point[1]
                reco_nue_dict['reco_leading_proton_start_point_z'] = reco_protons[0].start_point[2]
                reco_nue_dict['reco_leading_proton_energy'] = reco_protons[0].csda_ke
            else: 
                reco_nue_dict['reco_leading_proton_start_point_x'] = None
                reco_nue_dict['reco_leading_proton_start_point_y'] = None
                reco_nue_dict['reco_leading_proton_start_point_z'] = None
                reco_nue_dict['reco_leading_proton_energy'] = None
            if len(reco_electrons)>0:
                reco_nue_dict['reco_electron_start_point_x'] = reco_electrons[0].start_point[0]
                reco_nue_dict['reco_electron_start_point_y'] = reco_electrons[0].start_point[1]
                reco_nue_dict['reco_electron_start_point_z'] = reco_electrons[0].start_point[2]
                reco_nue_dict['reco_electron_energy'] = reco_electrons[0].calo_ke
                reco_nue_dict['reco_leading_electron_conversion_dist'] = reco_conversion_dist
            else:
                reco_nue_dict['reco_electron_start_point_x'] = None
                reco_nue_dict['reco_electron_start_point_y'] = None
                reco_nue_dict['reco_electron_start_point_z'] = None
                reco_nue_dict['reco_electron_energy'] = None
                reco_nue_dict['reco_leading_electron_conversion_dist'] = None
            
            reco_nue_dict['true_topology'] = true_topology
            reco_nue_dict['true_category'] = true_category
            reco_nue_dict['true_flash_time'] = true_inter.particles[0].parent_t
            reco_nue_dict['true_containment'] = true_inter.is_contained
            reco_nue_dict['true_neutrino_id'] = true_inter.nu_id
            reco_nue_dict['true_nu_energy'] = true_inter.energy_init

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
            else:
                reco_nue_dict['true_leading_electron_energy_deposit'] = None
                reco_nue_dict['true_leading_electron_energy_init'] = None
                reco_nue_dict['true_electron_start_point_x'] = None
                reco_nue_dict['true_electron_start_point_y'] = None
                reco_nue_dict['true_electron_start_point_z'] = None
                
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
            #job_id = 1
            # Append row to CSV
            self.append(f'log_part_t2r_v2', **reco_nue_dict)
            
    
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
        elif particles[1] == 1 and particles[2] == 0:
            if interaction.is_contained == False:
                category = 'uncontained'
            if particles[1] == 1 and particles[0] == 0 and particles[3] == 0 and particles[4] == 1 and interaction.is_contained == True:
                category = '1e1p'
            elif particles[1] == 1 and particles[0] == 0 and particles[3] == 0 and particles[4] ==0 and interaction.is_contained == True:
                category = '1e'
            elif particles[1] == 1 and particles[0] == 0 and  particles[3] == 0 and particles[4] >1 and interaction.is_contained == True:
                category = '1eNp'
            elif particles[1] == 1 and particles[0] == 0 and particles[3] == 1 and particles[4] ==1 and interaction.is_contained == True:
                category = '1e1pi1p'
            else: 
                category = 'Nue Other'
        else: 
            category = 'Not selected'
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
                    
