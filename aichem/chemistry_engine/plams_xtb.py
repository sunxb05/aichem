from plams_setting import gfnxtb_go,  list_AMSjobs
import os


working_directory = os.getcwd()
starting_structure = "the_new_molecules"           #name of the foleder with starting structures - to be optimized

ircCoordfile ="/home/x2sun/aichem/aichem/molecule_i.xyz"
fragDefinitions = [{'frag1': [5], 'frag2': [1,2,3,4,6]}]
bondDefinitions = [{'bond1': Distance(0, 4) , 'bond2': Distance(0, 5)}]


class ChemEngine:
   def __init__(self):
      self.observation_space   = ['B1', 'E']
      self.action_space        = ['bond1']

   def run(self, action):

      # Reads molecules from a folder
      molecules = read_molecules(os.path.join(working_directory,starting_structure))
      # define preoptimization engine
      geo_opt= gfnxtb_go(freq=False)
      # Defines list of jobs
      list_go_jobs = list_AMSjobs('opt',geo_opt,molecules)
      # Runs jobs and from results takes optimized coordinates
      opt_res = [optimisation.run() for optimisation in list_go_jobs]
      # Checks if geometry converged
      sp_mol_en = {res.job.molname: res.get_energy()  for res in opt_res}


      return bondlength_r1, bondlength_r2, energy_r1, energy_r2
