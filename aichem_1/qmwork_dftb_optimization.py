#provide computational engine for rl

import os ; os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
from qmworks import Settings, templates, run, molkit
from noodles import gather
from qmworks.packages.SCM import dftb
from qmworks.components import Distance
from qmworks.packages.PyFragModules import GetFragmentList
from qmworks import plams

# User define Settings
settings = Settings()
settings.specific.dftb.dftb.scc.ndiis = 4
settings.specific.dftb.dftb.scc.Mixing = 0.1
settings.specific.dftb.dftb.scc.iterations = 300
settings.specific.dftb.system.charge = -1
consset = Settings()

ircCoordfile ="/home/x2sun/aichem/aichem/molecule_i.xyz"
fragDefinitions = [{'frag1': [5], 'frag2': [1,2,3,4,6]}]
bondDefinitions = [{'bond1': Distance(0, 4) , 'bond2': Distance(0, 5)}]
bondlengthR1 = 2.993
endbondlengthR1= 2.242
steps = 20
bondlengthR2 = 1.819
energyR1 = -8.7115308
energyR2 = -8.7000666
job_list = []

def jobrun():
   for ircIndex, ircFrags in enumerate(GetFragmentList(ircCoordfile, fragDefinitions)):
       ircTag  = str(ircIndex+1).zfill(4)
       p_mol   = ircFrags['complex']
       bond1   = bondDefinitions[ircIndex]['bond1']
       bond2   = bondDefinitions[ircIndex]['bond2']
       steplength = (bondlengthR1 - endbondlengthR1)/steps


       for i in range(steps):
           name = str(i+1)
           bondlength = bondlengthR1
           bondlength -= steplength*(i+1)
           consset.constraint.update(bond1.get_settings(bondlength))
           p_dftb = dftb(templates.geometry.overlay(settings).overlay(consset), p_mol, job_name="p_DFTB" + name)
           job_list.append(gather(p_dftb))
   # Finalize and draw workflow
   wf = gather(*job_list)
   # Actual execution of the jobs
   results = run(wf, n_processes=20)
   finalresult = []
   for p_result in results:
       for r in  p_result:
           try:
           # Retrieve the molecular coordinates
               mol = r.molecule
               d1  = bond1.get_current_value(mol)
               d2  = bond2.get_current_value(mol)
               ene = r.energy
               bondlength_r1 = bondlengthR1 - d1
               bondlength_r2 = d2  - bondlengthR2
               energy_r1     = ene - energyR1
               energy_r2     = ene - energyR2
           except:
               # Not converge, mission failed.
               bondlength_r1, bondlength_r2, energy_r1, energy_r2    = 0, 0, 0, 0
           finalresult.append([bondlength_r1, bondlength_r2, energy_r1, energy_r2])
   return finalresult

fresult =jobrun()
print (fresult)
