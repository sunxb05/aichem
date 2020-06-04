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

ircCoordfile ="/home/x2sun/aichem/aichem/molecule_i_single.xyz"
fragDefinitions = [{'frag1': [5], 'frag2': [1,2,3,4]}]

def jobrun(step_name):
   for ircIndex, ircFrags in enumerate(GetFragmentList(ircCoordfile, fragDefinitions)):
       ircTag  = str(ircIndex+1).zfill(4)
       name    = step_name
       p_mol   = ircFrags['complex']
       # p_dftb = dftb(templates.singlepoint.overlay(settings), p_mol, job_name=name + "_p_DFTB")
       p_dftb = dftb(templates.singlepoint.overlay(settings), p_mol, job_name=name + "_p_DFTB")
       # Add the jobs to the job list
       job_list = []
       job_list.append(gather(p_dftb))
   # Finalize and draw workflow
   wf = gather(*job_list)
   # Actual execution of the jobs
   results = run(wf, n_processes=1)
   for p_result in results:
       for r in  p_result:
           try:
           # Retrieve the molecular coordinates
               mol = r.molecule
               ene = r.energy
           except:
               # Not converge, mission failed.
               ene = 'failed'
   return ene

energy =jobrun("fragment2")
print (energy)
