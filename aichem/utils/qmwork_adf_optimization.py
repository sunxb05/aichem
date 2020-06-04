import os ; os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
from qmworks import Settings, templates, run, molkit
from noodles import gather
from qmworks.packages.SCM import adf
from qmworks.components import Distance
from qmworks.packages.PyFragModules import GetFragmentList
from qmworks import plams

# User define Settings
settings = Settings()
settings.specific.adf.xc.GGA = "OLYP"
settings.specific.adf.basis.type = "TZP"
settings.specific.adf.CHARGE = -1
ircCoordfile ="/home/x2sun/aichem/aichem/utils/molecule_rc.xyz"
fragDefinitions = [{'frag1': [5], 'frag2': [1,2,3,4,6]}]

def jobrun(step_name):
   for ircIndex, ircFrags in enumerate(GetFragmentList(ircCoordfile, fragDefinitions)):
       ircTag  = str(ircIndex+1).zfill(4)
       name    = step_name
       p_mol   = ircFrags['complex']
       # p_dftb = dftb(templates.singlepoint.overlay(settings), p_mol, job_name=name + "_p_DFTB")
       p_adf = adf(templates.geometry.overlay(settings), p_mol, job_name=name + "_adf")
       # Add the jobs to the job list
       job_list = []
       job_list.append(gather(p_adf))
   # Finalize and draw workflow
   wf = gather(*job_list)
   # Actual execution of the jobs
   results = run(wf, n_processes=24)
   for p_result in results:
       for r in  p_result:
           try:
           # Retrieve the molecular coordinates
               mol = r.molecule
               ene = r.energy * 627.51
           except:
               # Not converge, mission failed.
               ene = 'failed'
   return ene

energy =jobrun("rc")
print (energy)
