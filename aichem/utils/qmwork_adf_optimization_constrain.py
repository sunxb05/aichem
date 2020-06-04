#provide computational engine for rl

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
ircCoordfile ="/home/x2sun/aichem/aichem/utils_I/molecule_rc.xyz"
# fragDefinitions = [{'frag1': [5], 'frag2': [1,2,3,4,6]}]
# consset = Settings()

# bondDefinitions = [{'bond1': Distance(0, 4) , 'bond2': Distance(0, 5)}]
# bondlengthR1 = 3.40
# endbondlengthR1= 1.844
# steps = 20
# bondlengthR2 = 1.844
# energyR1 = -589.67
# energyR2 = -589.67


fragDefinitions = [{'frag1': [6], 'frag2': [1,2,3,4,5]}]
consset = Settings()

bondDefinitions = [{'bond1': Distance(0, 5) , 'bond2': Distance(0, 4)}]
bondlengthR1 = 4.039
endbondlengthR1= 2.434
steps = 20
bondlengthR2 = 1.827
energyR1 = -578.48
energyR2 = -562.13


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
           p_adf = adf(templates.geometry.overlay(settings).overlay(consset), p_mol, job_name="ADF" + name)
           job_list.append(gather(p_adf))
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
               ene = r.energy * 627.51
               bondlength_r1 = bondlengthR1 - d1
               bondlength_r2 = d2  - bondlengthR2
               energy_r1     = ene - energyR1
               energy_r2     = ene - energyR2
           except:
               # Not converge, mission failed.
               bondlength_r1, bondlength_r2, energy_r1, energy_r2    = 0, 0, 0, 0
           finalresult.append([bondlength_r1, bondlength_r2, energy_r1, energy_r2])
   return finalresult

# def writeKey(file, value, pform=r'%7.3f', ljustwidth=16):
#    # write all data into a file.Keep 7 digits and 5 decimals and the width of each entry is 16
#    for val in value:
#       if val is None:
#          file.write(str.ljust('  ---', ljustwidth))
#       else:
#          if type(val) == float:
#             file.write(str.ljust(pform % (val), ljustwidth))
#          else:
#             file.write(str.ljust(str(val), ljustwidth))
#    file.write('\n')


# def WriteTable(tableValues, fileName):
#    energyfile  = open('pyfrag'+fileName+'.txt', "w")
#    for entry in tableValues:
#       sortedEntry = [entry[i] for i in headerlist]
#       writeKey(energyfile, sortedEntry)
#    energyfile.close()


fresult =jobrun()

print (fresult)
