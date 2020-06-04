from scm.plams import *

##########   settting     ###############################
def dftb3_3ob_GO(freq=False):
    settings = Settings()
    settings.input.AMS.Task = 'GeometryOptimization'
    settings.input.DFTB.ResourcesDir = 'DFTB.org/3ob-3-1'
    settings.input.DFTB.Model= 'DFTB3'
    settings.input.DFTB.DispersionCorrection = 'D3-BJ'
    if freq:
        settings.input.AMS.Properties.NormalModes = 'Yes'
    return settings

def dftb3_3ob_SP(charge=0):
    settings = Settings()
    settings.input.AMS.Task = 'SinglePoint'
    settings.input.DFTB.Model = 'DFTB3'
    settings.input.DFTB.ResourcesDir = 'DFTB.org/3ob-3-1'
    settings.input.DFTB.Occupation.Strategy = 'Aufbau'
    if charge != 0:
        settings.input.DFTB.SCC.Iterations = 700
        settings.input.AMS.System.Charge = charge
    return settings
def gfnxtb_go():
    settings = Settings()
    settings.input.ams.Task = 'SinglePoint'
    settings.input.DFTB.Model= 'GFN1-xTB'
    return settings
def adf_sp(basis='DZP', xc_fun=('gga','BLYP'), charge=0, tddft=False, scf_iter=150):
    settings = Settings()
    settings.input.SCF.Iterations = scf_iter
    settings.input.XC[xc_fun[0]] = xc_fun[1]
    settings.input.Basis.Type = basis
    settings.input.Basis.Core = 'None'
    settings.input.Basis.Createoutput = 'None'
    settings.input.NumericalQuality = 'Normal'
    settings.input.Symmetry = 'auto'
    settings.input.Occupations = 'IntegerAufbau'
    settings.input.Totalenergy = ''
    if charge != 0:
        settings.input.Unrestricted = ''
        settings.input.Charge = charge
    return settings
##############properties#######################
@add_to_class(AMSResults)
def get_orbitals(results):
    engine_name = results.engine_names()
    if results.job.ok() and engine_name[0]=='dftb':
        orbitals = Units.convert(results.readrkf('Orbitals', 'Energies(1)', engine_name[0]),'hartree','eV')
        occupations = results.readrkf('Orbitals', 'Occupations(1)', engine_name[0])
    else:
        orbitals, occupations = []
    return orbitals, occupations
###################utility#####################
def list_AMSjobs(job_name,job_sett, molecules):
  list_of_jobs = []
  for name, coords in molecules.items():
    job = AMSJob(name = job_name+name, molecule = coords, settings = job_sett)
    job.molname = name
    list_of_jobs.append(job)
  return list_of_jobs
def list_ADFjobs(job_name,job_sett, molecules):
  list_of_jobs = []
  for name, coords in molecules.items():
    job = ADFJob(name = job_name+name, molecule = coords, settings = job_sett)
    job.molname = name
    list_of_jobs.append(job)
  return list_of_jobs
