from plams import *
import re
import sys, math, os
import numpy as np
import subprocess

# path = os.getcwd()

def GetOutputTable(data):

   outputTable = {}
   for key, val in list(data.items()):
      if type(val) == list and len(val) != 1:
         for i in range(len(val)):
            outputTable[key+'_'+str(i+1)] = val[i]
      elif type(val) == list and len(val) == 1:
         outputTable[key] = val[0]
      else:
         outputTable[key] = val
   return outputTable

def PatternDetermine(mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted):
   if bool(mixList_3):
      return mixList_3, mixList_3_adjusted
   elif bool(mixList_2):
      return mixList_2, mixList_2_adjusted
   else:
      return mixList_1, mixList_1_adjusted

def GetHOMO(coefList):

   swpitem = coefList[0]
   for item in coefList:
      if abs(item['coef']) > abs(swpitem['coef']):
         swpitem = item
   return swpitem

def PyFragCycle(complexOrbital, t21File, mt21File=None, inputKeys=None, mixList=None, Mark_mix_1=False):

   complexResult    = KFFile(t21File)
   pyfragResult     = PyFragResult(complexResult, inputKeys, mt21file=mt21File)

   HOMO_orbitals    = ['HOMO', 'HOMO-1', 'HOMO-2', 'HOMO-3', 'HOMO-4', 'HOMO-5', 'HOMO-6', 'HOMO-7', 'HOMO-8', 'HOMO-9']
   LUMO_orbitals    = ['LUMO', 'LUMO+1', 'LUMO+2', 'LUMO+3', 'LUMO+4', 'LUMO+5', 'LUMO+6', 'LUMO+7', 'LUMO+8', 'LUMO+9']
   comp_unoccupoedMO= LUMO_orbitals


   # make sure number of the HOMOs
   frag1_elecnum, frag2_elecnum = pyfragResult.GetElectronNumber()
   frag1_orbnum      =  int(frag1_elecnum/2)
   frag2_orbnum      =  int(frag2_elecnum/2)
   comp_orbnum       =  frag1_orbnum + frag2_orbnum
   if frag1_orbnum <= 10:
      frag1_orbitals = HOMO_orbitals[0:frag1_orbnum] + LUMO_orbitals
   else:
      frag1_orbitals  = HOMO_orbitals + LUMO_orbitals

   if frag2_orbnum <= 10:
      frag2_orbitals = HOMO_orbitals[0:frag2_orbnum] + LUMO_orbitals
   else:
      frag2_orbitals = HOMO_orbitals + LUMO_orbitals

   if comp_orbnum  <= 10:
      comp_occupiedMO= HOMO_orbitals[0:comp_orbnum]
      comp_orbitals  = comp_occupiedMO + comp_unoccupoedMO
   else:
      comp_occupiedMO= HOMO_orbitals[0:10]
      comp_orbitals  = comp_occupiedMO + comp_unoccupoedMO

   if mixList is not None:
      for mixlist in mixList:
         # remove previous fragment orbitals
         for key, val in list(mixlist['comp_orb_1'].items()):
            if val['frag_type'] == 'frag1':
               frag1_orbitals.remove(val['frag_orb'])
            elif val['frag_type'] == 'frag2':
               frag2_orbitals.remove(val['frag_orb'])

      # for comp_orb, comp_val in list(mixlist.items()):
      #    for frag_orb, frag_val in list(mixlist['comp_orb_1'].items()):
      #       if frag_val['frag_type'] == 'frag2':
      #          coefi = comp_val[frag_orb]['frag_coef']




         # remove previous complex second orbital
         if len(mixlist) == 1:
            if pyfragResult.GetFrontIndex(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])['holu']   == 'HOMO':
               comp_occupiedMO.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])
            else:
               comp_unoccupoedMO.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])

         # remove previous complex second orbital
         if len(mixlist) == 2:

            if pyfragResult.GetFrontIndex(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])['holu']   == 'HOMO':
               comp_occupiedMO.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])
            else:
               comp_unoccupoedMO.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])

            if pyfragResult.GetFrontIndex(mixlist['comp_orb_2']['fragorb_2']['comp_orb'])['holu']   == 'HOMO':
               comp_occupiedMO.remove(mixlist['comp_orb_2']['fragorb_2']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_2']['fragorb_2']['comp_orb'])
            else:
               comp_unoccupoedMO.remove(mixlist['comp_orb_2']['fragorb_2']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_2']['fragorb_2']['comp_orb'])

         # remove previous complex third orbital
         elif len(mixlist) == 3:

            if pyfragResult.GetFrontIndex(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])['holu']   == 'HOMO':
               comp_occupiedMO.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])
            else:
               comp_unoccupoedMO.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_1']['fragorb_1']['comp_orb'])
            if pyfragResult.GetFrontIndex(mixlist['comp_orb_2']['fragorb_2']['comp_orb'])['holu'] == 'HOMO':
               comp_occupiedMO.remove(mixlist['comp_orb_2']['fragorb_2']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_2']['fragorb_2']['comp_orb'])
            else:
               comp_unoccupoedMO.remove(mixlist['comp_orb_2']['fragorb_2']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_2']['fragorb_2']['comp_orb'])

            if pyfragResult.GetFrontIndex(mixlist['comp_orb_3']['fragorb_3']['comp_orb'])['holu'] == 'HOMO':
               comp_occupiedMO.remove(mixlist['comp_orb_3']['fragorb_3']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_3']['fragorb_3']['comp_orb'])
            else:
               comp_unoccupoedMO.remove(mixlist['comp_orb_3']['fragorb_3']['comp_orb'])
               comp_orbitals.remove(mixlist['comp_orb_3']['fragorb_3']['comp_orb'])


   try:
      comp_orbitals.remove(complexOrbital)
      if complexOrbital == 'LUMO':
         comp_unoccupoedMO.remove('LUMO')
      else:
         comp_occupiedMO.remove(complexOrbital)
   except:
      pass


   # test = pyfragResult.ReadSFO({'type': "HOMO-1",   'frag': 'frag2',   'complex': 'HOMO'})
   # # # test = pyfragResult.ReadSFO({'type': "HOMO",   'frag': 'frag2',   'complex': 'HOMO'})
   # print ('test',test)
   # sys.exit()


   # read fragment orbitals coefficient to the complex orbital
   coefList    = []
   for fragorb in frag1_orbitals:
      orbitalCoeffcient = pyfragResult.ReadSFO({'type': fragorb, 'frag': 'frag1', 'complex': complexOrbital})
      coefList.append({'type': fragorb, 'frag': 'frag1', 'coef': orbitalCoeffcient})


   for fragorb in frag2_orbitals:
      orbitalCoeffcient = pyfragResult.ReadSFO({'type': fragorb, 'frag': 'frag2', 'complex': complexOrbital})
      coefList.append({'type': fragorb, 'frag': 'frag2', 'coef': orbitalCoeffcient})


   # choose highest orbital in 2 fragments
   fragHOMO = GetHOMO(coefList)

   # choose highest orbital in another fragment
   coefList_noHOMO = coefList
   coefList_noHOMO.remove(fragHOMO)

   fragList_2  = [item for item in coefList_noHOMO if item['frag'] != fragHOMO['frag']]
   fragHOMO_1  = GetHOMO(fragList_2)

   # choose third highest orbital
   coefList_no2HOMO = coefList_noHOMO
   coefList_no2HOMO.remove(fragHOMO_1)
   fragHOMO_2 = GetHOMO(coefList_no2HOMO)

   # choose two highest complex orbital
   mixList_1   = {}
   mixList_2   = {}
   mixList_3   = {}
   mixList_1_adjusted   = {}
   mixList_2_adjusted   = {}
   mixList_3_adjusted   = {}
   compLUMO    = ''
   compHOMO_1  = ''


   for comporb in comp_orbitals:
      fragHOMOCoef   = pyfragResult.ReadSFO({'type': fragHOMO['type'],   'frag': fragHOMO['frag'],   'complex': comporb})
      fragHOMO_1Coef = pyfragResult.ReadSFO({'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag'], 'complex': comporb})
      if abs(fragHOMOCoef) > 0.01 and abs(fragHOMO_1Coef) > 0.01 and (abs(fragHOMOCoef) + abs(fragHOMO_1Coef)) > 0.3:

         mixList_2 = {
         "comp_orb_1": {
                        'fragorb_1': {
                                    'frag_orb'   : fragHOMO['type'],
                                    'frag_type'  : fragHOMO['frag'],
                                    'comp_orb'   : complexOrbital,
                                    'frag_coef'  : fragHOMO['coef'],
                                    'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                    'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                    'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                    'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                           pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                    },

                        'fragorb_2': {
                                    'frag_orb'   : fragHOMO_1['type'],
                                    'frag_type'  : fragHOMO_1['frag'],
                                    'comp_orb'   : complexOrbital,
                                    'frag_coef'  : fragHOMO_1['coef'],
                                    'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                    'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                    'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                    'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                           pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                    }
                        },

         "comp_orb_2": {
                        'fragorb_1': {
                                    'frag_orb' : fragHOMO['type'],
                                    'frag_type'  : fragHOMO['frag'],
                                    'comp_orb'   : comporb,
                                    'frag_coef'  : fragHOMOCoef,
                                    'comp_energy': pyfragResult.ReadComporbEnergy({'type': comporb}),
                                    'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                    'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                    'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                           pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                    },

                        'fragorb_2': {
                                    'frag_orb' : fragHOMO_1['type'],
                                    'frag_type'  : fragHOMO_1['frag'],
                                    'comp_orb'   : comporb,
                                    'frag_coef'  : fragHOMO_1Coef,
                                    'comp_energy': pyfragResult.ReadComporbEnergy({'type': comporb}),
                                    'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                    'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                    'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                           pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                    }
                        }

         }

         mixList_2_adjusted = {
         "comp_orb_1": {
                        'fragorb_1': {
                                    'frag_orb'   : fragHOMO['type'],
                                    'frag_type'  : fragHOMO['frag'],
                                    'comp_orb'   : complexOrbital,
                                    'frag_coef'  : fragHOMO['coef'],
                                    'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                    'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                    'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO['type'],   'frag': fragHOMO['frag'],   'complex': complexOrbital}),
                                    'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                           pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                    },

                        'fragorb_2': {
                                    'frag_orb'   : fragHOMO_1['type'],
                                    'frag_type'  : fragHOMO_1['frag'],
                                    'comp_orb'   : complexOrbital,
                                    'frag_coef'  : fragHOMO_1['coef'],
                                    'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                    'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                    'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO_1['type'],   'frag': fragHOMO_1['frag'],   'complex': complexOrbital}),
                                    'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                           pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                    }
                        },

         "comp_orb_2": {
                        'fragorb_1': {
                                    'frag_orb' : fragHOMO['type'],
                                    'frag_type'  : fragHOMO['frag'],
                                    'comp_orb'   : comporb,
                                    'frag_coef'  : fragHOMOCoef,
                                    'comp_energy': pyfragResult.ReadComporbEnergy({'type': comporb}),
                                    'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                    'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO['type'],   'frag': fragHOMO['frag'],   'complex': comporb}),
                                    'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                           pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                    },

                        'fragorb_2': {
                                    'frag_orb' : fragHOMO_1['type'],
                                    'frag_type'  : fragHOMO_1['frag'],
                                    'comp_orb'   : comporb,
                                    'frag_coef'  : fragHOMO_1Coef,
                                    'comp_energy': pyfragResult.ReadComporbEnergy({'type': comporb}),
                                    'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                    'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO_1['type'],   'frag': fragHOMO_1['frag'],   'complex': comporb}),
                                    'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                           pyfragResult.GetOrbitalIndex(
                                                                              {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                    }
                        }

         }

         break

   # Mark_mix_1 = False

   # if mixList_2 == {} and complexOrbital == "HOMO-9":
   #    complexOrbital = comp_occupiedMO[0]
   #    Mark_mix_1 = True

   if complexOrbital == "HOMO" or complexOrbital == "LUMO" and mixList_2 == {}:
      Mark_mix_1 = True

   if Mark_mix_1:

      mixList_1 = {
      "comp_orb_1": {
                     'fragorb_1': {
                                 'frag_orb'   : fragHOMO['type'],
                                 'frag_type'  : fragHOMO['frag'],
                                 'comp_orb'   : complexOrbital,
                                 'frag_coef'  : fragHOMO['coef'],
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                 }

                     }

      }


      mixList_1_adjusted = {
      "comp_orb_1": {
                     'fragorb_1': {
                                 'frag_orb'   : fragHOMO['type'],
                                 'frag_type'  : fragHOMO['frag'],
                                 'comp_orb'   : complexOrbital,
                                 'frag_coef'  : fragHOMO['coef'],
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO['type'],   'frag': fragHOMO['frag'],   'complex': complexOrbital}),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                 }
                     }
      }






   # Find lowest LUMO has the 3 frag MOs
   for comporb in comp_unoccupoedMO:
      fragHOMOCoef   = pyfragResult.ReadSFO({'type': fragHOMO['type'],   'frag': fragHOMO['frag'],   'complex': comporb})
      fragHOMO_1Coef = pyfragResult.ReadSFO({'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag'], 'complex': comporb})
      fragHOMO_2Coef = pyfragResult.ReadSFO({'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag'], 'complex': comporb})
      if abs(fragHOMOCoef) > 0.01 and abs(fragHOMO_1Coef) > 0.01 and abs(fragHOMO_2Coef) > 0.01 and (abs(fragHOMOCoef) + abs(fragHOMO_1Coef) + abs(fragHOMO_2Coef)) > 0.3:

         compLUMO                  = comporb
         fragHOMOCoef_compLUMO     = fragHOMOCoef
         fragHOMO_1Coef_compLUMO   = fragHOMO_1Coef
         fragHOMO_2Coef_compLUMO   = fragHOMO_2Coef

         break

   for comporb in comp_occupiedMO:
      fragHOMOCoef   = pyfragResult.ReadSFO({'type': fragHOMO['type'],   'frag': fragHOMO['frag'],   'complex': comporb})
      fragHOMO_1Coef = pyfragResult.ReadSFO({'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag'], 'complex': comporb})
      fragHOMO_2Coef = pyfragResult.ReadSFO({'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag'], 'complex': comporb})
      if abs(fragHOMOCoef) > 0.01 and abs(fragHOMO_1Coef) > 0.01 and abs(fragHOMO_2Coef) > 0.01 and (abs(fragHOMOCoef) + abs(fragHOMO_1Coef) + abs(fragHOMO_2Coef)) > 0.3:

         compHOMO_1                  = comporb
         fragHOMOCoef_compHOMO_1     = fragHOMOCoef
         fragHOMO_1Coef_compHOMO_1   = fragHOMO_1Coef
         fragHOMO_2Coef_ccompHOMO_1  = fragHOMO_2Coef


         break


   if compLUMO and compHOMO_1:

      mixList_3 ={
      "comp_orb_1": {
                     'fragorb_1': {
                                 'frag_orb'   : fragHOMO['type'],
                                 'frag_type'  : fragHOMO['frag'],
                                 'comp_orb'   : complexOrbital,
                                 'frag_coef'  : fragHOMO['coef'],
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                 },

                     'fragorb_2': {
                                 'frag_orb'   : fragHOMO_1['type'],
                                 'frag_type'  : fragHOMO_1['frag'],
                                 'comp_orb'   : complexOrbital,
                                 'frag_coef'  : fragHOMO_1['coef'],
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 },

                     'fragorb_3': {
                                 'frag_orb'   : fragHOMO_2['type'],
                                 'frag_type'  : fragHOMO_2['frag'],
                                 'comp_orb'   : complexOrbital,
                                 'frag_coef'  : fragHOMO_2['coef'],
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']})),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 }


               },

      "comp_orb_2": {
                     'fragorb_1': {
                                 'frag_orb'   : fragHOMO['type'],
                                 'frag_type'  : fragHOMO['frag'],
                                 'comp_orb'   : compLUMO,
                                 'frag_coef'  : fragHOMOCoef_compLUMO,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compLUMO}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                 },

                     'fragorb_2': {
                                 'frag_orb'   : fragHOMO_1['type'],
                                 'frag_type'  : fragHOMO_1['frag'],
                                 'comp_orb'   : compLUMO,
                                 'frag_coef'  : fragHOMO_1Coef_compLUMO,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compLUMO}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 },

                     'fragorb_3': {
                                 'frag_orb'   : fragHOMO_2['type'],
                                 'frag_type'  : fragHOMO_2['frag'],
                                 'comp_orb'   : compLUMO,
                                 'frag_coef'  : fragHOMO_2Coef_compLUMO,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compLUMO}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']})),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 }

                     },

      "comp_orb_3": {
                     'fragorb_1': {
                                 'frag_orb'   : fragHOMO['type'],
                                 'frag_type'  : fragHOMO['frag'],
                                 'comp_orb'   : compHOMO_1,
                                 'frag_coef'  : fragHOMOCoef_compHOMO_1,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compHOMO_1}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                 },

                     'fragorb_2': {
                                 'frag_orb'   : fragHOMO_1['type'],
                                 'frag_type'  : fragHOMO_1['frag'],
                                 'comp_orb'   : compHOMO_1,
                                 'frag_coef'  : fragHOMO_1Coef_compHOMO_1,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compHOMO_1}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 },

                     'fragorb_3': {
                                 'frag_orb'   : fragHOMO_2['type'],
                                 'frag_type'  : fragHOMO_2['frag'],
                                 'comp_orb'   : compHOMO_1,
                                 'frag_coef'  : fragHOMO_2Coef_ccompHOMO_1,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compHOMO_1}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']})),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 }
                     }
      }

      mixList_3_adjusted ={
      "comp_orb_1": {
                     'fragorb_1': {
                                 'frag_orb'   : fragHOMO['type'],
                                 'frag_type'  : fragHOMO['frag'],
                                 'comp_orb'   : complexOrbital,
                                 'frag_coef'  : fragHOMO['coef'],
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO['type'],   'frag': fragHOMO['frag'],   'complex': complexOrbital}),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                 },

                     'fragorb_2': {
                                 'frag_orb'   : fragHOMO_1['type'],
                                 'frag_type'  : fragHOMO_1['frag'],
                                 'comp_orb'   : complexOrbital,
                                 'frag_coef'  : fragHOMO_1['coef'],
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO_1['type'],   'frag': fragHOMO_1['frag'],   'complex': complexOrbital}),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 },

                     'fragorb_3': {
                                 'frag_orb'   : fragHOMO_2['type'],
                                 'frag_type'  : fragHOMO_2['frag'],
                                 'comp_orb'   : complexOrbital,
                                 'frag_coef'  : fragHOMO_2['coef'],
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': complexOrbital}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO_2['type'],   'frag': fragHOMO_2['frag'],   'complex': complexOrbital}),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 }


               },

      "comp_orb_2": {
                     'fragorb_1': {
                                 'frag_orb'   : fragHOMO['type'],
                                 'frag_type'  : fragHOMO['frag'],
                                 'comp_orb'   : compLUMO,
                                 'frag_coef'  : fragHOMOCoef_compLUMO,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compLUMO}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO['type'],   'frag': fragHOMO['frag'],   'complex': compLUMO}),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                 },

                     'fragorb_2': {
                                 'frag_orb'   : fragHOMO_1['type'],
                                 'frag_type'  : fragHOMO_1['frag'],
                                 'comp_orb'   : compLUMO,
                                 'frag_coef'  : fragHOMO_1Coef_compLUMO,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compLUMO}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO_1['type'],   'frag': fragHOMO_1['frag'],   'complex': compLUMO}),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 },

                     'fragorb_3': {
                                 'frag_orb'   : fragHOMO_2['type'],
                                 'frag_type'  : fragHOMO_2['frag'],
                                 'comp_orb'   : compLUMO,
                                 'frag_coef'  : fragHOMO_2Coef_compLUMO,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compLUMO}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO_2['type'],   'frag': fragHOMO_2['frag'],   'complex': compLUMO}),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 }

                     },

      "comp_orb_3": {
                     'fragorb_1': {
                                 'frag_orb'   : fragHOMO['type'],
                                 'frag_type'  : fragHOMO['frag'],
                                 'comp_orb'   : compHOMO_1,
                                 'frag_coef'  : fragHOMOCoef_compHOMO_1,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compHOMO_1}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO['type'],   'frag': fragHOMO['frag'],   'complex': compHOMO_1}),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}))
                                 },

                     'fragorb_2': {
                                 'frag_orb'   : fragHOMO_1['type'],
                                 'frag_type'  : fragHOMO_1['frag'],
                                 'comp_orb'   : compHOMO_1,
                                 'frag_coef'  : fragHOMO_1Coef_compHOMO_1,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compHOMO_1}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO_1['type'],   'frag': fragHOMO_1['frag'],   'complex': compHOMO_1}),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_1['type'], 'frag': fragHOMO_1['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 },

                     'fragorb_3': {
                                 'frag_orb'   : fragHOMO_2['type'],
                                 'frag_type'  : fragHOMO_2['frag'],
                                 'comp_orb'   : compHOMO_1,
                                 'frag_coef'  : fragHOMO_2Coef_ccompHOMO_1,
                                 'comp_energy': pyfragResult.ReadComporbEnergy({'type': compHOMO_1}),
                                 'population' : pyfragResult.ReadPopulation(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']})),
                                 'frag_energy': pyfragResult.ReadFragorbEnergy_adjusted({'type': fragHOMO_2['type'],   'frag': fragHOMO_2['frag'],   'complex': compHOMO_1}),
                                 'overlap'    : pyfragResult.ReadOverlap(pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO['type'], 'frag': fragHOMO['frag']}),
                                                                        pyfragResult.GetOrbitalIndex(
                                                                           {'type': fragHOMO_2['type'], 'frag': fragHOMO_2['frag']}))
                                 }
                     }
      }

   return mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted


def PyFragDriver(t21File, inputKeys, mt21File=None):
   '''
   'overlap': [({'type': 'LUMO+2', 'frag': 'frag1'}, {'type': 'LUMO+2', 'frag': 'frag2'}),({'type': 'LUMO', 'frag': 'frag1'}, {'type': 'HOMO', 'frag': 'frag2'})]

   'orbitalcoeffcient': [{'type': 'LUMO+2', 'frag': 'frag1', 'complex': 'LUMO+2'}, {'type': 'LUMO+2', 'frag': 'frag2', 'complex': 'LUMO+2'}]

   'energy': 'energy': [{'type': 'LUMO+2'}, {'type': 'LUMO+1'}]

   'orbitalenergy': [{'type': 'LUMO+2', 'frag': 'frag1'}, {'type': 'LUMO+1', 'frag': 'frag1'}]
   '''

   resultsTable     = []
   outputData       = {}

   complexResult    = KFFile(t21File)
   pyfragResult     = PyFragResult(complexResult, inputKeys, mt21file=mt21File)

   HOMO_orbitals     = ['HOMO', 'HOMO-1', 'HOMO-2', 'HOMO-3', 'HOMO-4', 'HOMO-5', 'HOMO-6', 'HOMO-7', 'HOMO-8', 'HOMO-9']
   frag1_elecnum, frag2_elecnum = pyfragResult.GetElectronNumber()
   frag1_orbnum      =  int(frag1_elecnum/2)
   frag2_orbnum      =  int(frag2_elecnum/2)
   comp_orbnum       =  frag1_orbnum + frag2_orbnum

   if comp_orbnum <= 10:
      comp_occupiedMO= HOMO_orbitals[1:comp_orbnum]
   else:
      comp_occupiedMO= HOMO_orbitals[1:10]


   mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted= PyFragCycle('HOMO', t21File, mt21File, inputKeys)
   firstmixing, firstmixing_adjusted   = PatternDetermine(mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted)

   if len(firstmixing) != 1 and firstmixing['comp_orb_2']['fragorb_2']['comp_orb'] == 'LUMO':

      Mark = True
      HOMO_low_list     = []
      while Mark:
         # check if there is symmetry
         totalSymmetry    = complexResult.read('Symmetry', 'nsym')
         if totalSymmetry == 1:
            #Read 2 highest occupied MO
            HOMO_energy = pyfragResult.ReadComporbEnergy({'type': 'HOMO'})
            for comporb in comp_occupiedMO:
               comp_occupiedMO.remove(comporb)
               swpOrb_energy   = pyfragResult.ReadComporbEnergy({'type': comporb})
               if (HOMO_energy - swpOrb_energy) > 0.01:
                  HOMO_1 = comporb
                  HOMO_low_list.append(HOMO_1)
                  break
         else:
            HOMO_sym          = pyfragResult.GetCompOrbSym({'type': 'HOMO'})
            HOMO_energy       = pyfragResult.ReadComporbEnergy({'type': 'HOMO'})
            for comporb in comp_occupiedMO:
               comp_occupiedMO.remove(comporb)
               swpOrb_sym     = pyfragResult.GetCompOrbSym({'type': comporb})
               swpOrb_energy   = pyfragResult.ReadComporbEnergy({'type': comporb})

               if swpOrb_sym != HOMO_sym and (HOMO_energy - swpOrb_energy) > 0.01:
                  HOMO_1 = comporb
                  HOMO_low_list.append(HOMO_1)
                  break

         mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted = PyFragCycle(HOMO_1, t21File, mt21File, inputKeys, [firstmixing])
         thirdmixing, thirdmixing_adjusted = PatternDetermine(mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted)
         # comp_occupiedMO.remove(HOMO_1)
         if bool(thirdmixing):
            Mark = False
         elif len(comp_occupiedMO) == 0:
            mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted = PyFragCycle(HOMO_low_list[0], t21File, mt21File, inputKeys, [firstmixing, secondmixing],  Mark_mix_1=True)
            thirdmixing, thirdmixing_adjusted = PatternDetermine(mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted)
            Mark = False

      return [firstmixing, thirdmixing],[firstmixing_adjusted,thirdmixing_adjusted]


   else:


      mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted= PyFragCycle('LUMO', t21File, mt21File, inputKeys, [firstmixing])
      secondmixing, secondmixing_adjusted = PatternDetermine(mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted)

      HOMO_low_list     = []
      Mark = True

      while Mark:
         # check if there is symmetry
         totalSymmetry    = complexResult.read('Symmetry', 'nsym')
         if totalSymmetry == 1:
            #Read 2 highest occupied MO
            HOMO_energy = pyfragResult.ReadComporbEnergy({'type': 'HOMO'})
            for comporb in comp_occupiedMO:
               comp_occupiedMO.remove(comporb)
               swpOrb_energy   = pyfragResult.ReadComporbEnergy({'type': comporb})
               if (HOMO_energy - swpOrb_energy) > 0.01:
                  HOMO_1 = comporb
                  HOMO_low_list.append(HOMO_1)
                  break
         else:
            HOMO_sym          = pyfragResult.GetCompOrbSym({'type': 'HOMO'})
            HOMO_energy       = pyfragResult.ReadComporbEnergy({'type': 'HOMO'})
            for comporb in comp_occupiedMO:
               comp_occupiedMO.remove(comporb)
               swpOrb_sym     = pyfragResult.GetCompOrbSym({'type': comporb})
               swpOrb_energy   = pyfragResult.ReadComporbEnergy({'type': comporb})

               if swpOrb_sym != HOMO_sym and (HOMO_energy - swpOrb_energy) > 0.01:
                  HOMO_1 = comporb
                  HOMO_low_list.append(HOMO_1)
                  break

         mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted = PyFragCycle(HOMO_1, t21File, mt21File, inputKeys, [firstmixing, secondmixing])
         thirdmixing, thirdmixing_adjusted = PatternDetermine(mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted)
         # comp_occupiedMO.remove(HOMO_1)

         if bool(thirdmixing):
            Mark = False
         elif len(comp_occupiedMO) == 0:
            mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted = PyFragCycle(HOMO_low_list[0], t21File, mt21File, inputKeys, [firstmixing, secondmixing],  Mark_mix_1=True)
            thirdmixing, thirdmixing_adjusted = PatternDetermine(mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted)
            Mark = False

      return [firstmixing, thirdmixing, secondmixing],[firstmixing_adjusted,thirdmixing_adjusted,secondmixing_adjusted]

      # return [firstmixing, thirdmixing],[firstmixing_adjusted,thirdmixing_adjusted]
class PyFragResult:
   def __init__(self, complexResult, inputKeys, mt21file=None):
      self.complexResult        = complexResult
      self.fragOrb              = complexResult.read('SFOs', 'ifo')
      #symmetry for each orbital of fragments
      self.fragIrrep            = str(complexResult.read('SFOs', 'subspecies')).split()
      #the fragment label for each orbital
      self.orbFragment          = complexResult.read('SFOs', 'fragment')
      #energy for each orbital
      self.orbEnergy            = complexResult.read('SFOs', 'energy')
      #occupation of each orbitals which is either 0 or 2
      self.orbOccupation        = complexResult.read('SFOs', 'occupation')
      #number of orbitals for each symmetry for complex
      self.irrepOrbNumber       = complexResult.read('Symmetry', 'norb')
      #irrep label for symmetry of complex
      self.irrepType            = str(complexResult.read('Symmetry', 'symlab')).split()
      self.coreOrbNumber        = complexResult.read('Symmetry', 'ncbs')
      if mt21file is not None:
         self.complexResult_Modified    = KFFile(mt21file)
      else:
         self.complexResult_Modified    = self.complexResult


   def ConvertList(self, obj):
      #single number in adf t21 is number fommat which is need to convert list
      if type(obj) == list:
         return obj
      else:
         return [obj]

   def GetElectronNumber(self):

      frag1_elecnum = 0
      frag2_elecnum = 0
      for frag, elecnum in zip(self.orbFragment, self.orbOccupation):
         if frag == 1:
            frag1_elecnum += elecnum
         else:
            frag2_elecnum += elecnum

      return frag1_elecnum, frag2_elecnum


   def GetFaIrrep(self):
      #append complex irrep label to each orbital, if symmetry is A, convert self.irrepOrbNum which is float type into list
      irreporbNum = self.ConvertList(self.irrepOrbNumber)
      faIrrepone  = [[irrep for i in range(number)] for irrep, number in zip(self.irrepType, irreporbNum)]
      return  [irrep for sublist in faIrrepone for irrep in sublist]

   def GetOrbNum(self):
      # GetOrbNumbers including frozen core orbitals, this is necessary to read population, list like [3,4,5,12,13,14,15,23,24,26]
      # core orbital number corresponding to each irrep of complex symmetry
      coreOrbNum           = self.ConvertList(self.coreOrbNumber)
      irrepOrbNum          = self.ConvertList(self.irrepOrbNumber)
      orbNumbers = []
      orbSum = 0
      for nrShell, nrCore in zip(irrepOrbNum, coreOrbNum):
         orbSum += (nrShell + nrCore)
         orbNumbers.extend(list(range(orbSum - nrShell + 1, orbSum + 1)))
      return orbNumbers

   def ReadPopulation(self, index):
      orbNumbers =  self.GetOrbNum()
      #populations of all orbitals
      sfoPopul = self.complexResult.read('SFO popul', 'sfo_grosspop')
      return sfoPopul[orbNumbers[index] - 1]

   def GetFragOrbNum(self):
      # GetOrbNumbers including frozen core orbitals, this is necessary to read population, list like [3,4,5,12,13,14,15,23,24,26]
      #core orbital number corresponding to each irrep of complex symmetry
      coreOrbNum           = self.ConvertList(self.coreOrbNumber)
      irrepOrbNum          = self.ConvertList(self.irrepOrbNumber)
      orbNumbers = []
      for nrShell, nrCore in zip(irrepOrbNum, coreOrbNum):
         orbNumbers.extend(range(nrCore + 1, nrShell + nrCore + 1))
      return orbNumbers

   def GetFrontIndex(self, orbSign):
     #convert HOMO/LUMO/HOMO-1/LUMO+1/INDEX into dict {'holu': 'HOMO', 'num': -1}
      for matchString in [r'HOMO(.*)', r'LUMO(.*)']:
         matchObj = re.match(matchString, orbSign)
         if matchObj:
            holu = re.sub(r'(.[0-9]+)',"", matchObj.group())
            num = re.sub(r'([a-zA-Z]+)',"", matchObj.group())
            if num:
               return {'holu': holu, 'num': num}
            else:
               return {'holu': holu, 'num': 0}

   def GetFragNum(self, frag):
      #change frag type like 'frag1' into number like "1" recorded in t21
      #append fragmenttype(like 1 or 2) to each orbital
      fragType = str(self.complexResult.read('Geometry', 'fragmenttype')).split()
      return fragType.index(frag) + 1

   def GetOrbitalIndex(self, orbDescriptor):
      fragOrbnum = self.GetFragNum(orbDescriptor['frag'])
      orbIndex = 0
      if self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'HOMO':
         orbIndex = sorted(range(len(self.orbEnergy)), key = lambda x: self.orbEnergy[x] if (self.orbFragment[x] == fragOrbnum) and self.orbOccupation[x] != 0 else -1.0E+100, reverse=True)[-int(self.GetFrontIndex(orbDescriptor['type'])['num'])]

      elif self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'LUMO':
         orbIndex = sorted(range(len(self.orbEnergy)), key = lambda x: self.orbEnergy[x] if (self.orbFragment[x] == fragOrbnum) and self.orbOccupation[x] == 0 else +1.0E+100)[int(self.GetFrontIndex(orbDescriptor['type'])['num'])]
      return orbIndex

   def ReadPopulation(self, index):
      orbNumbers =  self.GetOrbNum()
      #populations of all orbitals
      sfoPopul = self.complexResult.read('SFO popul', 'sfo_grosspop')
      return sfoPopul[orbNumbers[index] - 1]

   def ReadOverlap(self, index_1, index_2):
      #orbital numbers according to the symmetry of the complex
      faOrb    = self.GetFragOrbNum()
      faIrrep  = self.GetFaIrrep()
      maxIndex = max(faOrb[index_1], faOrb[index_2])
      minIndex = min(faOrb[index_1], faOrb[index_2])
      index = maxIndex * (maxIndex - 1) / 2 + minIndex - 1
      if faIrrep[index_1] == faIrrep[index_2]:
         self.overlap_matrix = self.complexResult.read(faIrrep[index_1], 'S-CoreSFO')
         return abs(self.overlap_matrix[int(index)])
      else:
         return 0.000

   def ReadFragorbEnergy(self, index):
      return self.complexResult.read('Ftyp '+str(self.orbFragment[index])+self.fragIrrep[index], 'eps')[self.fragOrb[index]-1] * 27.2114


   def GetCompOrbSym(self, orbDescriptor):
      # get the complex HOMO LUMO etc symmetry
      faIrrep  = self.GetFaIrrep()
      orbIndex = 0
      ComporbOccupy = [self.complexResult.read(irreps,'froc_A') for irreps in self.irrepType]
      ComporbOccupy = [orbital for orbitList in ComporbOccupy for orbital in orbitList]
      ComporbEnergy = [self.complexResult.read(irreps,'eps_A') for irreps in self.irrepType]
      ComporbEnergy = [orbital for orbitList in ComporbEnergy for orbital in orbitList ]

      if self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'HOMO':
         orbIndex = sorted(range(len(ComporbEnergy)), key = lambda x: ComporbEnergy[x] if ComporbOccupy[x] != 0 else -1.0E+100, reverse=True)[-int(self.GetFrontIndex(orbDescriptor['type'])['num'])]

      elif self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'LUMO':
         orbIndex = sorted(range(len(ComporbEnergy)), key = lambda x: ComporbEnergy[x] if ComporbOccupy[x] == 0 else +1.0E+100)[int(self.GetFrontIndex(orbDescriptor['type'])['num'])]

      return faIrrep[orbIndex]


   def ReadComporbEnergy(self, orbDescriptor):
      orbIndex = 0
      ComporbOccupy = [self.complexResult.read(irreps,'froc_A') for irreps in self.irrepType]
      ComporbOccupy = [orbital for orbitList in ComporbOccupy for orbital in orbitList]
      ComporbEnergy = [self.complexResult.read(irreps,'eps_A') for irreps in self.irrepType]
      ComporbEnergy = [orbital for orbitList in ComporbEnergy for orbital in orbitList ]

      if self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'HOMO':
         orbIndex = sorted(range(len(ComporbEnergy)), key = lambda x: ComporbEnergy[x] if ComporbOccupy[x] != 0 else -1.0E+100, reverse=True)[-int(self.GetFrontIndex(orbDescriptor['type'])['num'])]

      elif self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'LUMO':
         orbIndex = sorted(range(len(ComporbEnergy)), key = lambda x: ComporbEnergy[x] if ComporbOccupy[x] == 0 else +1.0E+100)[int(self.GetFrontIndex(orbDescriptor['type'])['num'])]
      return ComporbEnergy[orbIndex] * 27.2114


   def GetFragOrbNum_coef(self):
      # GetOrbNumbers including frozen core orbitals, this is necessary to read population, list like [3,4,5,12,13,14,15,23,24,26]
      # core orbital number corresponding to each irrep of complex symmetry
      coreOrbNum           = self.ConvertList(self.coreOrbNumber)
      irrepOrbNum          = self.ConvertList(self.irrepOrbNumber)
      orbNocores      = []
      orbWicores      = []
      orbNumbers      = []
      orbSumWicores   = 0
      orbSum          = 0
      for nrShell, nrCore in zip(irrepOrbNum, coreOrbNum):

         orbSumWicores += nrShell
         orbWicores.append(orbSumWicores)
         orbNocores.append(nrShell + nrCore)
         orbSum        += (nrShell + nrCore)
         orbNumbers.extend(list(range(orbSum - nrShell + 1, orbSum + 1)))
      return orbNumbers, orbNocores, orbWicores


   def ReadOrbCoef(self, orbDescriptor):

      orbIndex      = 0
      ComporbOccupy = [self.complexResult.read(irreps,'froc_A')   for irreps in self.irrepType]
      ComporbOccupy = [orbital for orbitList in ComporbOccupy for orbital in orbitList]
      ComporbEnergy = [self.complexResult.read(irreps,'eps_A') for irreps in self.irrepType]
      ComporbEnergy = [orbital for orbitList in ComporbEnergy for orbital in orbitList]
      ComporbCoef   = [self.complexResult.read(irreps,'Eig-CoreSFO_A') for irreps in self.irrepType]
      ComporbCoef   = [orbital for orbitList in ComporbCoef for orbital in orbitList]

      if self.GetFrontIndex(orbDescriptor['complex'])['holu'] == 'HOMO':
         orbIndex = sorted(range(len(ComporbEnergy)), key = lambda x: ComporbEnergy[x] if ComporbOccupy[x] != 0 else -1.0E+100, reverse=True)[-int(self.GetFrontIndex(orbDescriptor['complex'])['num'])]

      elif self.GetFrontIndex(orbDescriptor['complex'])['holu'] == 'LUMO':
         orbIndex = sorted(range(len(ComporbEnergy)), key = lambda x: ComporbEnergy[x] if ComporbOccupy[x] == 0 else +1.0E+100)[int(self.GetFrontIndex(orbDescriptor['complex'])['num'])]

      fragorbNumber                     =  self.GetOrbitalIndex(orbDescriptor)
      faIrrep  = self.GetFaIrrep()
      orbNumbers, orbNocores, orbWicores=  self.GetFragOrbNum_coef()
      # select irrep orbital number less than orbIndex
      irrepOrbnum_list        =  [irrepNum for irrepNum in orbWicores if irrepNum <= orbIndex]
      # in case orbIndex is smaller than first number
      irrepOrbnum             =  [irrepNum for irrepNum in orbWicores[0:1] if irrepOrbnum_list==[]]
      irrepOrbnum.extend(irrepOrbnum_list)

      coefIndex   = 0
      i = 0
      for indexNumber in orbWicores:
         if len(irrepOrbnum)==1:
            coefIndex  = 0

         elif len(irrepOrbnum)>1 and indexNumber <= irrepOrbnum[-1]:
            coefIndex += orbNocores[i]*self.irrepOrbNumber[i]
            i+=1

      orbCount   = sum(orbNocores[0:i])


      if faIrrep[fragorbNumber] == faIrrep[orbIndex]:

         coefIndex += orbNocores[i]*(orbIndex-orbCount) + fragorbNumber-orbCount
         return ComporbCoef[coefIndex]
      else:
         return 0.00


   def ReadSFO(self, orbDescriptor):

      ###############
      #coeffent part#
      #####
      #####
      #####
      #####
      #####
      ##
      ##
      ###
      ###
      ###
      ###############

      orbIndex      = 0
      ComporbOccupy = [self.complexResult.read(irreps,'froc_A')   for irreps in self.irrepType]
      ComporbOccupy = [orbital for orbitList in ComporbOccupy for orbital in orbitList]
      ComporbEnergy = [self.complexResult.read(irreps,'eps_A') for irreps in self.irrepType]
      ComporbEnergy = [orbital for orbitList in ComporbEnergy for orbital in orbitList]
      ComporbCoef   = [self.complexResult.read(irreps,'Eig-CoreSFO_A') for irreps in self.irrepType]
      ComporbCoef   = [orbital for orbitList in ComporbCoef for orbital in orbitList]

      if self.GetFrontIndex(orbDescriptor['complex'])['holu'] == 'HOMO':
         orbIndex = sorted(range(len(ComporbEnergy)), key = lambda x: ComporbEnergy[x] if ComporbOccupy[x] != 0 else -1.0E+100, reverse=True)[-int(self.GetFrontIndex(orbDescriptor['complex'])['num'])]

      elif self.GetFrontIndex(orbDescriptor['complex'])['holu'] == 'LUMO':
         orbIndex = sorted(range(len(ComporbEnergy)), key = lambda x: ComporbEnergy[x] if ComporbOccupy[x] == 0 else +1.0E+100)[int(self.GetFrontIndex(orbDescriptor['complex'])['num'])]

      fragorbNumber                     =  self.GetOrbitalIndex(orbDescriptor)
      orbNumbers, orbNocores, orbWicores=  self.GetFragOrbNum_coef()

      # select irrep orbital number less than orbIndex
      irrepOrbnum       =  [irrepNum for irrepNum in orbWicores if irrepNum < (orbIndex+1)]

      coefIndex   = 0
      i = 0
      if len(irrepOrbnum)==0:
         coefIndex  = 0
         orbCount   = 0
      else:
         for indexNumber in irrepOrbnum:
            coefIndex += orbNocores[i]*self.irrepOrbNumber[i]
            i+=1
         orbCount  = sum(orbNocores[0:i])

      coefIndex_1 = coefIndex
      coefIndex_2 = coefIndex
      print (coefIndex,'coefIndex')
      print (orbCount,'orbCount')

      faIrrep  = self.GetFaIrrep()
      if faIrrep[fragorbNumber] == faIrrep[orbIndex]:
         coefIndex   += orbNocores[i]*(orbIndex-orbCount) + fragorbNumber-orbCount
         coefIndex_1 += orbNocores[i]*(orbIndex-orbCount)
         coefIndex_2 += orbNocores[i]*(orbIndex-orbCount+1)

      else:
         return 0.00

      ##############
      #overlap part#
      #####
      ####
      ###
      ##
      #
      ###
      ##
      #
      ##############


      #orbital numbers according to the symmetry of the complex
      overlap=[]
      faOrb    = self.GetFragOrbNum()
      orbNumber = faOrb[fragorbNumber]
      overlap_matrix = self.complexResult.read(faIrrep[fragorbNumber], 'S-CoreSFO')
      number_overlap_matrix  = len(overlap_matrix)

      index_horizontal_start = (orbNumber) * (orbNumber - 1)/2
      index_horizontal_end   = (orbNumber) * (orbNumber + 1)/2

      orbitalnumber = int((math.sqrt(1+8*number_overlap_matrix)-1)/2)
      index_vertical = [int((orbNumber+i+1)*(orbNumber+i+2)/2-i-1-1) for i in range(orbitalnumber-orbNumber)]

      overlap.extend(overlap_matrix[int(index_horizontal_start):int(index_horizontal_end)])
      overlap.extend([overlap_matrix[i] for i in index_vertical])

      return sum(np.array(ComporbCoef[coefIndex_1:coefIndex_2])*np.array(overlap))*ComporbCoef[coefIndex]



   '''
   def ReadFragorbEnergy_adjusted(self, orbDescriptor):


      orbIndex           = 0
      ComporbOccupy      = [self.complexResult.read(irreps,'froc_A')   for irreps in self.irrepType]
      ComporbOccupy      = [orbital for orbitList in ComporbOccupy for orbital in orbitList]
      ComporbEnergy_list = [self.complexResult.read(irreps,'eps_A') for irreps in self.irrepType]
      ComporbEnergy      = [orbital for orbitList in ComporbEnergy_list for orbital in orbitList]
      ComporbCoef        = [self.complexResult.read(irreps,'Eig-CoreSFO_A') for irreps in self.irrepType]
      ComporbCoef        = [orbital for orbitList in ComporbCoef for orbital in orbitList]

      if self.GetFrontIndex(orbDescriptor['complex'])['holu'] == 'HOMO':
         orbIndex = sorted(range(len(ComporbEnergy)), key = lambda x: ComporbEnergy[x] if ComporbOccupy[x] != 0 else -1.0E+100, reverse=True)[-int(self.GetFrontIndex(orbDescriptor['complex'])['num'])]

      elif self.GetFrontIndex(orbDescriptor['complex'])['holu'] == 'LUMO':
         orbIndex = sorted(range(len(ComporbEnergy)), key = lambda x: ComporbEnergy[x] if ComporbOccupy[x] == 0 else +1.0E+100)[int(self.GetFrontIndex(orbDescriptor['complex'])['num'])]

      fragorbNumber                     =  self.GetOrbitalIndex(orbDescriptor)
      orbNumbers, orbNocores, orbWicores=  self.GetFragOrbNum_coef()
      # select irrep orbital number less than orbIndex
      irrepOrbnum_list        =  [irrepNum for irrepNum in orbWicores if irrepNum <= orbIndex]
      # in case orbIndex is smaller than first number
      irrepOrbnum             =  [irrepNum for irrepNum in orbWicores[0:1] if irrepOrbnum_list==[]]
      irrepOrbnum.extend(irrepOrbnum_list)
      coefIndex   = 0


      i = 0
      for indexNumber in orbWicores:
         if len(irrepOrbnum)==1:
            coefIndex  = 0

         elif len(irrepOrbnum)>1 and indexNumber <= irrepOrbnum[-1]:
            coefIndex += orbNocores[i]*self.irrepOrbNumber[i]
            i+=1

      # get all complex energy belonging to the symmetry irreducible presentation
      compenergy  = ComporbEnergy_list[i]
      # print (compenergy)

      orbCount     = sum(orbNocores[0:i])
      compcoeflist = []
      coefIndex_1 = 0
      coefIndex_2 = 0
      # get all complex coeffcient belonging to the symmetry irreducible presentation related to the fragment orbital
      for comorbnum in range(orbNocores[i]):
         coefIndex_0 = coefIndex + orbNocores[i]*(comorbnum) + fragorbNumber-orbCount
         coefIndex_1 = coefIndex + orbNocores[i]*(comorbnum)
         coefIndex_2 = coefIndex + orbNocores[i]*(comorbnum+1)

         compcoeflist.append({'start' :coefIndex_1,
                              'end'   :coefIndex_2,
                              'middle':coefIndex_0})



      # get all fragment orbitals overlap belonging to the symmetrical irreducible presentation related to the one fragment orbital
      overlap=[]
      faOrb    = self.GetFragOrbNum()
      orbNumber = faOrb[fragorbNumber]
      faIrrep  = self.GetFaIrrep()
      overlap_matrix = self.complexResult.read(faIrrep[fragorbNumber], 'S-CoreSFO')
      number_overlap_matrix  = len(overlap_matrix)

      index_horizontal_start = (orbNumber) * (orbNumber - 1)/2
      index_horizontal_end   = (orbNumber) * (orbNumber + 1)/2

      orbitalnumber = int((math.sqrt(1+8*number_overlap_matrix)-1)/2)

      index_vertical = [int((orbNumber+i+1)*(orbNumber+i+2)/2-i-1-1) for i in range(orbitalnumber-orbNumber)]

      overlap.extend(overlap_matrix[int(index_horizontal_start):int(index_horizontal_end)])
      overlap.extend([overlap_matrix[i] for i in index_vertical])

      # SFO contibution = coef * overlap
      # adjusted energy = sum[(SFO contibution) * (complex energy)]
      SFOlist = []
      for compcoef in compcoeflist:
         SFO=sum(np.array(ComporbCoef[compcoef['start']:compcoef['end']])*np.array(overlap))*ComporbCoef[compcoef['middle']]
         SFOlist.append(SFO)

      return sum(np.array(SFOlist)*np.array(compenergy)) * 27.2114
   '''


   # def ReadFragorbEnergy_adjusted(self, orbDescriptor):

   #    index      = self.GetOrbitalIndex(orbDescriptor)


   #    return self.energyList[index] * 27.2114
   #    # return self.energyList[index]


   def ReadFragorbEnergy_adjusted(self, orbDescriptor):

      energyList = []
      # when total symmetry, not a list but a number
      for symLab, symNum in zip(self.ConvertList(self.irrepType), self.ConvertList(self.irrepOrbNumber)):
         enerList    = []
         orbenerList = []
         orbList = map(lambda x: int(x*(x+1)/2-1), range(1,symNum+1))
         # deal with situation like "E1:1, E1:2"
         sym = symLab.split(':')[0]
         # read the matrix from the 0 literation, the rest is the same
         enerList = self.complexResult_Modified.read('SFO_Fock', sym)
         orbenerList = [enerList[i] for i in orbList]
         energyList.extend(orbenerList)

      index      = self.GetOrbitalIndex(orbDescriptor)

      return energyList[index] * 27.2114

   def GetOutputData(self, outputData, inputKeys):
      #collect user defined data
      for key, val in list(inputKeys.items()):

         if key == 'population':
            outputData[key] = [self.ReadPopulation(self.GetOrbitalIndex(od)) for od in val]

         elif key == 'orbitalenergy':
            outputData[key] = [self.ReadFragorbEnergy(self.GetOrbitalIndex(od)) for od in val]

         elif key == 'overlap':
            outputData[key] = [self.ReadOverlap(self.GetOrbitalIndex(od1), self.GetOrbitalIndex(od2)) for od1, od2 in val]

         elif key == 'energy':
            outputData[key] = [self.ReadComporbEnergy(od) for od in val]

         elif key == 'orbitalcoeffcient':
            outputData[key] = [self.ReadOrbCoef(od) for od in val]

      return GetOutputTable(outputData)
