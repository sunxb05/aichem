from plams import *
import re
import sys, math, os
import numpy as np


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

def PyFragCycle(complexOrbital, t21File, mt21File=None, mixList=None, Mark_mix_1=False):

   complexResult    = KFFile(t21File)
   pyfragResult     = PyFragResult(complexResult, mt21file=mt21File)

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
