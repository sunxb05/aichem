from pyfragresult import PyFragResult
import plams
from pyfragcycle import PyFragCycle, PatternDetermine

def PyFragDriver(t21File, mt21File=None):
   '''
   'overlap': [({'type': 'LUMO+2', 'frag': 'frag1'}, {'type': 'LUMO+2', 'frag': 'frag2'}),({'type': 'LUMO', 'frag': 'frag1'}, {'type': 'HOMO', 'frag': 'frag2'})]

   'orbitalcoeffcient': [{'type': 'LUMO+2', 'frag': 'frag1', 'complex': 'LUMO+2'}, {'type': 'LUMO+2', 'frag': 'frag2', 'complex': 'LUMO+2'}]

   'energy': 'energy': [{'type': 'LUMO+2'}, {'type': 'LUMO+1'}]

   'orbitalenergy': [{'type': 'LUMO+2', 'frag': 'frag1'}, {'type': 'LUMO+1', 'frag': 'frag1'}]
   '''

   resultsTable     = []
   outputData       = {}

   complexResult    = KFFile(t21File)
   pyfragResult     = PyFragResult(complexResult, mt21file=mt21File)

   HOMO_orbitals     = ['HOMO', 'HOMO-1', 'HOMO-2', 'HOMO-3', 'HOMO-4', 'HOMO-5', 'HOMO-6', 'HOMO-7', 'HOMO-8', 'HOMO-9']
   frag1_elecnum, frag2_elecnum = pyfragResult.GetElectronNumber()
   frag1_orbnum      =  int(frag1_elecnum/2)
   frag2_orbnum      =  int(frag2_elecnum/2)
   comp_orbnum       =  frag1_orbnum + frag2_orbnum

   if comp_orbnum <= 10:
      comp_occupiedMO= HOMO_orbitals[1:comp_orbnum]
   else:
      comp_occupiedMO= HOMO_orbitals[1:10]


   mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted= PyFragCycle('HOMO', t21File, mt21File)
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
                  # HOMO_low_list.append(HOMO_1)
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
                  # HOMO_low_list.append(HOMO_1)
                  break

         mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted = PyFragCycle(HOMO_1, t21File, mt21File, [firstmixing])
         thirdmixing, thirdmixing_adjusted = PatternDetermine(mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted)
         HOMO_low_list.append([thirdmixing,thirdmixing_adjusted])
         if bool(thirdmixing):
            Mark = False
         elif len(comp_occupiedMO) == 0:
            # mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted = PyFragCycle(HOMO_low_list[0], t21File, mt21File, [firstmixing, secondmixing],  Mark_mix_1=True)
            # thirdmixing, thirdmixing_adjusted = PatternDetermine(mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted)
            population = 0
            mix_choose = HOMO_low_list[0]

            for mixing in HOMO_low_list[1:]:
               population_1 = mix_choose[0]["comp_orb_1"]['fragorb_1']['population']
               population   = mixing[0]["comp_orb_1"]['fragorb_1']['population']
               if population > 1 and population_1 > 1:
                  # population donated
                  if population < population_1:
                     mix_choose = mixing
               elif opulation < 1 and population_1 < 1:
                  # population accepted
                  if population > population_1:
                     mix_choose = mixing
               elif population > 1 and population_1 < 1:
                  if (2-population) > population_1:
                     mix_choose = mixing
               else:
                  if population > (2 - population_1):
                     mix_choose = mixing
            thirdmixing          = mix_choose[0]
            thirdmixing_adjusted = mix_choose[1]

            Mark = False

      return [firstmixing, thirdmixing],[firstmixing_adjusted,thirdmixing_adjusted]

   else:

      mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted= PyFragCycle('LUMO', t21File, mt21File, [firstmixing])
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
                  # HOMO_low_list.append(HOMO_1)
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
                  # HOMO_low_list.append(HOMO_1)
                  break

         mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted = PyFragCycle(HOMO_1, t21File, mt21File, [firstmixing, secondmixing])
         thirdmixing, thirdmixing_adjusted = PatternDetermine(mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted)
         HOMO_low_list.append([thirdmixing,thirdmixing_adjusted])
         if bool(thirdmixing):
            Mark = False
         elif len(comp_occupiedMO) == 0:
            # mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted = PyFragCycle(HOMO_low_list[0], t21File, mt21File, [firstmixing, secondmixing],  Mark_mix_1=True)
            # thirdmixing, thirdmixing_adjusted = PatternDetermine(mixList_1, mixList_2, mixList_3, mixList_1_adjusted, mixList_2_adjusted, mixList_3_adjusted)

            population = 0
            mix_choose = HOMO_low_list[0]

            for mixing in HOMO_low_list[1:]:
               population_1 = mix_choose[0]["comp_orb_1"]['fragorb_1']['population']
               population   = mixing[0]["comp_orb_1"]['fragorb_1']['population']
               if population > 1 and population_1 > 1:
                  # population donated
                  if population < population_1:
                     mix_choose = mixing
               elif opulation < 1 and population_1 < 1:
                  # population accepted
                  if population > population_1:
                     mix_choose = mixing
               elif population > 1 and population_1 < 1:
                  if (2-population) > population_1:
                     mix_choose = mixing
               else:
                  if population > (2 - population_1):
                     mix_choose = mixing
            thirdmixing          = mix_choose[0]
            thirdmixing_adjusted = mix_choose[1]

            Mark = False

      return [firstmixing, thirdmixing, secondmixing],[firstmixing_adjusted,thirdmixing_adjusted,secondmixing_adjusted]
