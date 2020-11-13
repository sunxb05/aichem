class PyFragResult:
   def __init__(self, complexResult, mt21file=None):
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
