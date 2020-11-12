# -*- coding: utf-8 -*-
from energydiagram import ED
import matplotlib.pyplot as plt
import re


def GetFrontIndex(orbSign):
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


def draw_orblevels(mixList,filename):

  energylist = []
  for mixlist in mixList:

    for key, val in list(mixlist['comp_orb_1'].items()):
      energylist.append(val['frag_energy'])

    for key, val in list(mixlist.items()):
      energylist.append(val['fragorb_1']['comp_energy'])

  energygap = abs(max(energylist) - min(energylist))

  orbital = ED(energy_gap=energygap)


  # add levels
  message   = ''
  orbnumber = 0


  if len(mixList) == 1:
    mixlist = mixList[0]
    # add frag1 level
    frag1number = 0
    for key, val in list(mixlist['comp_orb_1'].items()):
      if val['frag_type'] == 'frag1':
        frag1number+=1
        label   = val['frag_orb']
        energy  = val['frag_energy']
        population= val['population']
        message = str(label) + '\n' + '{0:.1f}'.format(energy) + '\n' + '{0:.1f}'.format(population) + ' e'
        orbital.add_level(energy, position='last', top_text=message, color='b')

        # add electrons

        if GetFrontIndex(label)['holu'] == 'HOMO':
          orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=2)
        else:
          orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=0)
        orbnumber+=1


    # add complex level
    i = 0
    for key, val in list(mixlist.items()):
      label   = val['fragorb_1']['comp_orb']
      energy  = val['fragorb_1']['comp_energy']
      message = str(label) + '\n' + '{0:.1f}'.format(energy)
      if i == 0:
        orbital.add_level(energy,                  right_text=message, color='b')
      else:
        orbital.add_level(energy, position='last', right_text=message, color='b')
      i+=1
      # add electrons

      if GetFrontIndex(label)['holu'] == 'HOMO':
        orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=2)
      else:
        orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=0)
      orbnumber+=1


    # add frag2 levels
    i = 0
    for key, val in list(mixlist['comp_orb_1'].items()):
      if val['frag_type'] == 'frag2':
        label   = val['frag_orb']
        energy  = val['frag_energy']
        population= val['population']
        message = str(label) + '\n' + '{0:.1f}'.format(energy) + '\n' + '{0:.1f}'.format(population) + ' e'
        # move one step right
        if i == 0:
          orbital.add_level(energy, top_text=message, color='b')
        else:
          orbital.add_level(energy, position='last', top_text=message, color='b')
        i+=1

        # add electrons

        if GetFrontIndex(label)['holu'] == 'HOMO':
          orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=2)
        else:
          orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=0)
        orbnumber+=1

    # add frag1 and complex link
    mixnumber = len(mixlist)
    i = 0
    for frag1_orb, frag1_val in list(mixlist['comp_orb_1'].items()):
      if frag1_val['frag_type'] == 'frag1':
        k = 0
        for comp_orb, comp_val in list(mixlist.items()):
          coefi = comp_val[frag1_orb]['frag_coef']
          orbital.add_link(i, k+frag1number, color='b', coeficient=coefi, reverse=True)
          k+=1
        i+=1

    # add complex and frag2 link
    i = 0
    for comp_orb, comp_val in list(mixlist.items()):
      k = 0
      for frag2_orb, frag2_val in list(mixlist['comp_orb_1'].items()):
        if frag2_val['frag_type'] == 'frag2':
          coefi = comp_val[frag2_orb]['frag_coef']
          orbital.add_link(i+frag1number, k+frag1number+mixnumber, color='b', coeficient=coefi, reverse=True)
          k+=1
      i+=1

  elif len(mixList) == 2:

    j = 0
    mix1_frag1number = 0
    mix2_frag1number = 0
    for mixlist in mixList:

      # add frag1 level

      for key, val in list(mixlist['comp_orb_1'].items()):
        if val['frag_type'] == 'frag1':
          if j == 0:
            mix1_frag1number+=1
          else:
            mix2_frag1number+=1
          label   = val['frag_orb']
          energy  = val['frag_energy']
          population= val['population']
          message = str(label) + '\n' + '{0:.1f}'.format(energy) + '\n' + '{0:.1f}'.format(population) + ' e'
          if j == 0:
            orbital.add_level(energy, position='last', top_text=message, color='b')
          else:
            orbital.add_level(energy, position='last', top_text=message, color='r')


          if GetFrontIndex(label)['holu'] == 'HOMO':
            orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=2)
          else:
            orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=0)

          orbnumber+=1

      j+=1

    j = 0
    for mixlist in mixList:

      # add complex level
      i = 0
      for key, val in list(mixlist.items()):
        label   = val['fragorb_1']['comp_orb']
        energy  = val['fragorb_1']['comp_energy']
        message = str(label) + '\n' + '{0:.1f}'.format(energy)
        if   j == 0 and i == 0:
          orbital.add_level(energy,                  right_text=message, color='b')
        elif j == 0 and i != 0:
          orbital.add_level(energy, position='last', right_text=message, color='b')
        else:
          orbital.add_level(energy, position='last', right_text=message, color='r')
        i+=1


        if GetFrontIndex(label)['holu'] == 'HOMO':
          orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=2)
        else:
          orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=0)
        orbnumber+=1
      j+=1


    j = 0
    mix1_frag2number = 0
    for mixlist in mixList:
      # add frag2 levels
      i = 0
      for key, val in list(mixlist['comp_orb_1'].items()):
        if val['frag_type'] == 'frag2':
          if j == 0:
            mix1_frag2number += 1
          label   = val['frag_orb']
          energy  = val['frag_energy']
          population= val['population']
          message = str(label) + '\n' + '{0:.1f}'.format(energy) + '\n' + '{0:.1f}'.format(population) + ' e'
          # move one step right
          if   j == 0 and i == 0:
            orbital.add_level(energy,                  top_text=message, color='b')
          elif j == 0 and i != 0:
            orbital.add_level(energy, position='last', top_text=message, color='b')
          else:
            orbital.add_level(energy, position='last', top_text=message, color='r')
          i+=1


          if GetFrontIndex(label)['holu'] == 'HOMO':
            orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=2)
          else:
            orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=0)
          orbnumber+=1
      j+=1


    mix1len = len(mixList[0])
    mix2len = len(mixList[1])

    j = 0
    for mixlist in mixList:
      # add frag1 and complex link
      i = 0
      for frag1_orb, frag1_val in list(mixlist['comp_orb_1'].items()):
        if frag1_val['frag_type'] == 'frag1':
          k = 0
          for comp_orb, comp_val in list(mixlist.items()):
            coefi = comp_val[frag1_orb]['frag_coef']
            if j == 0:
              orbital.add_link(i, k+mix1_frag1number+mix2_frag1number, color='b', coeficient=coefi, reverse=True)
            else:
              orbital.add_link(i+mix1_frag1number, k+mix1_frag1number+mix2_frag1number+mix1len, color='r', coeficient=coefi, reverse=True)
            k+=1
          i+=1
      j+=1


    j = 0
    for mixlist in mixList:
        # add complex and frag2 link
      i = 0
      for comp_orb, comp_val in list(mixlist.items()):
        k = 0
        for frag2_orb, frag2_val in list(mixlist['comp_orb_1'].items()):
          if frag2_val['frag_type'] == 'frag2':
            coefi = comp_val[frag2_orb]['frag_coef']
            if j == 0:
              orbital.add_link(i+mix1_frag1number+mix2_frag1number, k+mix1_frag1number+mix2_frag1number+mix1len+mix2len, color='b', coeficient=coefi, reverse=True)
            else:
              orbital.add_link(i+mix1_frag1number+mix2_frag1number+mix1len, k+mix1_frag1number+mix2_frag1number+mix1len+mix2len+mix1_frag2number, color='r', coeficient=coefi, reverse=True)
            k+=1
        i+=1
      j+=1


  else:
    # j is the order of the mixlist (0-mix2,1-mix3,2-mix1)
    j = 0
    # order is mix2 mix3 mix1
    mix2_frag1number = 0
    mix3_frag1number = 0
    mix1_frag1number = 0

    for mixlist in mixList:

      # add frag1 level

      for key, val in list(mixlist['comp_orb_1'].items()):
        if val['frag_type'] == 'frag1':
          if   j == 0:
            mix2_frag1number+=1
          elif j == 1:
            mix3_frag1number+=1
          elif j == 2:
            mix1_frag1number+=1

          label   = val['frag_orb']
          energy  = val['frag_energy']
          population= val['population']
          message = str(label) + '\n' + '{0:.1f}'.format(energy) + '\n' + '{0:.1f}'.format(population) + ' e'

          if   j == 0:
            orbital.add_level(energy, position='last', top_text=message, color='b')
          elif j == 1:
            orbital.add_level(energy, position='last', top_text=message, color='r')
          elif j == 2:
            orbital.add_level(energy, position='last', top_text=message, color='g')


          if GetFrontIndex(label)['holu'] == 'HOMO':
              orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=2)
          else:
              orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=0)

          orbnumber+=1

      j+=1


    # j is the order of the mixlist (0-mix2,1-mix3, 2-mix1)
    j = 0
    i = 0
    for mixlist in mixList:

      # add complex level
      # i is to locate the position of the first level

      for key, val in list(mixlist.items()):
        label   = val['fragorb_1']['comp_orb']
        energy  = val['fragorb_1']['comp_energy']
        message = str(label) + '\n' + '{0:.1f}'.format(energy)

        if   j == 0 and i == 0:
          orbital.add_level(energy,                  right_text=message, color='b')
        elif j == 0 and i != 0:
          orbital.add_level(energy, position='last', right_text=message, color='b')
        elif j == 1 and i != 0:
          orbital.add_level(energy, position='last', right_text=message, color='r')
        elif j == 2 and i != 0:
          orbital.add_level(energy, position='last', right_text=message, color='g')
        i+=1


        if GetFrontIndex(label)['holu'] == 'HOMO':
          orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=2)
        else:
          orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=0)
        orbnumber+=1
      j+=1


    j = 0
    i = 0
    mix2_frag2number = 0
    mix3_frag2number = 0
    for mixlist in mixList:
      # add frag2 levels
      for key, val in list(mixlist['comp_orb_1'].items()):
        if val['frag_type'] == 'frag2':
          if   j == 0:
            mix2_frag2number += 1
          elif j == 1:
            mix3_frag2number += 1

          label     = val['frag_orb']
          energy    = val['frag_energy']
          population= val['population']
          message = str(label) + '\n' + '{0:.1f}'.format(energy) + '\n' + '{0:.1f}'.format(population) + ' e'

          # move one step right
          if   j == 0 and i == 0:
            orbital.add_level(energy,                  right_text=message, color='b')
          elif j == 0 and i != 0:
            orbital.add_level(energy, position='last', right_text=message, color='b')
          elif j == 1 and i != 0:
            orbital.add_level(energy, position='last', right_text=message, color='r')
          elif j == 2 and i != 0:
            orbital.add_level(energy, position='last', right_text=message, color='g')
          i+=1


          if GetFrontIndex(label)['holu'] == 'HOMO':
            orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=2)
          else:
            orbital.add_electronbox(level_id=orbnumber,boxes=1, electrons=0)
          orbnumber+=1
      j+=1


    mix2len = len(mixList[0])
    mix3len = len(mixList[1])
    mix1len = len(mixList[2])

    j = 0
    for mixlist in mixList:
      # add frag1 and complex link
      i = 0
      for frag1_orb, frag1_val in list(mixlist['comp_orb_1'].items()):
        if frag1_val['frag_type'] == 'frag1':
          k = 0
          for comp_orb, comp_val in list(mixlist.items()):
            coefi = comp_val[frag1_orb]['frag_coef']
            if   j == 0:
              orbital.add_link(i,                                  k+mix1_frag1number+mix2_frag1number+mix3_frag1number,             color='b', coeficient=coefi, reverse=True)
            elif j == 1:
              orbital.add_link(i+mix2_frag1number,                 k+mix1_frag1number+mix2_frag1number+mix3_frag1number+mix2len,             color='r', coeficient=coefi, reverse=True)
            elif j == 2:
              orbital.add_link(i+mix2_frag1number+mix3_frag1number,k+mix1_frag1number+mix2_frag1number+mix3_frag1number+mix2len+mix3len,     color='g', coeficient=coefi, reverse=True)
            k+=1
          i+=1
      j+=1


    j = 0
    for mixlist in mixList:
        # add complex and frag2 link
      i = 0
      for comp_orb, comp_val in list(mixlist.items()):
        k = 0
        for frag2_orb, frag2_val in list(mixlist['comp_orb_1'].items()):
          if frag2_val['frag_type'] == 'frag2':
            coefi = comp_val[frag2_orb]['frag_coef']
            if    j == 0:
              orbital.add_link(i+mix1_frag1number+mix2_frag1number+mix3_frag1number,                 k+mix1_frag1number+mix2_frag1number+mix3_frag1number+mix1len+mix2len+mix3len,                                   color='b', coeficient=coefi, reverse=True)
            elif  j == 1:
              orbital.add_link(i+mix1_frag1number+mix2_frag1number+mix3_frag1number+mix2len,         k+mix1_frag1number+mix2_frag1number+mix3_frag1number+mix1len+mix2len+mix3len+mix2_frag2number,                  color='r', coeficient=coefi, reverse=True)
            elif  j == 2:
              orbital.add_link(i+mix1_frag1number+mix2_frag1number+mix3_frag1number+mix2len+mix3len, k+mix1_frag1number+mix2_frag1number+mix3_frag1number+mix1len+mix2len+mix3len+mix2_frag2number+mix3_frag2number, color='g', coeficient=coefi, reverse=True)
            k+=1
        i+=1
      j+=1



  orbital.plot(show_IDs=True)
  plt.savefig(filename+'level.png',dpi=300)

  plt.show()

  plt.close()
