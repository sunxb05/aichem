import subprocess
from os.path import join as opj
import os

path = os.getcwd()
runscript=opj(path, 'adfreport.sh')
webeditscript_1mixing=opj(path, 'web_1mixing.sh')
webeditscript_2mixing=opj(path, 'web_2mixing.sh')
webeditscript_3mixing=opj(path, 'web_3mixing.sh')

# take second element for sort
def takeSecond(elem):
    return elem[1]


def draw_orbitals(mixList,t21File,fileName):

  frag1orbitallist   = []
  complexorbtiallist = []
  frag2orbitallist   = []
  for mixlist in mixList:
    frag1orbitals    = []
    complexorbtials  = []
    frag2orbitals    = []
    for frag_orb, frag_val in list(mixlist['comp_orb_1'].items()):
      if frag_val['frag_type'] == 'frag1':
        frag1orbitals.append((frag_val['frag_orb'],frag_val['frag_energy']))
      if frag_val['frag_type'] == 'frag2':
        frag2orbitals.append((frag_val['frag_orb'],frag_val['frag_energy']))

    for comp_orb, comp_val in list(mixlist.items()):
      complexorbtials.append((comp_val['fragorb_1']['comp_orb'],comp_val['fragorb_1']['comp_energy']))


    complexorbtials.sort(key=takeSecond,reverse=True)
    complexorbtials = [orbital[0] for orbital in complexorbtials]
    complexOrbtials = str(complexorbtials).replace(',' , ' ').lstrip('[').rstrip(']')
    complexOrbtiallist = str(complexorbtials).replace(',' , ';').lstrip('[').rstrip(']')
    complexorbtiallist.append(complexOrbtiallist)

    if frag1orbitals != []:
      frag1orbitals.sort(key=takeSecond,reverse=True)
      frag1orbitals   = [orbital[0] for orbital in frag1orbitals]
      frag1Orbitals   = str(frag1orbitals).replace(',' , ' ').lstrip('[').rstrip(']')
      frag1Orbitallist= str(frag1orbitals).replace(',' , ';').lstrip('[').rstrip(']')
      frag1orbitallist.append(frag1Orbitallist)
    else:
      frag1Orbitals   = 'None'
      frag1orbitallist.append('None')

    if frag2orbitals != []:
      frag2orbitals.sort(key=takeSecond,reverse=True)
      frag2orbitals   = [orbital[0] for orbital in frag2orbitals]
      frag2Orbitals   = str(frag2orbitals).replace(',' , ' ').lstrip('[').rstrip(']')
      frag2Orbitallist= str(frag2orbitals).replace(',' , ';').lstrip('[').rstrip(']')
      frag2orbitallist.append(frag2Orbitallist)
    else:
      frag2Orbitals   = 'None'
      frag2orbitallist.append('None')

    command = ['bash', runscript, t21File, frag1Orbitals, complexOrbtials, frag2Orbitals, fileName]
    subprocess.call(command)


  if len(mixList)   == 1:
    command = ['bash', webeditscript_1mixing, t21File, frag1orbitallist[0], complexorbtiallist[0], frag2orbitallist[0], fileName]
    subprocess.call(command)

  elif len(mixList) == 2:
    command = ['bash', webeditscript_2mixing, t21File, frag1orbitallist[0], frag1orbitallist[1], complexorbtiallist[0], complexorbtiallist[1], frag2orbitallist[0], frag2orbitallist[1], fileName]
    subprocess.call(command)

  elif len(mixList) == 3:
    command = ['bash', webeditscript_3mixing, t21File, frag1orbitallist[0], frag1orbitallist[1], frag1orbitallist[2], complexorbtiallist[0], complexorbtiallist[1], complexorbtiallist[2], frag2orbitallist[0], frag2orbitallist[1], frag2orbitallist[2], fileName]
    subprocess.call(command)
