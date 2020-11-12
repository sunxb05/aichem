from plams import *
from orbital_levels import draw_orblevels
from adfreport import draw_orbitals
import argparse as ag
from readt21  import PyFragDriver
import shutil,sys,os

parser = ag.ArgumentParser(description='Print user defined values')

parser.add_argument("--population", type=str, nargs='*',action='append', help='print population for fragment orbital')
parser.add_argument("--overlap", type=str, nargs='*',action='append', help='print overlap between two fragment orbitals')
parser.add_argument("--orbitalenergy", type=str, nargs='*',action='append', help='print orbital energy')
parser.add_argument("--energy", type=str, nargs='*',action='append', help='print complex orbital energy')
parser.add_argument("--orbitalcoeffcient", type=str, nargs='*',action='append', help='print complex orbital energy')
parser.add_argument("--path", help="Provide path to t21 file.")
parser.add_argument("--mpath", help="Provide path to modified t21 file.")
parser.add_argument("--name", help="Provide path and name to generated data file and figures.")
parser.add_argument("--draworbital", help="Provide path and name to generated data file and figures.", default="No")

inputKeys = {}
for key, val in vars(parser.parse_args()).items():
   if val != None:
      inputValue = []

      if key == 'overlap':
         for term in val:
            inputValue.append(({'type':term[1],'frag':term[0]},{'type':term[3],'frag':term[2]}))
         inputKeys[key] = inputValue

      elif key == 'population':
         for term in val:
            inputValue.append({'type':term[1],'frag':term[0]})
         inputKeys[key] = inputValue

      elif key == 'orbitalenergy':
         for term in val:
            inputValue.append({'type':term[1],'frag':term[0]})
         inputKeys[key] = inputValue

      elif key == 'energy':
         for term in val:
            inputValue.append({'type':term[0]})
         inputKeys[key] = inputValue

      elif key == 'orbitalcoeffcient':
         for term in val:
            inputValue.append({'type':term[2],'frag':term[1],'complex':term[0]})
         inputKeys[key] = inputValue

      else:
         inputKeys[key] = [term for term in val]


args = parser.parse_args()
t21File     = args.path
mt21File    = args.mpath
drawOrbital = args.draworbital


dirname = os.path.split(t21File)[0]
jobname = os.path.split(dirname)[1]
figuredir = os.path.split(dirname)[0]
figurepath= os.path.join(figuredir,'figures')
figurename= os.path.join(figuredir,'figures',jobname)

if not os.path.exists(figurepath):
  os.mkdir(figurepath)


mixlist, mixlist_adjusted = PyFragDriver(t21File, inputKeys, mt21File)

if drawOrbital=="yes" or drawOrbital=="Yes" or drawOrbital=="Y" or drawOrbital=="y":
   # if os.path.exists(figurepath):
   #    if os.listdir(figurepath):
   #       sys.exit("The figures directory is not empty, delete this figures directory or clean the content in the figures directory to continue.")

   ## If you are sure you will not delete important data, uncomment this and comment the above part to delete figures directory automatically, otherwise be careful to use this delete function.
   shutil.rmtree(figurepath)
   os.mkdir(figurepath)

   draw_orbitals(mixlist,t21File,figurename)

draw_orblevels(mixlist, figurename)
# draw_orblevels(mixlist_adjusted,figurename+'adjusted')
