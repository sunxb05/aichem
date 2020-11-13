from plams import *
from orbital_levels import draw_orblevels
from adfreport import draw_orbitals
import argparse as ag
from pyfragdriver  import PyFragDriver
import shutil,sys,os

parser = ag.ArgumentParser(description='Print user defined values')

parser.add_argument("--path", help="Provide path to t21 file.")
parser.add_argument("--mpath", help="Provide path to modified t21 file.")
parser.add_argument("--draworbital", help="Provide path and name to generated data file and figures.", default="No")
parser.add_argument("--adjustorbital", help="Provide path and name to generated data file and figures.", default="No")

args = parser.parse_args()
t21File       = args.path
mt21File      = args.mpath
drawOrbital   = args.draworbital
adjustorbital = args.adjustorbital

dirname   = os.path.split(t21File)[0]
jobname   = os.path.split(dirname)[1]
figuredir = os.path.split(dirname)[0]
figurepath= os.path.join(figuredir,'figures')
figurename= os.path.join(figuredir,'figures',jobname)

if not os.path.exists(figurepath):
  os.mkdir(figurepath)

mixlist, mixlist_adjusted = PyFragDriver(t21File, mt21File)

if drawOrbital=="yes" or drawOrbital=="Yes" or drawOrbital=="Y" or drawOrbital=="y":
   # if os.path.exists(figurepath):
   #    if os.listdir(figurepath):
   #       sys.exit("The figures directory is not empty, delete this figures directory or clean the content in the figures directory to continue.")

   ## If you are sure you will not delete important data, uncomment this and comment the above part to delete figures directory automatically, otherwise be careful to use this delete function.
   shutil.rmtree(figurepath)
   os.mkdir(figurepath)

   draw_orbitals(mixlist,t21File,figurename)

draw_orblevels(mixlist, figurename)

draw_orblevels(mixlist_adjusted,figurename+'_adjusted')

if adjustorbital=="yes" or adjustorbital=="Yes" or adjustorbital=="Y" or adjustorbital=="y":
   draw_orblevels(mixlist_adjusted,figurename+'_adjusted_1')

if mt21File=="yes" or mt21File=="Yes" or mt21File=="Y" or mt21File=="y":
   draw_orblevels(mixlist_adjusted,figurename+'_adjusted_2')
