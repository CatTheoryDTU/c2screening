import os,sys
import ase.calculators.vasp as vasp_calculator
import numpy as np
from ase import *
from ase.visualize import view
from ase.io import read, write
from ase.build import molecule
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo

mol = 'H2'

if mol in ['CO','H2']:
    geo = 'linear'
else:
    geo = 'nonlinear'

atoms=molecule(mol)
atoms.pbc = [1,1,1]
atoms.center(vacuum=7)
#view(atoms)

write('init.traj', atoms)


calc = vasp_calculator.Vasp(
           #Performance
           istart=1,  #Start from WAVECAR
           #ncore=40,  #Parallelization
           #kpar=1,    #Parallelization over nodes
           npar =4,

           #DFT parameters
           encut=500,  #PW cutoff
#           kpts  = (4,3,1),  #kpts
           kpts = (1,1,1),
#           gamma = False,    #Enforce gamma point
           xc='PBE',         #xc for pseudopotential
           gga='BF',         #actual xc
           luse_vdw=True,
           zab_vdw=-1.8867,
           lbeefens=True,
           ismear=-1,         #Fermi smearing
           sigma = 0.1,       #Fermi temperature
           ldipol=False,
           #dipol=(0.5,0.5,0.5),
           #idipol=3,
           #Structure optimization
           ibrion=2,   #optimization method 2=cg

           #Convergence
           algo = 'Normal',
           ediffg=-0.025,
           ediff=1e-4,
           prec='Accurate',
           nelm=350,
           nsw=200,  #Maximunm number of ionic steps

           #OUTPUT
#           lvhar=False,       #write hartree potential
#           lcharg=True,       #write CHARGECAR

           #VASPSOL
#           lsol=True,         #Vaspsol
#           lambda_d_k=3.0,    #Debye screening length
#           tau=0,             #surface tension (apparently creates problem if nonzero)
#           nelect= nelecs0-2*charge,      #Total number of electrons
    #       lasph=True,
    #       lorbit=11,
                   )

#atoms.rattle()
atoms.set_calculator(calc)

#atoms.get_potential_energy()
#atoms.calc.set(xc = 'BEEF-vdW')
#atoms.get_potential_energy()
#from ase.optimize import BFGS
#qn = BFGS(atoms,trajectory = 'relax.traj',logfile = 'relax.log')
#qn.run(fmax=0.03)
#E1=atoms.get_potential_energy()

potentialenergy = atoms.get_potential_energy()
write(mol+'_relaxed.traj',atoms)
sys.exit()
os.chdir('vibs')

vib = Vibrations(atoms)
vib.run()
vib_energies = vib.get_energies()

thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=potentialenergy,
                        atoms=atoms,
                        geometry='linear',
                        symmetrynumber=1, spin=0)
G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325.)
