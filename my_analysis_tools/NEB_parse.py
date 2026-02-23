import numpy as np
import os,sys
from general_tools import get_cellsize,get_reference_energies
from ase.io import read
from workfunction_tools import *
from ase.units import C
from  analysis_tools import *

e2muC = 1/(C*1e-6)
sizecoeff = 1e-16
Capacity=20
zero_volt=4.6

def get_NEB_results(basehome,home,element,facet,E_bind,E_bind_pot,pzc_data,slab_wf0,Ef_shift_slab,
        slab_energies=None,mirror=False,cellsize=None,potential=None,linear_fit=True,ads='OC-CO',
        gas_references={'C':'CO'}):

    found_NEB=False
    for dire in os.listdir(home):
        if dire[:3] == 'NEB':
            found_NEB=True
            print('-'*10+'\nChecking NEB: %s '%(dire[4:]))
            os.chdir(dire)
            home2=os.getcwd()
            allbarriers=[]
            allfrontbarriers=[]
            allbackbarriers=[]
            allNEB_energies={}
            allimages={}
            allFS=[]

            #In case the slab is a different on than in the thermocalcs
            #get  its energies
            if slab_energies is None:
                    slabbase=home2+'/clean_slab/'
            elif isinstance(slab_energies,str):
                    slabbase = slab_energies
            else:
                slabbase=None

            if slabbase:
                slab_energies=[]
                for infile in os.listdir(slabbase):
                    if len(infile.split('_')) > 5:
                        if  infile.split('_')[3] == 'relaxed':
                            charge = float(infile.split('_')[-1][:-5])
                            en = read(slabbase+infile).get_potential_energy()
                            slab_energies.append([charge,en])
                slab_energies=np.array(slab_energies)


            barr_indexes={}
            #Get NEB total energies
            for qdire in os.listdir(home2):
                if os.path.isdir(home2+'/'+qdire) and qdire[:2] == 'q_':
                    #First get the  cellsize for charge to charge density
                    if cellsize is None:
                        #cellsize = np.product(np.diag(images[0].cell[:2,:2]))#get_cellsize(read(home2+'/'+qdire+'/climb/00/POSCAR'))
                        cellsize = get_cellsize(read(home2+'/'+qdire+'/00/POSCAR'))

                    nimg=0
                    #Get the number of images
                    for imgdir in os.listdir(home2+'/'+qdire):
#                        print(imgdir[0])
                        if imgdir[0] in ['0','1'] and os.path.isdir(home2+'/'+qdire+'/'+imgdir):
                            if int(imgdir) > nimg:
                                nimg=int(imgdir)
                    NEB_energies=[]
                    images=[]
                    for imgdir in range(nimg+1):
                        #oszilines = open(home2+'/'+qdire+'/climb/%02d/OSZICAR'%imgdir,'r').readlines()
                        #try:
                        #    images.append(read(home2+'/'+qdire+'/climb/%02d/OUTCAR'%imgdir))
                        #except FileNotFoundError:
                        images.append(read(home2+'/'+qdire+'/%02d/OUTCAR'%imgdir))

                        #oszilines = open(home2+'/'+qdire+'/climb/%02d/OSZICAR'%imgdir,'r').readlines()
                        #for line in oszilines:
                        #    if line.split()[1] == 'F=':
                        #        en=float(line.split()[4])
                        NEB_energies.append(images[-1].get_potential_energy())


                    NEB_energies=np.array(NEB_energies)
                    NEB_eb = NEB_energies.copy()
                    allNEB_energies[float(qdire[2:])] = NEB_energies-NEB_energies[0]
                    allimages[float(qdire[2:])] = images



                    for slab_dat in slab_energies:
                        if abs(float(slab_dat[0]) - float(qdire[2:])) < 1e-3:
                            NEB_eb -= float(slab_dat[1])

                    #print(NEB_eb)
                    NEB_eb -= get_reference_energies(ads.replace('-',''),code='VASP',references=gas_references)
                    if mirror:
                        NEB_eb -= get_reference_energies(ads.replace('-',''),code='VASP',references=gas_references)
                        NEB_eb /= 2.
                    #print(NEB_eb)
                    #das

                    #Find the image that is the transition state
                    barr_index = np.argsort(NEB_eb)[-1]
                    #print('NEB binding energies:', NEB_eb,barr_index)
                    barr_indexes[qdire] = barr_index
                    #barr_index=5

                    #allbarriers.append([float(qdire[2:])*e2muC/cellsize,NEB_eb[barr_index]])
                    allbarriers.append([float(qdire[2:]),NEB_eb[barr_index]])
                    allfrontbarriers.append([float(qdire[2:]),NEB_energies[barr_index] - NEB_energies[0]])
                    allbackbarriers.append([float(qdire[2:]),NEB_energies[barr_index] - NEB_energies[-1]])

                    #Get NEB FS binding energy in order to detect barrierless behavior
                    allFS.append([float(qdire[2:]),NEB_eb[-1]])


            print('Number of images in the band:', nimg)
            print('Index of the barrier:', barr_indexes)

            charges=np.linspace(2,-2,21)
            barr_wf,barr_ef_shift = get_workfunction('','dummy',
                    system=ads,charges=charges,verbose=False,barr_index=barr_index)


            #Add Fermi shift difference correction and change charge to charge density
            for chandE in allbarriers:
                chandE[1]+=float(chandE[0])*(barr_ef_shift[float(chandE[0])]-Ef_shift_slab[float(chandE[0])])
                chandE[0]*=e2muC/cellsize

            from .analysis_tools import charge_to_pot
            charge_results = charge_to_pot(
                element,ads,'dummy',slab_wf0,
                barr_wf,allbarriers,cellsize,Capacity=Capacity,
                ref_charge=max(barr_wf))

            if charge_results is not None:
                barr_vs_pot,barr_coeffs_vs_pot = charge_results
            #print(barr_vs_pot,barr_coeffs_vs_pot)

            #Do  the same for the FS in order to detect barrierless behavior
            FS_wf,FS_ef_shift = get_workfunction('','dummy',
                    system=ads,charges=charges,verbose=False,barr_index=len(NEB_eb)-1)

            for chandE in allFS:
                chandE[1]+=float(chandE[0])*(FS_ef_shift[float(chandE[0])]-Ef_shift_slab[float(chandE[0])])
                chandE[0]*=e2muC/cellsize

            charge_results = charge_to_pot(
                element,ads,'dummy',slab_wf0,
                FS_wf,allFS,cellsize,Capacity=Capacity,
                ref_charge=max(FS_wf))

            if charge_results is not None:
                FS_vs_pot,FS_coeffs_vs_pot = charge_results

            #Get Eb vs q
            allbarriers=np.array(allbarriers)
            cbarr,d = curve_fit(lin_fun,allbarriers[:,0],allbarriers[:,1])
            allFS = np.array(allFS)
            cFS,d =  curve_fit(lin_fun,allFS[:,0],allFS[:,1])

            if cbarr[1] <  cFS[1]:
                print('At q=0 the reaction is barrierless')
                E_bind[ads] = cFS[1]
            else:
                E_bind[ads] = cbarr[1]
            #print('OC-CO barrier at q=0: ',cbarr[1]-NEB_eb[0])

            #Get Eb vs pot
            if potential is not None:
                E_bind_barr_pot = lin_fun(potential+zero_volt,*barr_coeffs_vs_pot)
                E_bind_FS_pot = lin_fun(potential+zero_volt,*FS_coeffs_vs_pot)
                    #print(NEB_eb)

            print('OC-CO binding energy at q=0: ',E_bind[ads])

            #Check if the reaction is still barrierless at the given potential
            if E_bind_barr_pot >= E_bind_FS_pot:
                E_bind_pot[ads]=E_bind_barr_pot
            else:
                print('Reaction is still barrierless at %1.2fV'%potential)
                E_bind_pot[ads]=E_bind_FS_pot
            print('OC-CO binding energy at %1.2fV:'%(potential), E_bind_pot['OC-CO'])

            plot_Ebind_vs_q(basehome,element,facet,ads,{'NEB':allbarriers},
                    {'NEB':cbarr},linear_fit=linear_fit,slab_wf0=slab_wf0[0.0],
                    Capacity=Capacity,zero_volt=zero_volt)
            plot_NEB_band(basehome,element,facet,allimages,system=ads)

            #Plot front and back barrier
            names=['_barrier','_reverse_barrier']
            for name_index,f_and_b in enumerate([allfrontbarriers,allbackbarriers]):
                f_and_b=np.array(f_and_b)
                f_and_b[:,0]*=e2muC/cellsize
                cfrontbarrier,d = curve_fit(lin_fun,f_and_b[:,0],f_and_b[:,1])
                plot_Ebind_vs_q(basehome,element,facet,ads+names[name_index],{'NEB':f_and_b},
                    {'NEB':cfrontbarrier},linear_fit=linear_fit,slab_wf0=slab_wf0[0.0],
                    Capacity=Capacity,zero_volt=zero_volt)


            #plot_Ebind_vs_q(basehome,element,facet,ads,{'NEB':allbarriers},
            #        {'NEB':cbarr},linear_fit=linear_fit,slab_wf0=slab_wf0[0.0])
            #plot_front_and_back_barrier(basehome,element,facet,allimages,system=ads)

    if not found_NEB:
        print('No NEB calculations found')

