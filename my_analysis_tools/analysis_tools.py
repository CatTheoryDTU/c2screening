from ase.io import read
#from adsorb_functions import get_gas_reference
import os,sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from adsorb_functions import *
from .plot_tools import *
from .catmap_tools import *
from .workfunction_tools import *
from general_tools import *
from .stable_adsorbate_parser import *
from .NEB_parse import *
#from check_calculations import get_raw_OCCO_energies

e2muC = 1/(C*1e-6)
sizecoeff = 1e-16
#potential=-1
zero_volt=4.6
#Change the units of charge for e/A**2 tp muC/cm**2
from ase.units import C
CRED = '\033[91m'
CEND = '\033[0m'

matplotlib.rcParams['axes.titlesize'] = 24
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 14
matplotlib.rcParams['ytick.labelsize'] = 14
CRED = '\033[91m'
CEND = '\033[0m'

Capacity = 20. #muF/cm2
verbose=False

home = os.getcwd()

scaling_data = {}
scaling_data_neg1V = {}
scaling_data_NEB = {}
scaling_data_NEB_neg1V = {}
wf_data = {}
neutral_pzc_binding_data = {}
change_in_work_function_data = {}
change_in_work_function_data_CO = {}
OCCO_binding_at_neg1V_vs_pzc = {}
OCCO_binding_vs_sigma = {}
CO_coeffs={}
OCCO_coeffs={}

alldicts=[scaling_data,scaling_data_neg1V,wf_data,neutral_pzc_binding_data,change_in_work_function_data,OCCO_binding_at_neg1V_vs_pzc,OCCO_binding_vs_sigma,change_in_work_function_data_CO,CO_coeffs,OCCO_coeffs]



def find_considered_sites(tmplist,facet):
    sites=[]
    for dirs in tmplist:
        if os.path.isdir(dirs) and dirs[:2] != 'qe':
            if dirs[:-6] not in sites and dirs[-3:] != 'old':
                    sites.append(dirs)
    return sites

def _check_CO_distance(atoms):
                C_pos = []
                for i,pos in enumerate(atoms.positions):
                    if atoms.get_chemical_symbols()[i] == 'C':
                        C_pos.append(pos)

                nOCCO=0
                for i,C in enumerate(C_pos):
                    for j,C2 in enumerate(C_pos):
                        zdist=C[2]-C2[2]
                        if zdist < 2 and zdist > 0.0:
                            for nx in [-1,0,1]:
                             for ny in [-1,0,1]:
                                xycell = np.diag(atoms.cell.copy())

                                xydist = np.linalg.norm(C-(C2+np.array([nx,ny,0])*xycell))
                                #print(zdist,nx, ny, xydist)
                                if xydist < 3:
                                    nOCCO += 1

                return nOCCO

def get_M_OCCO_total_energies(tmplist,element,facet,charges=np.linspace(-2.0,2.0,21),plot_energies=False,mirror = True,linear_fit=True,adsorbate='OCCO'):
    sites = find_considered_sites(tmplist,facet)
    OCCO_energies={}
    atoms=None
    if len(sites) == 0:
        print(CRED+'No sites for %s binding have been found'%adsorbate+CEND)
    for site in sites:
        OCCO_energies[site] = []
        for charge in charges:
            filename = '%s/%s_%s_relaxed_q_%1.2f.traj'%(site,element,adsorbate,charge)
            try:
                atoms = read(filename)
            except:
                #print( filename+ ' not found!')
                continue
            #print(site,'%s/%s_%s_relaxed_q_%1.2f.traj'%(site,element,adsorbate,charge))

            if adsorbate == 'OCCO':
                nOCCO = _check_CO_distance(atoms)
                #print(nOCCO)
                nOCCO_max = 1
                if mirror:
                    nOCCO_max = 2

                if nOCCO >= nOCCO_max:#not nsplit:
                    OCCO_energies[site].append([charge,atoms.get_potential_energy()])
            else:
                    OCCO_energies[site].append([charge,atoms.get_potential_energy()])


        if len(OCCO_energies[site]) < 2:
            del OCCO_energies[site]
        #else:
#           OCCO_energies[site] = np.array(OCCO_energies[site])
#           if plot_energies:
#               if linear_fit:
#                   coeffs,dummy = curve_fit(lin_fun,OCCO_energies[site][:,0],OCCO_energies[site][:,1])
#                   plt.plot(OCCO_energies[site][:,0]+[0.0],
#                      lin_fun(OCCO_energies[site][:,0]+[0.0],
#                      *coeffs))
#               else:
#                   coeffs,dummy = curve_fit(quad_fun,OCCO_energies[site][:,0],OCCO_energies[site][:,1])
#                   plt.plot(OCCO_energies[site][:,0]+[0.0],
#                       quad_fun(OCCO_energies[site][:,0]+[0.0],
#                       *coeffs))
#               plt.plot(OCCO_energies[site][:,0],OCCO_energies[site][:,1],'+')
    #print('MOCCO',OCCO_energies)
    #most_stable_site = find_most_stable_site(OCCO_energies)

    #print ('Most stable binding site:',most_stable_site)
    #return OCCO_energies[most_stable_site],most_stable_site,get_cellsize(atoms)
    if atoms is not None:
        return OCCO_energies,get_cellsize(atoms)
    else:
        return None

def get_M_OCCO_binding_energy(tmplist,element,facet,slab_energies,home,
        charges=np.around(np.linspace(-2.0,6.0,41),2),
        plot_energies=False,mirror = True,verbose=False,
        slab_wf0=None,slab_wf1=None,Ef_shift_slab=None,Ef_shift_slab_charged=None,
        Capacity=None,linear_fit=True,adsorbate='OCCO',gas_references={'C':'CO'}):

    results = get_M_OCCO_total_energies(tmplist,element,facet,
              charges=charges,plot_energies=plot_energies,
              mirror = mirror,linear_fit=linear_fit,adsorbate=adsorbate)

    if results is None:
        print(CRED+adsorbate+' total energies could not be retrieved'+CEND)
        return None
    else:
        E_OCCOa,cellsize = results
#    print(results)
    if 1:
    #try:
        E_bind_OCCO ={}
        Eb_0={}
        OCCO_coeffs2={}
        slab_wf_OCCO={}
        Eb_OCCO_v_pot={}
        OCCO_coeffs_vs_pot={}
        #print('slab energies',slab_energies)
        for OCCO_site in E_OCCOa.keys():
         #print(OCCO_site)
         Eb = []
         for OCCO_en in E_OCCOa[OCCO_site]:
            for j, slab_en in enumerate(slab_energies):
                if abs(float(OCCO_en[0]) - float(slab_en[0])) < 1e-5:
                     Eb.append([OCCO_en[0],OCCO_en[1]-float(slab_en[1])])
         Eb = np.array(Eb)
         if len(Eb) == 0:
             print(CRED+OCCO_site+'does not have all the ingredients for the binding energies'+CEND)
             continue
#         print(Eb)
         if mirror:
            adname=adsorbate.replace('O2','OO')
            Eb[:,1] -= 2*get_reference_energies(adname,code='VASP',references=gas_references)   #E_OCCO
            Eb[:,1] /= 2.
         else:
            Eb[:,1] -= get_reference_energies(adname,code='VASP',references=gas_references)   #E_OCCO


         #Add Fermi shift correction to binding energy

         try:
                slab_wf_OCCO[OCCO_site],Ef_shift_OCCO = get_workfunction('',facet,system=adsorbate+' %s'%(OCCO_site.split('_')[-1]),charges=charges,site=OCCO_site,verbose=verbose)
         except:
                 print(CRED+'Workfunction of %s could not be retrieved'%adsorbate)
                 continue
         else:
             # Apply Fermi level shift correction
             # E = E + q*(Ef_shift_ads - Ef_shift_slab)
                for chandE in Eb:
                    chandE[1]+=float(chandE[0])*(Ef_shift_OCCO[float(chandE[0])]-Ef_shift_slab[float(chandE[0])])


         #Find charge closest to zero in both slab and OCCO calculations
         charge_index=100

         for charge in slab_wf0.keys():
             if charge in Eb[:,0]:
         #        print(charge)
                 if abs(charge) < abs(charge_index):
                     charge_index=charge
         #print('refcharge',OCCO_site,charge_index*e2muC/cellsize,charge_index)

         Eb[:,0] *= e2muC
         Eb[:,0] /= cellsize


         #Convert charged results into potential
         Eb_OCCO_v_pot[OCCO_site], OCCO_coeffs_vs_pot[OCCO_site] = \
                 charge_to_pot(element,adsorbate,facet,slab_wf0,
                         slab_wf_OCCO[OCCO_site],Eb, cellsize,
                         Capacity=Capacity,
                         ref_charge=float(charge_index))
        #OCCO_binding_vs_sigma[element] = E_bind_OCCO

    #Extrapolate to zero charge
    #try:
         #print(Eb)
         if len(Eb) > 1:
            if linear_fit:
                OCCO_coeffs2[OCCO_site],dummy = curve_fit(lin_fun,Eb[:,0],Eb[:,1])
                Eb_0[OCCO_site] = lin_fun(0.0,*OCCO_coeffs2[OCCO_site])
            else:
                OCCO_coeffs2[OCCO_site],dummy = curve_fit(quad_fun,E_bind_OCCO[:,0],E_bind_OCCO[:,1])
                Eb_0[OCCO_site] = quad_fun(0.0,*OCCO_coeffs2[OCCO_site])

         E_bind_OCCO[OCCO_site]=Eb
        if len(Eb_0)>0:
            site,Eb_OCCO_0 = sorted(Eb_0.items(), key=lambda kv: kv[1])[0]
        else:
            return 0

        print(adsorbate+' binding energy and site at q=0:', Eb_OCCO_0, site)
        if verbose:
        #if 1:
            print (adsorbate+' binding energies')
            print( E_bind_OCCO)
    #except:
    #    print(CRED+'Calculating OCCO binding energies failed!'+CEND)
    ##    continue
    return Eb_OCCO_0,E_bind_OCCO,OCCO_coeffs2,Eb_OCCO_v_pot, OCCO_coeffs_vs_pot,site,cellsize

def find_most_stable_site(OCCO_energies):
    lowest=0
    for site in OCCO_energies.keys():
        for charge in range(len(OCCO_energies[site])):
            if OCCO_energies[site][charge][0] + 1.00 < 1e5:
                energy = OCCO_energies[site][charge][1]
                if energy < lowest:
                    lowest = energy
                    lowest_site = site
    return lowest_site


def get_slab_energies(tmplist):
    slab_energies = []
    cellsize=None
    for dirs in tmplist:
            filename = dirs
            if all(x in dirs.split('_') for x in ['clean','slab','relaxed']):
                try:
                    atoms = read(dirs)
                except:
                    print(CRED+'Somethings wrong with the slab %s'%(dirs)+CEND)
                else:
                    cellsize = get_cellsize(atoms)
                    slab_energies.append( [dirs.split('_')[-1][:-5],atoms.get_potential_energy()])
    return np.array(slab_energies),cellsize

def charge_to_pot(element,adsorbate,facet,pzc_slab,pzc_ads,E,cellsize,Capacity=20,ref_charge=0,method='Capacity'):

    #Relate charge to potential via just capacity
    Eb_vs_pot = []
    if method == 'Capacity':
        for q_and_E in E:
    #    Eb_vs_pot.append([(q_and_E[0]-ref_charge)/Capacity + (pzc_slab+pzc_ads)/2.,q_and_E[1]])
        #Eb_vs_pot.append([(q_and_E[0]-ref_charge)/Capacity + pzc_slab,q_and_E[1]])
         if 0:
          if adsorbate == 'OCCO':
            print('ref charge',ref_charge*e2muC/cellsize,ref_charge)
            print('q_and_E',E)
            print('pzc_slab[ref_charge]',pzc_slab)

         #Eb_vs_pot.append([(q_and_E[0]-ref_charge*e2muC/cellsize)/Capacity + pzc_slab[ref_charge],q_and_E[1]])
         Eb_vs_pot.append([q_and_E[0]/Capacity + pzc_slab[0.0],q_and_E[1]])
    elif  method ==  'Capacity_with_mean_onset':
        for q_and_E in E:
            if ref_charge in pzc_ads.keys():
                Eb_vs_pot.append([(q_and_E[0]-ref_charge*e2muC/cellsize)/Capacity + (pzc_slab[ref_charge]+pzc_ads[ref_charge])/2.,
                                   q_and_E[1]])
    elif method == 'Mean_potential':
        for i,q_and_E in enumerate(E):
            charge=np.around(q_and_E[0]/(e2muC/cellsize),1)
            Eb_vs_pot.append([(pzc_slab[charge]+pzc_ads[charge])/2.,q_and_E[1]])

    Eb_vs_pot = np.array(Eb_vs_pot)
    if 0:
     if adsorbate == 'OCCO':
        print(Eb_vs_pot)
        da
    if len(Eb_vs_pot) > 1:
            coeffs,dummy = curve_fit(lin_fun,Eb_vs_pot[:,0],Eb_vs_pot[:,1])
    else:
            print ('Could not calulate the slope of the %s binding energy with potential. Setting it to 0'%(adsorbate))
            #coeffs = np.array([0.0,Eb_vs_pot[0,1]])
            coeffs = np.array([0,0])
    return Eb_vs_pot,coeffs


def get_OCCO_results(E_bind,E_bind_1V,pzc_data,OCCO_binding_at_neg1V_vs_pzc,
        scaling_data_neg1V, dpzc_data,home,home2,element,facet,slab_energies,
        verbose,slab_wf0,Ef_shift_slab,
        Capacity,mirror=True,alloy=False,termination=0,linear_fit=True,x_range=None,
        potential=potential,adsorbate='OCCO',gas_references={'C':'CO'}):
    try:
        os.chdir(home2+'/'+adsorbate+'/')
    except:
        print(CRED+adsorbate+' folder doesnt exist'+CEND)

    if alloy:
        tmplist_all = os.listdir()
        tmplist=[]
        for tmp in tmplist_all:
             if 'OCCO-fcc'+facet in tmp.split('_'):
                if termination == int(tmp.split('_')[0]):
                    tmplist.append(tmp)
    else:
        tmplist = os.listdir()

    if len(tmplist) == 0:
        print(CRED+'No folders have been found for the given facet and termination combination')

    if 1:
    #try:
            OCCO_results=\
            get_M_OCCO_binding_energy(tmplist,element,facet,slab_energies,home,
                    verbose=verbose,slab_wf0=slab_wf0, Ef_shift_slab=Ef_shift_slab,
                    Capacity=Capacity,mirror=mirror,linear_fit=linear_fit,adsorbate=adsorbate,
                    gas_references=gas_references)

            if OCCO_results:
                E_bind[adsorbate],Eb_OCCO_v_q, OCCO_coeffs_vs_q,Eb_OCCO_v_pot, all_OCCO_coeffs_vs_pot,site,cellsize =  \
                OCCO_results
            else:
                print(CRED+adsorbate+' binding energy could not be calculated'+CEND)
                return False
            charges=np.around(np.linspace(-2.0,5.0,36),2)
            slab_wf_OCCO,Ef_shift_OCCO = get_workfunction(tmplist,facet,system=adsorbate,charges=charges
                    ,site=site,verbose=False)
    #except:
    #        pzc_data['OCCO'] = 0
    #        print(CRED+'OCCO binding energy could not be calculated'+CEND)
    #        return False
    #else:
            #print('Binding energy vs q:',Eb_OCCO_v_q)
            #Find the most stable site at the given potential
            Ebs_at_pot = [[cof_site,lin_fun(zero_volt+potential,*all_OCCO_coeffs_vs_pot[cof_site])]
                    for cof_site in all_OCCO_coeffs_vs_pot.keys()]
            Ebs_at_pot = np.array(Ebs_at_pot)

            print (adsorbate+' binding sites and strength at %1.2fV:'%potential,Ebs_at_pot)
            #print(zero_volt+potential)
            #print(Ebs_at_pot)
            print('Coefficients with potential',all_OCCO_coeffs_vs_pot)
            minEbs_ind=np.argsort([float(i) for i in Ebs_at_pot[:,1]])[0]
            E_bind_1V[adsorbate]=float(Ebs_at_pot[minEbs_ind,1])
            print(adsorbate+' binding site and energy  at %1.2fV:'%potential, Ebs_at_pot[minEbs_ind])


            pzc_data[adsorbate] = E_bind[adsorbate]#b_OCCO_0
            dpzc_data['CO'] = slab_wf0[max(slab_wf_OCCO)]
            dpzc_data[adsorbate] = slab_wf_OCCO[max(slab_wf_OCCO)]-slab_wf0[max(slab_wf_OCCO)]

            if alloy:
                plot_Ebind_vs_q(home,element,facet+'_%s'%termination,adsorbate,
                    Eb_OCCO_v_q,OCCO_coeffs_vs_q,linear_fit=linear_fit,
                    x_range=x_range,slab_wf0=slab_wf0[0.0],
                    functiontxt=True,Capacity=Capacity,zero_volt=zero_volt)
                #print(Eb_OCCO_v_q)
                #ds
            else:
                plot_Ebind_vs_q(home,element,facet,adsorbate,Eb_OCCO_v_q,
                    OCCO_coeffs_vs_q,linear_fit=linear_fit,x_range=x_range,
                    slab_wf0=slab_wf0[0.0],functiontxt=True,Capacity=Capacity,
                    zero_volt=zero_volt)
            if alloy:
                scaling_data_neg1V[element][facet][termination] = E_bind_1V
                OCCO_binding_at_neg1V_vs_pzc[element][facet][termination] = [slab_wf0[0.0]-zero_volt,Ebs_at_pot[minEbs_ind,1]]
            else:
                scaling_data_neg1V[element] = E_bind_1V
                OCCO_binding_at_neg1V_vs_pzc[element] = [slab_wf0[0.0]-zero_volt,Ebs_at_pot[minEbs_ind,1]]
            #change_in_work_function_data[element][facet] = dpzc_data
            return E_bind,slab_wf_OCCO,pzc_data,OCCO_binding_at_neg1V_vs_pzc,scaling_data_neg1V,dpzc_data,OCCO_coeffs_vs_q


def get_slab_results(home,home2,facet,verbose,alloy=False,termination=0,potential=-1):
    os.chdir(home2+'/clean_slab')
    if alloy:
        for slabdir in os.listdir(os.getcwd()):
          if os.path.isdir(os.getcwd()+'/'+slabdir) and slabdir not in ['backup']:
             if termination == int(slabdir.split('_')[-1]):
                if facet in slabdir.split('_'):
                    os.chdir(os.getcwd()+'/'+slabdir)
    else:
            pass

    slab_energies=None
    tmplist = os.listdir()
    if 1:
    #try:
        slab_energies,cellsize = get_slab_energies(tmplist)
        #cellsize *= 1e-16

    if len(slab_energies) == 0:
    #except:
        print(CRED+'slab structures create a problem'+CEND)
        os.chdir(home)
        return False#continue

    #if 1:

    try:
        slab_wf0,Ef_shift_slab = get_workfunction(tmplist,facet,system='slab',charges=np.around(np.linspace(-2.0,6.0,41),2))
    except:
        print (CRED+"Slab neutral workfunction not found"+CEND)
        os.chdir(home)
    else:
 #   try:
 #       slab_wf1,Ef_shift_slab_charged = get_workfunction(tmplist,facet,system='slab',charges=[0.])
 #   except:
 #       print ("Slab charged workfunction not found")
    #    os.chdir(home)
        print('Potential of zero charge: ',slab_wf0[0.0]-zero_volt)
        print ('Charge density (muC/cm2) at 0V: %1.3f (%1.4f e/cell)' %(Capacity*(-(slab_wf0[0.0]-zero_volt)),Capacity*(-(slab_wf0[0.0]-zero_volt))*(cellsize/e2muC)))
        print ('Charge density (muC/cm2) at %1.2fV (%1.2feV): %1.3f (%1.4f e/cell)'%(potential, zero_volt+potential, Capacity*(potential-(slab_wf0[0.0]-zero_volt)),Capacity*(potential-(slab_wf0[0.0]-zero_volt))*(cellsize/e2muC)))


    if len(slab_energies) < 3:
        try:
            print('Slab energies only contain the following charges:',slab_energies[:,0])
        except:
            print (CRED+'No slab energies have been found'+CEND)
            return False
    if verbose:
        print('Slab energies:', slab_energies)
    os.chdir(home2)
    return slab_energies,slab_wf0,Ef_shift_slab

