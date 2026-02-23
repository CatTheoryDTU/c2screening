from ase.io import read
#from adsorb_functions import get_gas_reference
import os,sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from adsorb_functions import *
from plot_tools import *
from catmap_tools import *
from workfunction_tools import *
from my_analysis_tools.analysis_tools import *
from general_tools import *
from ase.units import C
from initialize_stuff import *
#from check_calculations import get_raw_OCCO_energies


def get_CO_results(home,home2,element,facet,ads,slab_energies,slab_wf0,
        Ef_shift_slab,E_bind,E_bind_1V,pzc_data,Capacity,CO_coeffs_vs_q,alloy=False,termination=0,
        mirror=True,linear_fit=True,dpzc_data_CO=None,potential=-1,x_range=[-30,0],zero_volt=4.6,
        gas_references={'C':'CO'}):

          from .analysis_tools import charge_to_pot

          try:
              os.chdir(home2+'/'+ads)
          except FileNotFoundError:
              print(CRED+ads+' directory not found'+CEND)
              return False
          else:
              Eb_vs_pot={}
          if alloy:
            tmplist_all = os.listdir()
            tmplist=[]
            for tmp in tmplist_all:
                  if ads+'-fcc'+facet in tmp.split('_') and tmp.split('_')[0] != 'neglect':
                      if termination == int(tmp.split('_')[0]):
                          tmplist.append(tmp)
          else:
                tmplist = os.listdir()

          #if 1:
          #try:
          CO_results = get_CO_binding_energy(tmplist,element,facet, slab_energies,Ef_shift_slab,adsorbate=ads,mirror=mirror,linear_fit=linear_fit,gas_references=gas_references)

          if not CO_results:
                print(CRED+ads+' binding energy could not be calculated'+CEND)

          else:
              E_bind[ads],Eb_CO_charged,CO_coeffs,site,cellsize = CO_results
              if ads=='CO':
                  for bindingsite in CO_coeffs.keys():
                    CO_coeffs_vs_q[element][facet][bindingsite]=CO_coeffs[bindingsite][0]
              Eb_CO_vs_pot={}
              CO_coeffs_vs_pot={}
              for COsite in Eb_CO_charged.keys():
                CO_ads_pzc,Ef_shift_CO=get_workfunction(tmplist,facet,system=ads+' %s'%COsite.split('_')[-1],site=COsite,charges=Eb_CO_charged[COsite][:,0]/(e2muC/cellsize))
                #print(CO_ads_pzc)
                if len(CO_ads_pzc) > 0:
                    Eb_CO_vs_pot[COsite],CO_coeffs_vs_pot[COsite] = charge_to_pot(element,'CO',facet,slab_wf0,CO_ads_pzc,Eb_CO_charged[COsite],
                            cellsize,Capacity=Capacity)


              Eb_1V = np.array([[ckeys,lin_fun(zero_volt+potential,*CO_coeffs_vs_pot[ckeys])]
                    for ckeys in CO_coeffs_vs_pot.keys()])
              if len(Eb_1V) > 0:
                minEb_ind = np.argsort([float(i) for i in Eb_1V[:,1]])[0]
                E_bind_1V[ads] = float(Eb_1V[minEb_ind,1])
                print (ads+' binding site and strength at %1.2fV:'%potential,Eb_1V[minEb_ind])
              else:
                  print(CRED+ads+' does  not have enough charges'+CEND)

              if alloy:
                plot_Ebind_vs_q(home,element,facet+'_'+str(termination),ads,Eb_CO_charged,CO_coeffs, x_range=x_range,linear_fit=linear_fit,slab_wf0=slab_wf0[0.0],Capacity=Capacity,zero_volt=zero_volt)
              else:
                plot_Ebind_vs_q(home,element,facet,ads,Eb_CO_charged,CO_coeffs, x_range=x_range,linear_fit=linear_fit,slab_wf0=slab_wf0[0.0],Capacity=Capacity,zero_volt=zero_volt)
              #print(ads+' binding energy at %1.2fV: %s'%(potential,E_bind_1V[ads]))
              if ads == 'CO':
                pzc_data['CO'] = slab_wf0[0.0]-zero_volt
              if isinstance(dpzc_data_CO,dict) and ads == 'CO':
                  dpzc_data_CO['CO'] = slab_wf0[0.0]
                  dpzc_data_CO['OCCO'] = CO_ads_pzc[0.0] - slab_wf0[0.0]


def get_CO_binding_energy(tmplist,element,facet,slab_energies,Ef_shift_slab,E_CO=None,mirror=True,GPAW=False,adsorbate='CO',linear_fit=True,gas_references={'C':'CO'}):

    if E_CO == None and GPAW:
        #E_CO = get_gas_reference('CO')
        E_C = get_gas_reference_GPAW('CH4')-2*get_gas_reference_GPAW('H2')
        E_O = get_gas_reference_GPAW('H2O')-get_gas_reference_GPAW('H2')
        E_CO = get_gas_reference_GPAW('CO')
    elif E_CO == None:
        E_CO = get_reference_energies(adsorbate,code='VASP',references=gas_references)

    MCO_toten=[]
    MCO_charged_energies={}
    charges={}
    for dirs in tmplist:
        if os.path.isdir(dirs) and dirs[:4] not in ['qe12','negl'] and dirs[-4:]  not in ['_old']:
         MCO_charged_energies[dirs] = []
         #Get the calculated charges
         for trajfile in os.listdir(os.getcwd()+'/'+dirs+'/'):
             if dirs not in charges.keys():
                 charges[dirs]=[]
             try:
                 if trajfile.split('_')[-3] == 'relaxed':
                    almost_charge=trajfile.split('_')[-1]
                    charges[dirs].append(float(almost_charge.split('.traj')[0]))
             except IndexError:
                 continue

         for charge in charges[dirs]:
            #print( dirs+'/'+element+'_'+adsorbate+'_relaxed_q_%1.2f.traj'%(charge))
            #print(dirs)
            try:
                    #print(os.getcwd()+'/'+dirs+'/'+element+'_'+adsorbate+'_relaxed_q_%1.2f.traj'%(charge))
                    atoms = read(os.getcwd()+'/'+dirs+'/'+element+'_'+adsorbate+'_relaxed_q_%1.2f.traj'%(charge))
            except:
                    #print (dirs+'/'+element+'_'+adsorbate+'_relaxed_q_%1.2f.traj'%(charge))
                    continue
            E = atoms.get_potential_energy()
            if abs(charge) < 0.01:
                MCO_toten.append([dirs,E])
            MCO_charged_energies[dirs].append([charge,E])
         #print(element,adsorbate,MCO_toten)

    #MCO_charged_energies = np.array(MCO_charged_energies)
    #print( MCO_charged_energies)
    #print(MCO_toten)
    MCO_toten = np.array(MCO_toten)
    if len(MCO_toten) > 0:
        #lowest_E_ind = np.argsort(MCO_toten[:,1])[-1]
        lowest_E_ind = np.argsort([float(i) for i in MCO_toten[:,1]])[0]
    else:
        return 0

    #das
    #TODO: Most stable site should also be checked at higher charge!
    print('Most stable %s binding site at q=0:'%adsorbate,MCO_toten[lowest_E_ind,0])

    for slaben in slab_energies:
        if abs(float(slaben[0])) < 0.01:
            E_slab0 = float(slaben[1])

    if GPAW or not  mirror:
            bind_E = (float(MCO_toten[lowest_E_ind,1]) - (E_slab0 + E_CO))
    else:
            #print(E_slab0)
            bind_E = (float(MCO_toten[lowest_E_ind,1]) - (E_slab0 + 2*(E_CO)))/2.

    print( '%s binding energy at q=0:'%(adsorbate),bind_E)

    #charged_tot_energies = MCO_charged_energies[MCO_toten[lowest_E_ind,0]]
    #if len(charges) > 1:
    if 1:
        allEb_CO_charged ={}
        for ch_sites in MCO_charged_energies.keys():
            Eb_CO_charged=[]
            for i,CO_en in enumerate(MCO_charged_energies[ch_sites]):
                for j, slab_en in enumerate(slab_energies):
                    if abs(float(CO_en[0]) - float(slab_en[0])) < 1e-5:
                         Eb_CO_charged.append([CO_en[0],CO_en[1]-float(slab_en[1])])
            if len(Eb_CO_charged) < 2:
                print(CRED+adsorbate+' %s does not seem to be finished'%ch_sites+CEND)
                continue

            Eb_CO_charged = np.array(Eb_CO_charged)
            if GPAW or not mirror:
                Eb_CO_charged[:,1] -= E_CO
            #Eb_CO_charged[:,1] /= 2.
            else:
                Eb_CO_charged[:,1] -= 2*E_CO
                Eb_CO_charged[:,1] /= 2.

            #Add Fermi shift correction

            CO_ads_pzc,Ef_shift_CO=get_workfunction(tmplist,facet,system=adsorbate+' %s'%ch_sites.split('_')[-1],site=ch_sites)

            if len(Ef_shift_CO) > 0:
                for chandE in Eb_CO_charged:
                    chandE[1]+=float(chandE[0])*(Ef_shift_CO[0.0]-Ef_shift_slab[0.0])
            else:
                print(ch_sites,' does not seem to have a neutral calculation')


            Eb_CO_charged[:,0] *= e2muC
            Eb_CO_charged[:,0] /= get_cellsize(atoms)
            allEb_CO_charged[ch_sites] = Eb_CO_charged
        allcoeffs={}
        for csite in allEb_CO_charged.keys():
            Eb_CO_charged=allEb_CO_charged[csite]
            if len(Eb_CO_charged) > 1:
                if linear_fit:
                    allcoeffs[csite],dummy = curve_fit(lin_fun,Eb_CO_charged[:,0],Eb_CO_charged[:,1])
                else:
                    allcoeffs[csite],dummy = curve_fit(quad_fun,Eb_CO_charged[:,0],Eb_CO_charged[:,1])
            else:
                print('Slope of %s binding energy with charge at site %s could not be determined.'%(adsorbate,csite)+
                     'Setting it to 0')
                #print(bind_E)
                allcoeffs[csite]=np.array([0.0,bind_E])
        return bind_E,allEb_CO_charged, allcoeffs,MCO_toten[lowest_E_ind,0],get_cellsize(atoms)#cellsize
            #return bind_E,Eb_CO_charged, np.array([0.0,bind_E]),MCO_toten[lowest_E_ind,0]

        #print( 'CO binding energy at q!=0:',Eb_CO_charged)#, lin_fun(0.0,*coeffs))
    else:
        return bind_E,get_cellsize(atoms)

def get_Cu_reference(adsorbate,facet='100',path='/home/cat/geokast/CO2R/OCCO_volcano/VASP/alloys/Cu_reference',
        gas_references={'C':'CO'}):
    actual_path=path+'/'+facet
    slab_ref = read(actual_path+'/clean_slab/Relaxed.traj').get_potential_energy()
    #Check adsorbed structures
    ads_tot_en={}
    E_bind=[]
    for dirs in os.listdir(actual_path+'/'+adsorbate):
        try:
            ads_tot_en[dirs] = read(actual_path+'/'+adsorbate+'/'+dirs+'/Relaxed.traj').get_potential_energy()
        except FileNotFoundError:
            print(dirs+' does not seem to be finished')
            continue
        else:
            E_bind.append(ads_tot_en[dirs] - (slab_ref+get_reference_energies(adsorbate,references=gas_references)))

    if len(E_bind) > 0:
        return np.min(np.array(E_bind))
    else:
        return None






