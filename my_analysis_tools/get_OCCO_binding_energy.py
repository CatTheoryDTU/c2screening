from ase.io import read
from adsorb_functions import get_gas_reference
import os,sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from adsorb_functions import *
from analysis_tools import *
#from check_calculations import get_raw_OCCO_energies

matplotlib.rcParams['axes.titlesize'] = 24
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 14
matplotlib.rcParams['ytick.labelsize'] = 14
verbose=False
CRED = '\033[91m'
CEND = '\033[0m'
#print( tmplist)

Capacity = 20.

#OCCO (g)
E_C = get_gas_reference('CH4')-2*get_gas_reference('H2')
E_O = get_gas_reference('H2O')-get_gas_reference('H2')
#E_CO = get_gas_reference('CO')
E_CO = E_C + E_O
E_OCCO = 2*(E_CO)
print('E_OCCO',E_OCCO)


facet = sys.argv[1]
#if facet not in ['100']:
home = os.getcwd()
#path_info = home.split('/')
#element,facet = path_info[-2],path_info[-1]
elements = ['Au','Cu','Ir','Rh','Pt','Ni','Ag','Pd']
#elements = ['Au','Cu','Rh','Pt','Ni','Ag']
#elements = ['Au','Cu','Pt','Ni']
#elements = ['Au','Cu','Rh','Pt','Ni','Ag','Pd']
#elements = ['Ir','Rh','Pt','Ni','Ag','Pd']
#elements = ['Au','Cu','Ir','Rh','Ni','Ag']
#elements = ['Au','Cu','Rh','Pt','Ni','Ag','Pd']
#elements = ['Ir']
#elements = ['Cu']

scaling_data = {}
scaling_data_neg1V = {}
wf_data = {}
neutral_pzc_data = {}
change_in_work_function_data = {}
OCCO_binding_at_neg1V_vs_pzc = {}
OCCO_binding_vs_sigma = {}

#Change the units of charge for e/A**2 tp muC/cm**2
sizecoeff = 1e-16
from ase.units import C
unitcoeff = 1/(C*1e-6)

for element in elements:
    print('-------\n', element)
    try:
        os.chdir(home+'/'+element+'/'+facet+'/')
    except FileNotFoundError:
        print(CRED+element+' has not been calculated'+CEND)
        continue

    home2 = os.getcwd()

    #OCCO*
    #if facet not in ['100']:
    os.chdir(home2+'/OCCO')
    tmplist = os.listdir()
    #if 1:
    try:
        OCCO_energies,site = get_M_OCCO_total_energies(tmplist,element,facet)
    except:
        print(CRED+'OCCO structures create a problem'+CEND)
        os.chdir(home)
        continue

    #print ('OCCO energies:',OCCO_energies)
    if 1:
    #try:
        slab_wf_OCCO = get_workfunction(tmplist,facet,system='OCCO',charge=-0.80,site=site)
    #except:
    #    print ("OCCO workfunction not found")
    #    os.chdir(home)
    #    continue
    os.chdir(home2)

    #slab
    os.chdir(home2+'/clean_slab')
    tmplist = os.listdir()
    if 1:
    #try:
        slab_energies,cellsize = get_slab_energies(tmplist)
    #except:
    #    print(CRED+'slab structures create a problem'+CEND)
    #    os.chdir(home)
    #    continue
    cellsize *= 1e-16

    try:
        slab_wf0 = get_workfunction(tmplist,facet,system='slab',charge=0.00)
    except:
        print ("Neutral workfunction not found")
        os.chdir(home)
    try:
        slab_wf1 = get_workfunction(tmplist,facet,system='slab',charge=-0.80)
    except:
        print ("Charged workfunction not found")
        os.chdir(home)

    if len(slab_energies) < 3:
        try:
            print('Slab energies only contain the following charges:',slab_energies[:,0])
        except:
            print (CRED+'No slab energies have been found'+CEND)
            continue

    #print('Slab energies:', slab_energies)
    os.chdir(home2)


    #Calculate_binding energy of OCCO and extrapolate to 0 charge

    #if 1:
    try:
        E_bind_OCCO =[]
        #print('slab energies',slab_energies)
        for i, OCCO_en in enumerate(OCCO_energies):
            for j, slab_en in enumerate(slab_energies):
                #print ('calc bind', OCCO_en,slab_en)
                if abs(float(OCCO_en[0]) - float(slab_en[0])) < 1e-5:
                     #print ('calc bind2', type(OCCO_en[1]), type(slab_en[1]))
                     E_bind_OCCO.append([OCCO_en[0],OCCO_en[1]-float(slab_en[1])])
        E_bind_OCCO = np.array(E_bind_OCCO)
        E_bind_OCCO[:,1] -= 2*E_OCCO
        E_bind_OCCO[:,1] /= 2.
        E_bind_OCCO[:,0] *= unitcoeff
        E_bind_OCCO[:,0] /= cellsize

        Eb_OCCO_v_pot, OCCO_coeffs_vs_pot = charge_to_pot(element,'OCCO',facet,slab_wf0,E_bind_OCCO,Capacity=Capacity)
        OCCO_binding_vs_sigma[element] = E_bind_OCCO
        if verbose:
            print ('OCCO binding energies')
            print( E_bind_OCCO)
    except:
        print(CRED+'Calculating OCCO binding energies failed!'+CEND)
        continue

    #Extrapolate to zero charge
    try:
        OCCO_coeffs,dummy = curve_fit(lin_fun,E_bind_OCCO[:,0],E_bind_OCCO[:,1])
        Eb_OCCO_0 = lin_fun(0.0,*OCCO_coeffs)
        print('OCCO binding energy at q=0:', Eb_OCCO_0)
    except:
        print(CRED+'Could not extrapolate OCCO binding energy to zero charge'+CEND)
        continue
    else:
        plt.clf()
        plt.plot(E_bind_OCCO[:,0],E_bind_OCCO[:,1],'+')
        plt.plot([-25,0.0],lin_fun(np.array([-25,0.0]),*OCCO_coeffs))
        plt.text(-5,E_bind_OCCO[-1,1],'%1.4f*x+%1.4f'%(OCCO_coeffs[0],OCCO_coeffs[1]))
        plt.savefig(home+'/results/OCCO_binding_vs_q_%s%s.png'%(element,facet))
        plt.clf()

    #Get CO binding energy
    os.chdir(home2+'/CO')
    tmplist = os.listdir()
    if 1:
    #try:
        Eb_CO = get_CO_binding_energy(tmplist,element,slab_energies,cellsize,E_CO=E_CO,charges=[0.0])
    #except:
    #    print(CRED+'CO binding energy create problems'+CEND)
    #    os.chdir(home)
    #    continue

    if 1:
        Eb_CO,Eb_CO_charged,CO_coeffs = get_CO_binding_energy(tmplist,element,slab_energies,cellsize)
        Eb_CO_v_pot,CO_coeffs_v_pot = charge_to_pot(element,'CO',facet,slab_wf0,Eb_CO_charged,Capacity=Capacity)

    if 0:
        plt.clf()
        plt.plot(Eb_vs_pot[:,0],Eb_vs_pot[:,1],'+')
        plt.plot([3.2,5.2],lin_fun(np.array([3.2,5.2]),*coeffs))
        plt.savefig(home+'/results/Eb_%s_vs_pot_%s_%s.png'%(adsorbate,element,facet))
        plt.clf()
#    print(Eb_CO_charged)
#    das
    scaling_data[element]=[Eb_CO,Eb_OCCO_0]
    neutral_pzc_data[element] = [slab_wf0-4.4,Eb_OCCO_0]
    wf_data[element] = [get_vacuum_wf(element,facet),Eb_OCCO_0]
    change_in_work_function_data[element] = [slab_wf0,slab_wf_OCCO-slab_wf1]


    if 1:
        scaling_data_neg1V[element] = [lin_fun(3.4,*CO_coeffs_v_pot),lin_fun(3.4,*OCCO_coeffs_vs_pot)]
        OCCO_binding_at_neg1V_vs_pzc[element] = [slab_wf0-4.4,lin_fun(3.4,*OCCO_coeffs_vs_pot)]
    os.chdir(home)

write_pckl_data(home+'/results/scaling_data_%s.pckl'%facet,scaling_data)
write_pckl_data(home+'/results/pzc_data_%s.pckl'%facet,neutral_pzc_data)
write_pckl_data(home+'/results/scaling_data_-1V_%s.pckl'%facet,scaling_data_neg1V)

plt.clf()
plot_scaling_line(scaling_data,facet,home,title = '%s@q=0'%facet)
plot_scaling_line(neutral_pzc_data,facet,home,basename='pzc_vs_bindE',xname='$\Phi_{pzc} [V]$')
plot_scaling_line(wf_data,facet,home,basename='wf_vs_bindE',xname='$\Phi_{vac} [V]$')
plot_scaling_line(change_in_work_function_data,facet,home,basename='dwf_vs_pzc',xname='$\Phi_{pzc} [V]$',yname='$\Phi_{*OCCO}-\Phi_{*}$')
if 1:
    plot_scaling_line(scaling_data_neg1V,facet,home,basename = 'scaling_-1V', title = '%s@-1V'%facet)
    plot_scaling_line(OCCO_binding_at_neg1V_vs_pzc,facet,home,basename='OCCO_binding_at_neg1V_vs_pzc',xname='$\Phi_{pzc} [V]$',title = '%s@-1V'%facet)

create_catmap_input(home,element,facet,scaling_data)
create_catmap_input(home,element,facet,scaling_data_neg1V,basename='/catmap/catmap_input_-1V')



