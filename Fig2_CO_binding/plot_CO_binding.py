#!/usr/bin/python3
import os,sys
from ase.io import read
from general_tools import get_reference_energies,get_crystal_from_info, get_reduced_alloy_name
import numpy as np
import plot_env
import matplotlib
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
sys.path.append(os.getcwd()+f"/../my_analysis_tools/")
from plot_tools import plot_alloy_dictionaries,sort_alloys_for_plot
from interesting_alloys import *
import pickle
import ase.db
from PIL import Image
#from check_calculations import get_raw_OCCO_energies

matplotlib.rcParams['axes.titlesize'] = 24
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 14
matplotlib.rcParams['ytick.labelsize'] = 14
fontsize=25

home = os.getcwd()

#adsorbates=['CO','H','OH']
adsorbates=['CO']
selected_alloys=['InNi','InPd','In3Pd2']
#PTM_for_plot = ['Zn','Cd','Hg','Ga','In','Tl','Ge','Sn','Pb','Sb','Bi']
PTM_all = ['Zn','Cd','Hg','Ga','In','Tl','Ge','Sn','Pb','Sb','Bi']
NM = ['Rh','Ir','Ni','Pd','Pt','Cu','Ag','Au']
NM_for_plot = ['Rh','Ir','Ni','Pd','Pt','Ag','Au']
CO_free_en_shift=0.45

for adsorbate in adsorbates:
    infile = open(r'../data/%s_binding_results.pckl'%adsorbate,'rb')
    data = pickle.load(infile)
    #CHANGE TO THE FOLLOWING IN THE FUTURE
    #import json
    #data = json.load(open('../CO_binding_results.json','r'))
    E_bind_min_all=data[0]
    lattices_all=data[-1]

    E_bind_min={}
    E_bind_min_nonoble={}
    E_bind_min_onlynoble={}
    E_bind_min_int={}
    E_bind_min_really={}
    E_bind_min_sel={}
    E_bind_min_for_exp={}
    lattices,lattices_int,lattices_really,lattices_sel={},{},{},{}
    lattices_nonoble,lattices_onlynoble={},{}

    #Clean binding energy dictionary from empy entries
    delcomps=[]
    for comp in E_bind_min_all.keys():
        for facet in E_bind_min_all[comp].keys():
            if len(E_bind_min_all[comp][facet]):
                break
        else:
            delcomps.append(comp)
            print(comp+'does not contain binding energies. Deleting')

    for comp in delcomps: del E_bind_min_all[comp]

    compsorted,nperTM=sort_alloys_for_plot(PTM_all,E_bind_min_all,PTM_all)
    print(compsorted)

    ci=0
    for elem in compsorted:
      for comp in compsorted[elem]:
        ci+=1
        if comp == 'Hg1Pd1':
            del E_bind_min_all[comp]['112']
            del E_bind_min_all[comp]['001']

        E_bind_min[comp] = E_bind_min_all[comp]
        lattices[comp] = lattices_all[comp]
        #Noble metal only alloys
        #E_bind_min_onlynoble['Cu'] =  E_bind_min_all['Cu']
        #lattices_onlynoble['Cu'] = lattices_all['Cu']
        #E_bind_min_nonoble['Cu'] =  E_bind_min_all['Cu']
        #lattices_nonoble['Cu'] = lattices_all['Cu']
        nonumbercomp=comp
        for i in range(10):
            nonumbercomp = nonumbercomp.replace(str(i),'')
        if (nonumbercomp[:2] not in PTM_all and
            nonumbercomp[2:] not in PTM_all):
                    print(comp)
                    E_bind_min_onlynoble[comp] =  E_bind_min_all[comp]
                    lattices_onlynoble[comp] = lattices_all[comp]
        #Exclude Noble metal only alloys
        else:
                    E_bind_min_nonoble[comp] =  E_bind_min_all[comp]
                    lattices_nonoble[comp] = lattices_all[comp]
        #Interesting alloys
        if get_reduced_alloy_name(comp) in interesting_comps or\
                comp  in interesting_comps:
                    E_bind_min_int[comp] =  E_bind_min_all[comp]
                    lattices_int[comp] = lattices_all[comp]
        #Really intersting alloys
        if get_reduced_alloy_name(comp) in really_interesting_comps or\
                comp  in really_interesting_comps and comp not in ['Pd1Zn1', 'Ga1Ni1']:
                    E_bind_min_really[comp] =  E_bind_min_all[comp]
                    lattices_really[comp] = lattices_all[comp]
        if get_reduced_alloy_name(comp) in selection_for_experimentalists or\
                comp  in selection_for_experimentalists:# and comp not in ['Pd1Zn1', 'Ga1Ni1']:
                    E_bind_min_for_exp[comp] =  E_bind_min_all[comp]
        if len(selected_alloys):
            if get_reduced_alloy_name(comp) in selected_alloys or\
                    comp  in selected_alloys:
                        E_bind_min_sel[comp] =  E_bind_min_all[comp]
                        lattices_sel[comp] = lattices_all[comp]

    E_bind_Cu=E_bind_min['Cu']['100']

    kwargs={    'CO_gas_line':False,
            'binding_on_Cu':E_bind_Cu,
            'PTMs':PTM_all,
            'visualize': False,
            'constant_shift':0.4,
            'ylabel':r'$\Delta$G$_{*\mathrm{CO}}$ / eV'}


    if 1:
        fig,ax=plt.subplots(1,1,figsize=(20,6))
        y_range = [-2,0.5]
        plt.axhspan(0.2, 10 ,facecolor='r',alpha=0.2, zorder=-10)
        plt.axhspan(-0.7, -10 ,facecolor='r',alpha=0.2, zorder=-10)
        plot_alloy_dictionaries(E_bind_min,'%s_binding_energy_all_markers.pdf'%adsorbate,
                **kwargs,
                y_range = y_range,
                figsize=(20,6),
                figax=(fig,ax),
                lattices=lattices,
                show_markers=True,
                markersize=15,
                separators=nperTM,
                fontsize=11)#,xsize=2)
#        img = Image.open("../Color_legend.png").convert("RGBA")
#        ax2=fig.add_axes(plt.axes([0.05, -0.55, 0.5, 0.5], anchor='NE', zorder=1))
#        ax[1].imshow(np.asarray(img))
#        ax[1].axis('off')
#        ax[0].set_position([0.10, 0.3, 0.85, 0.6])
#        ax[1].set_position([0.05, 0.0, 0.85, 0.2])
        #plt.show()

    if 1:
        fig,ax=plt.subplots(1,1,figsize=(10,6))
        y_range = [-1.25,0.5]
        compsorted,nperTM=sort_alloys_for_plot(PTM_all,E_bind_min_int,PTM_all)
        plt.axhspan(0.2, 10 ,facecolor='r',alpha=0.2, zorder=-10)
        plt.axhspan(-0.7, -10 ,facecolor='r',alpha=0.2, zorder=-10)
        plot_alloy_dictionaries(E_bind_min_int,'%s_binding_energy_interesting.pdf'%adsorbate,
                **kwargs,
                y_range = y_range,
    #            figsize=(figsize[0]*10,figsize[1]),
                figax=(fig,ax),
                separators=nperTM,
                markersize=12,
                lattices=lattices_int,
                show_markers=True,
                fontsize=13)
    if 1:
        y_range = [-2,0.5]
        fig,ax=plt.subplots(1,1,figsize=(20,6))
        compsorted,nperTM=sort_alloys_for_plot(PTM_all,E_bind_min_nonoble,PTM_all)
        plot_alloy_dictionaries(E_bind_min_nonoble,'%s_binding_energy_nonoble.pdf'%adsorbate,
                **kwargs,
                y_range = y_range,
                figax=(fig,ax),
                show_markers=True,
                markersize=15,
                lattices=lattices_nonoble,
                separators=nperTM,
                fontsize=fontsize)#,xsize=2)
    if 1:
        y_range = [-2,0.5]
        fig,ax=plt.subplots(1,1,figsize=(10,6))
        plot_alloy_dictionaries(E_bind_min_onlynoble,'%s_binding_energy_onlynoble.pdf'%adsorbate,
                **kwargs,
                y_range = y_range,
                figax=(fig,ax),
                show_markers=False,
                markersize=15,
                lattices=lattices_onlynoble,
                fontsize=fontsize)#,xsize=2)
