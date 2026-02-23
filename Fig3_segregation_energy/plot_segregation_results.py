import sys,os
#from matplotlib import pyplot as plt
import pickle
#from ase.io import read,write
#import numpy as np
#from general_tools import *
import plot_env
#from ase.db import connect
from interesting_alloys import interesting_comps
sys.path.append('../my_analysis_tools')
from analysis_tools import *
#from pylab import rcParams
import mpl_toolkits.axisartist as axisartist


def setup_axes(fig, pos):
    ax = fig.add_subplot(pos, axes_class=axisartist.Axes)
    return ax

#home = os.getcwd()

figsize=(10,6.3)
fontsize=14
y_range=[-1.5,4]
textshift=0.0


lattices={}
#db = connect('segregation.db')

#for comp in interesting_comps:
# if not db.count(study='segregation alloys',
#                 name=comp):
#     print(comp+' not in database')
#     lattices[comp]='tetragonal'
# for row in db.select(name=comp
#     lattices[comp] = row.lattice
#     break

with open('segregation_results.pckl','rb') as pklin:
    data_vac = pickle.load(pklin)
with open('segregation_with_CO_binding_results.pckl','rb') as pklin:
    data_with_CO = pickle.load(pklin)


intdat_vac= {key: data_vac[key] for key in interesting_comps if key in data_vac.keys()}
intdat_with_CO= {key: data_with_CO[key] for key in interesting_comps if key in data_with_CO.keys()}

del intdat_vac['Pt4Sb8'] # remove this one, as it is empty
del intdat_with_CO['Pt4Sb8'] # remove this one, as it is empty

fig,ax=plt.subplots(2,1,figsize=figsize)
ax[0].axhspan(0.0, -10 ,facecolor='r',alpha=0.2, zorder=-10)
PTM_all = ['Zn','Cd','Hg','Ga','In','Tl','Ge','Sn','Pb','Sb','Bi']

kwargs= {'y_range': y_range,
         'figsize': figsize,
         'PTMs': PTM_all,
         'fontsize': fontsize,
         'show_markers': True,
         'markersize': 10,
         'zeroline': True,
         'lattices': lattices,
         'visualize': False,
         'alpha': 1,
        }

compsorted,nperTM=sort_alloys_for_plot(PTM_all,intdat_vac,PTM_all)
# Find alloys that are only in one of the
# two datasets and remove them from both
keys_vac=set(intdat_vac.keys())
keys_with_CO=set(data_with_CO.keys())
keys_to_remove=keys_vac.symmetric_difference(keys_with_CO)
print('Removing alloys not in both datasets: ',keys_to_remove)
for key in keys_to_remove:
    if key in intdat_vac:
        del intdat_vac[key]
    if key in data_with_CO:
        del data_with_CO[key]

# Make sure both datasets are sorted the same way
co_dat={}
for i in intdat_vac.keys():
    co_dat[i]=data_with_CO[i]

#print(intdat_vac.keys())
#print(co_dat.keys())
#dd
print(intdat_vac['Hg1Pd1'])
print(co_dat['Hg1Pd1'])
del intdat_vac['Hg1Pd1']['101']

plot_alloy_dictionaries(intdat_vac,
#plot_alloy_dictionaries(compsorted,
                'segregation_energy_vacuum.pdf',
                ylabel=r'$\Delta$E$_{seg}^{vac}$ / eV',
                figax=[fig,ax[0]],
                separators=nperTM, #.values(),
                **kwargs)
plot_alloy_dictionaries(co_dat,
                'segregation_energy_CO.pdf',
                ylabel=r'$\Delta$E$_{seg}^{CO}$ / eV',
                figax=[fig,ax[1]],
                separators=nperTM.values(),
                **kwargs)
print('Number of alloys plotted: '+str(len(intdat_vac.keys())),intdat_vac)
print('Number of alloys with CO plotted: '+str(len(intdat_with_CO.keys())),intdat_with_CO)

plotted_el=[get_reduced_alloy_name(elem) for elem in intdat_vac.keys()]
plel=[]
for el in plotted_el:
    print(el)
    if el == 'AuCu2Pd': plel.append(r'AuCu$_2$Pd')
    elif len(el.split('_'))==1: plel.append(' '*3+el+' '*3)
    elif len(el.split('_')) == 2: plel.append(' '*3+el+' '*3)
    elif len(el.split('_')) == 3: plel.append(' '*2+el+' '*2)
    else: plel.append(' '*2+el+' '*4)


ax[0].set_xticks(np.arange(len(plotted_el)),plel,rotation='vertical',fontsize=fontsize)
ax[1].axhspan(0.0, -10 ,facecolor='r',alpha=0.2, zorder=-10)

ax[0].xaxis.set_tick_params(pad=1)
ax[1].xaxis.set_tick_params(top=True)
ax[1].xaxis.set_tick_params(bottom=False)
ax[1].xaxis.set_tick_params(labelbottom=False)
plt.subplots_adjust(left=0.08, right=0.95, top=0.93, bottom=0.05, wspace=-0.2, hspace=0.47)
# reduce space between subplots
plt.savefig('Fig3_segregation_energies.pdf',bbox_inches='tight')
plt.show()
