import sys,os
import pickle
import plot_env
from ase.db import connect
sys.path.append(os.getcwd()+r'/../../my_analysis_tools')
from analysis_tools import *
import mpl_toolkits.axisartist as axisartist


figsize=(12,6.5)
fontsize=14
y_range=[-5,0]
textshift=0.0

with open(r'../../data/d_band_centers_vac.pckl','rb') as pklin:
    data_vac = pickle.load(pklin)

lattices={}
db = connect(r'../../data/CO_binding.db')
for comp in data_vac.keys():
 if not db.count(#study='CO_binding',
                 name=comp+'slab'):
     print(comp+' not in database')
     lattices[comp]='cubic'
 for row in db.select(name=comp+'slab'):
     lattices[comp] = row.lattice
     break
lattices['Bi2Rh2']='hexagonal'

intdat_vac= {key: data_vac[key] for key in data_vac.keys()}

del intdat_vac['In'] # remove this one, as it is empty


fig,ax=plt.subplots(1,1,figsize=figsize)
ax.axhspan(-2.2, -2.4 ,facecolor='g',alpha=0.2, zorder=-10)

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
plot_alloy_dictionaries(intdat_vac,
                'segregation_energy_vacuum.pdf',
                ylabel=r'd-band center (eV)',
                figax=[fig,ax],

                separators=nperTM, #.values(),
                **kwargs)
print('Number of alloys plotted: '+str(len(intdat_vac.keys())),intdat_vac)

plotted_el=[get_reduced_alloy_name(elem) for elem in intdat_vac.keys()]
plel=[]
plel=plotted_el

ax.set_xticks(np.arange(len(plotted_el)),plel,rotation='vertical',fontsize=fontsize)

plt.savefig('SI_dband_centers.pdf',bbox_inches='tight')
plt.show()
