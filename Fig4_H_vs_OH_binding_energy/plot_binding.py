import sys,os
sys.path.append(os.getcwd()+'/../fig_tools/')
from plot_env import *
import pickle
from interesting_alloys import interesting_comps_after_segr
from my_analysis_tools.analysis_tools import *
import mpl_toolkits.axisartist as axisartist



figsize=(10,6.3)
fontsize=10
y_range=[-1.5,1]
textshift=0.0
H_free_en=0.09
OH_free_en=0.3

OH_min,OH_all,OH_tot,slab,lattices_OH = pickle.load(open('../data/OH_binding_results.pckl','rb'))
H_min,H_all,H_tot,slab_H,lattices_H = pickle.load(open('../data/H_binding_results.pckl','rb'))

H= {key: H_min[key] for key in interesting_comps_after_segr if key in H_min.keys()}
OH= {key: OH_min[key] for key in interesting_comps_after_segr if key in OH_min.keys()}
H_lattices={key: lattices_H[key] for key in interesting_comps_after_segr if key in lattices_H.keys()}
OH_lattices={key: lattices_OH[key] for key in interesting_comps_after_segr if key in lattices_OH.keys()}

fig,ax=plt.subplots(2,1,figsize=figsize)
ax[0].axhspan(0.0, -10 ,facecolor='r',alpha=0.2, zorder=-10)
PTM_all = ['Zn','Cd','Hg','Ga','In','Tl','Ge','Sn','Pb','Sb','Bi']

kwargs= {#'y_range': y_range,
         'figsize': figsize,
         'PTMs': PTM_all,
         'fontsize': fontsize,
         'show_markers': True,
         'markersize': 10,
         'zeroline': False,
         #'lattices': lattices,
         'visualize': False,
        }

compsorted,nperTM=sort_alloys_for_plot(PTM_all,H,PTM_all)
plot_alloy_dictionaries(H,
                'H_binding.pdf',
                ylabel=r'$\Delta$G$_{*\mathrm{H}}$ / eV',
                y_range = [-0.5,1],
                figax=[fig,ax[0]],

                separators=nperTM, #.values(),
                lattices=H_lattices,
                constant_shift=H_free_en,
                **kwargs)
plot_alloy_dictionaries(OH,
                'OH_binding.pdf',
                ylabel=r'$\Delta$G$_{*\mathrm{OH}}$ / eV',
                y_range = [-1.0,1.7],
                figax=[fig,ax[1]],
                lattices=OH_lattices,
                separators=nperTM.values(),
                        constant_shift=OH_free_en,
                **kwargs)

plotted_el=[get_reduced_alloy_name(elem) for elem in H.keys()]
plel=[]
for el in plotted_el:
    print(el)
    if el == 'AuCu2Pd': plel.append(r'AuCu$_2$Pd')
    elif len(el.split('_'))==1: plel.append(' '*3+el+' '*1)
    elif len(el.split('_')) == 2: plel.append(' '*3+el+' '*1)
    elif len(el.split('_')) == 3: plel.append(' '*2+el+' '*0)
    else: plel.append(' '*2+el+' '*4)


ax[0].set_xticks(np.arange(len(plotted_el)),plel,rotation='vertical',fontsize=14)
ax[1].axhspan(-0.5, -10 ,facecolor='r',alpha=0.2, zorder=-10)
#ax[0].text(-2.5,0.9,'(a)',fontsize=16)
#ax[1].text(-2.5,1.5,'(b)',fontsize=16)

ax[0].xaxis.set_tick_params(pad=1)
ax[1].xaxis.set_tick_params(top=True)
ax[1].xaxis.set_tick_params(bottom=False)
ax[1].xaxis.set_tick_params(labelbottom=False)
plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.05, wspace=-0.5, hspace=0.33)
# reduce space between subplots
plt.savefig('Fig_4_H_OH_binding.pdf',bbox_inches='tight')
plt.show()









