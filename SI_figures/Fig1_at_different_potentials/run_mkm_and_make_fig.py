#! /usr/bin/env python3
import sys, os
from catmap import analyze
from plot_env import *
import pickle as pkl
#from interesting_alloys import really_interesting_comps
import json
sys.path.append(os.getcwd()+'/../../fig_tools/')
from mkm_tools import *

home=os.getcwd()
facet='100' #sys.argv[1]
potential=-1.50 #sys.argv[2]
potentials=[-1.50, -0.5]
free_en_corr = {'CO': 0.45, 'OCCO': 0.8} #eV
pot_state = {'CO': 'vacuum', 'OCCO': 'pot'} #Use the result at q=0 ('q0') or at the potential ('pot')
if len(sys.argv) < 2:
    print(f'Facet and potential not given assuming facet ({facet}) and potential {potential}V vs SHE')
else:
    facet=sys.argv[1]
    potential=sys.argv[2]

setup_file = r'../../fig_tools/mkm_file.mkm'
mkm_out_pklfile=f'mkm_results_-1.50V_{facet}.pkl'

cov_and_rate={'coverage':['CO_'+facet],'production_rate': ['C2H4_g']}
potentialscale='SHE'
dattype='production_rate'
run_mkm=1

def main():
    for potential in potentials:
        all_data = json.load(open(f'../../data/all_scaling_data_{float(potential):1.2f}V.json','r'))
        all_data['vacuum'] = json.load(open(f'../../data/vacuum_binding_results_tm.json','r'))
        if run_mkm:
            model=run_catmap(facet,setup_file,home)
            steps=[i for i,label in enumerate(model.output_labels[dattype])
               if label in cov_and_rate[dattype]]
            data,phs,pots=extract_data_from_mkm(model,steps,dattype)
            run_catmaps_own_analysis(model)
        else:
            data,phs,pots=read_data(mkm_out_pklfile)
            steps=[0]
        fig,ax=plt.subplots(len(steps),figsize=(10,7))
        plot_heatplot(data,phs,pots,steps,fig,ax,dattype=dattype)
    #    all_data = json.load(open(f'../all_scaling_data_-1.00V.json','r'))
        add_scaling_line(ax,all_data,facet,pot_state,free_en_corr)
    #    all_data = json.load(open(f'../all_scaling_data_-1.50V.json','r'))
    #    add_scaling_line(ax,all_data)
    #    add_alloy_points(ax,pots,phs)

        ax.set_ylim([np.array(phs).min(),np.array(phs).max()])
        ax.set_xlim([np.array(pots).min(),np.array(pots).max()])
        #plt.legend(loc='upper left',bbox_to_anchor=(1,1),fontsize=10)
        #plt.savefig(f'results/{dattype}_{facet}_{potential}V.pdf')
        plt.show()


if __name__ == "__main__":
    main()
