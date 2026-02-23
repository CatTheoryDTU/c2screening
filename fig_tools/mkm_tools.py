#! /usr/bin/env python3
import sys, os
from catmap import analyze
from plot_env import *
import pickle as pkl
#from interesting_alloys import really_interesting_comps
import json
sys.path.append(os.path.abspath('../fig_tools/'))
from general_tools import get_reduced_alloy_name

def add_scaling_line(ax,all_data,facet,pot_state,free_en_corr,plot_line=True):
    #COdat=all_data[facet][pot_state['CO']]
    if pot_state['CO']=='vacuum':
        basedat=all_data['vacuum']['CO']['E_bind_min']
        COdat={}
        for dat in basedat:
            if len(dat) == 2:
                COdat[dat] = basedat[dat]
    else:
        COdat=all_data[facet][pot_state['CO']]
    OCCOdat=all_data[facet][pot_state['OCCO']]
    scaling_fit_data=[]

    #Plot the transition metal points
    for elem in OCCOdat:
        if elem[-2:] in ['_g']: continue
        if elem in ['gas']: continue
        color='k'
        if elem=='Cu': color='r'
        if pot_state['CO'] == 'vacuum':
            codat,occodat=COdat[elem][facet][elem]+free_en_corr['CO'],OCCOdat[elem]['OCCO']+free_en_corr['OCCO']
        else:
            codat,occodat=COdat[elem]['CO']+free_en_corr['CO'],OCCOdat[elem]['OCCO']+free_en_corr['OCCO']
        ax.plot(codat,occodat,'o',color=color)
        scaling_fit_data.append([codat,occodat])
        ax.annotate(f'{elem}',[codat,occodat],color=color).draggable()

    from scipy.optimize import curve_fit
    from general_tools import lin_fun
    scaling_fit_data=np.array(scaling_fit_data)
    coeff,d=curve_fit(lin_fun,scaling_fit_data[:,0],scaling_fit_data[:,1])
    mins,maxs=scaling_fit_data[:,0].min(),scaling_fit_data[:,0].max()
    if plot_line:
        ax.plot([mins,-1.1],[coeff[0]*mins+coeff[1],-1.1*coeff[0]+coeff[1]],'k-')
        ax.plot([-0.35,0.5],[coeff[0]*-0.35+coeff[1],0.5*coeff[0]+coeff[1]],'k-')
        ax.annotate('$\Delta$G$_{*OCCO}$=%1.2f$\Delta$G$_{*CO}$+%1.2feV'%(coeff[0],coeff[1]),
                    #(-1.11,-1.43),rotation=34.0,color='k',fontsize=15).draggable()
                    (-1.11,-1.43),rotation=32.0,color='k',fontsize=15).draggable()

        ax.axvline(x=0.0, color='k', linestyle='dashed')
        CO_OCCO_eq_line=np.array([[-1.5,-3],[-0.75,-1.5]])
        ax.plot(CO_OCCO_eq_line[:,0],CO_OCCO_eq_line[:,1],'k--')
        #ax.annotate('$\Delta$G$_{*OCCO}$=2$\Delta$G$_{*CO}$',(-0.73,-1.46),rotation=43)
        ax.annotate('$\Delta$G$_{*OCCO}$=2$\Delta$G$_{*CO}$',(-0.73,-1.46),rotation=40).draggable()
        CO_OCCO_eq_line=np.array([[-0.28,-0.56],[1,2]])
        ax.plot(CO_OCCO_eq_line[:,0],CO_OCCO_eq_line[:,1],'k--')

def plot_Cu_points(ax,all_data,pot_state,free_en_corr, facets=['100','111','211'], plot_line=True):
    #COdat=all_data[facet][pot_state['CO']]
    for facet in facets:
        elem,color='Cu','r'
        if pot_state['CO']=='vacuum':
            basedat=all_data['vacuum']['CO']['E_bind_min']
            COdat={}
            for dat in basedat:
                if len(dat) == 2:
                    COdat[dat] = basedat[dat]
        else:
            COdat=all_data[facet][pot_state['CO']]
        OCCOdat=all_data[facet][pot_state['OCCO']]
        scaling_fit_data=[]

        #Plot the transition metal points
        if pot_state['CO'] == 'vacuum':
                codat,occodat=COdat[elem][facet][elem]+free_en_corr['CO'],OCCOdat[elem]['OCCO']+free_en_corr['OCCO']
        else:
                codat,occodat=COdat[elem]['CO']+free_en_corr['CO'],OCCOdat[elem]['OCCO']+free_en_corr['OCCO']
        ax.plot(codat,occodat,'o',color=color, markeredgecolor='k', markersize=10)
        ax.annotate(f'{elem}({facet})',[codat,occodat],color=color).draggable()

def add_alloy_points(ax, all_data, pot_state, free_en_corr, really_interesting_comps,tms=None,lattices=None):
    if pot_state['CO']=='vacuum':
        basedat=all_data['vacuum']['CO']['E_bind_min']
        COscadat={}
        for dat in basedat:
            if len(dat) > 2 or dat in ['Cu']:
                COscadat[dat] = basedat[dat]
    else:
        COscadat=all_data['alloys'][pot_state['CO']] #pkl.load(open(COscalinepkl,'rb'))
    scadat=all_data['alloys'][pot_state['OCCO']] #.load(open(scalinepkl,'rb'))
#    scadat.update(all_data['100'][pot_state['OCCO']])
    scaling_fit_data=[]
    #colors=cm.nipy_spectral(np.linspace(0,1,len(scadat.keys())))
    colors=cm.nipy_spectral(np.linspace(0,1,len(really_interesting_comps)))
    markers={'100':'s','111':'h','110':'D','011':'D','001':'s',
                    '11-20':'^','2-1-10':'v','0001':'>',
                 '010':'P','101': 's','121':'X','221':'*','021':'*'}
    #for iel,elem in enumerate(scadat):
    for iel,elem in enumerate(really_interesting_comps):
        if elem[-2:] in ['_g']: continue
        if elem in ['gas']: continue
        color='k'
        if elem=='Cu':
            color='r'
        if elem not in really_interesting_comps: continue
        for facet in scadat[elem].keys():
            if not len(scadat[elem][facet]):
                print(f'No data for {elem} {facet} skipping')
                continue

            ax.plot([],[],markers[facet],color=colors[iel],markeredgecolor='k',markersize=10,label=f'{get_reduced_alloy_name(elem)}({facet})')
            ## When more than one termination exists, choose which to plot - this is a manual hack
            terminations = ['0']
            if elem+facet in ['In1Pd1100','Pd1Zn1001']:
                terminations = ['1']
            elif elem+facet in ['Ga1Ni1100', 'Ge3Pd60001','Ga8Pt4100']:
                terminations = ['0','1']
            ###
            for termination in terminations:
                if pot_state['CO'] == 'vacuum':
                    if elem+facet in ['Ga1Ni1100']:
                        site = 'Ni'
                        if termination == '0':
                            site = 'Ga'
                    elif elem+facet in ['Ga8Pt4100']:

                        site = 'Pt'
                        if termination == '1':
                            site = 'Ga'
                    else:
                      for sit in COscadat[elem][facet].keys():
                        if sit in tms:
                            site = sit
                            break
                    print(f'Using site {site} for {elem} {facet} termination: {termination} for CO (vacuum)')
                    codat = COscadat[elem][facet][site] + free_en_corr['CO']
                else:
                    codat = COscadat[elem][facet][termination]['CO'] + free_en_corr['CO']
                occodat = scadat[elem][facet][termination]['OCCO'] + free_en_corr['OCCO']
                ax.plot(codat,occodat,markers[facet],color=colors[iel],markeredgecolor='k',markersize=14)
                if elem+facet in ['Ga1Ni1100', 'Ga8Pt4100']:
                    ax.annotate(site,[codat,occodat],color='w', ha='center', va='center',fontsize=11).draggable()

def add_annotations(ax):
    ax.annotate('Cu(%s)'%facet,(-2.0,1),ha='left',va='top',color='w',fontsize=40)

def run_catmap(facet,setup_file, runbasedir=os.getcwd()):
    from catmap import ReactionModel
    model = ReactionModel(setup_file = setup_file)
    model.output_variables+=['production_rate', 'free_energy','coverage']
    model.run()
    return model

def plot_heatplot(data,phs,pots,steps,fig,ax,dattype='coverage'):
    X=np.array(sorted(pots))
    Y=np.array(sorted(phs))
    nsteps=len(steps)

    for col in range(1):
     for istep in steps:
        R,S = get_rate_and_selectivity(col,istep,data,nsteps,X,Y)
        plot_it(R,S,fig,ax,col,istep,X,Y,nsteps)

        ax.set_ylabel('$\Delta$G$_{*\mathrm{OCCO}}$ / eV')
        ax.set_xlabel('$\Delta$G$_{*\mathrm{CO}}$ / eV')

def get_rate_and_selectivity(col,istep,data,steps,X,Y):
    Selectivity=np.ones((len(X),len(Y)))*0.5
    rate=np.ones((len(X),len(Y)))*0.5
    for ix,x in enumerate(X):
       for iy,y in enumerate(Y):
        try:
            if col == 1:
                Selectivity[ix][iy]=data[x][y][istep]/np.sum(data[x][y][:nsteps])
            else:
                rate[ix][iy]=data[x][y][istep]
                if rate[ix][iy] <= 1e-10: rate[ix][iy]=1e-8
        except:
            Selectivity[ix][iy]=np.nan
            rate[ix][iy]=1e-20
    return rate, Selectivity

def plot_it(R,S,fig,ax,col,istep,X,Y,nsteps):
        vmin=1e6
        vmax=1e13

        if nsteps == 1:
            thisax=ax
        else:
            thisax=ax[istep]

        for ir in range(R.shape[0]):
            for ic in range(R.shape[1]):
                if R[ir,ic] < 1e-3:
                    R[ir,ic]=1e-3

        b = thisax.imshow(R.T,
                interpolation='bicubic',
                cmap=cm.RdYlGn, #Greens, #jet,
                   origin='lower', extent=[X.min(), X.max(), Y.min(), Y.max()],
                    alpha=0.5,
                    norm=LogNorm(#,
                    vmin=vmin,
                    vmax=vmax),#)
                    aspect='auto')



def extract_data_from_mkm(model,steps,dattype='coverage'):
    data={}
    pots,phs=[],[]

    if dattype=='coverage':
        datin=model.coverage_map
    elif dattype=='production_rate':
        datin=model.production_rate_map

    for dat in datin:
        pot,ph=np.around(dat[0][0],3),np.around(dat[0][1],3)
        if pot not in data:
            data[pot] = {}
        data[pot][ph] = dat[1]
        if pot not in pots:
            pots.append(pot)
        if ph not in phs:
            phs.append(ph)
    return data,phs,pots

def read_data(infile,dattype='production_rate'):
    data_in = pkl.load(open(infile,'rb'),encoding='latin1')
    data={}
    pots,phs=[],[]
    for dat in data_in[f'{dattype}_map']:
        pot,ph=np.around(dat[0][0],3),np.around(dat[0][1],3)
        if pot not in data:
            data[pot] = {}
        data[pot][ph] = dat[1]
        if pot not in pots:
            pots.append(pot)
        if ph not in phs:
            phs.append(ph)
    return data,phs,pots

def run_catmaps_own_analysis(model):
        if not os.path.exists('output'):
            os.mkdir('output')

        vm = analyze.VectorMap(model)
        vm.plot_variable = 'production_rate'
        vm.log_scale = True
        vm.colorbar = True
        vm.min = 1e-5
        vm.max = 1e+2
        fig = vm.plot(save=False)
        fig.savefig('output/production_rate_catmap.pdf')

        vm = analyze.VectorMap(model)
        vm.plot_variable = 'coverage'
        vm.log_scale = True
        vm.colorbar = True
        vm.min = 1e-5
        vm.max = 1e0
        fig = vm.plot(save=False)
        fig.savefig('output/coverage_catmap.pdf')

if __name__ == "__main__":
    main()
