from ase.io import read
import os,sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit
from adsorb_functions import *
from .analysis_tools import *
from general_tools import *
#from check_calculations import get_raw_OCCO_energies

#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#plt.rcParams["font.family"] = "Times New Roman"
#plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
#plt.rcParams['xtick.labelsize'] = 18
#plt.rcParams['ytick.labelsize'] = 18
#plt.rcParams['figure.figsize'] = (7,5)
#markersize=10

def lin_fun(x,a,b):
    return a*x+b

def quad_fun(x,a,b,c):
    return a*x*x+b*x+c

#def get_cellsize(atoms):
#    return np.product(np.diag(atoms.cell[:2,:2].copy()))*sizecoeff

def plot_scaling_line(scaling_data,facet,home,basename='scaling_line',
        xname=None,yname=None, title = None,alloy=False,
        return_plt=False,plotter=False,fit_line=True,txt=True,pointcolor='.r',
        fontsize=15,legend=False,different_colors=False,markersymbol='o',descriptors=['CO','OCCO'],
        single_marker_per_alloy=False,separate_facets=True,single_legend_per_alloy=False,
        xrange=[-2.5,0.2],yrange=[-3,0.5],figsize=None):
    if xname is None:
        xname='$\Delta$E$_{%s}$ [eV]'%descriptors[0]
    if yname is None:
        yname='$\Delta$E$_{%s}$ [eV]'%descriptors[1]

    if figsize is not None:
        from pylab import rcParams
        rcParams['figure.figsize'] = figsize
    alldata=[]
    #print('scaling_data',scaling_data)
    if plotter:
        plt = plotter
    else:
        from matplotlib import pyplot as plt
    #colors=['r','b','g','y','m','c','sienna','teal','coral','sandybrown','silver','darkseagreen','indigo','maroon']
    colors=cm.gist_ncar(np.linspace(0,1,len(scaling_data.keys())+1))
        #legend=['Elements','Scaling line']
    if alloy:
        name = basename+'.pdf'
        counter=0
        markers={'100':'s','111':'h','110':'D','011':'D','001':'s',
                    '11-20':'^','2-1-10':'v','0001':'>',
                 '010':'P','101': 's','121':'X','221':'*','021':'*'}
        for i,comp in enumerate(scaling_data.keys()):

            compname=get_reduced_alloy_name(comp)
            #if comp in ['Au1Cu1']:
            #    markers=['h','s','D']
            #else:
            allfacetsdata=[]
            for j,facet in enumerate(scaling_data[comp].keys()):
                thisfacetdata=[]
                for term in scaling_data[comp][facet].keys():
                    if descriptors[0] in scaling_data[comp][facet][term].keys():
                     if descriptors[1] in scaling_data[comp][facet][term].keys():
                        #alldata.append([scaling_data[comp][facet][term]['CO'],scaling_data[comp][facet][term]['OCCO']])
                        alldata.append([scaling_data[comp][facet][term][descriptors[0]],scaling_data[comp][facet][term][descriptors[1]]])
                        thisfacetdata.append([scaling_data[comp][facet][term][descriptors[0]],scaling_data[comp][facet][term][descriptors[1]]])

                        #if txt or not legend:
                        if txt:
                            ha=['left','right']
                            if counter%2:
                                plt.text(alldata[-1][0]-0.03,alldata[-1][1],comp+'('+facet+')',
                                    fontweight='bold',color=colors[counter%len(colors)],fontsize=fontsize,ha=ha[counter%2])
                            else:
                                plt.text(alldata[-1][0]+0.03,alldata[-1][1],comp+'('+facet+')',
                                    fontweight='bold',color=colors[counter%len(colors)],fontsize=fontsize,ha=ha[counter%2])
                        #elif isinstance(legend,list) and compname+'('+facet+')' not in legend:
                        #    legend.append(compname+'('+facet+')')

                allfacetsdata.extend(thisfacetdata)
                thisfacetdata=np.array(thisfacetdata)
                if separate_facets:
                    if len(thisfacetdata) > 0:
                        if single_marker_per_alloy:
                            marker='o'
                        elif comp+facet == 'Ni2Sb4001':
                            marker='P'
                        else:
                            marker=markers[facet]
                        plt.plot(thisfacetdata[:,0],thisfacetdata[:,1],
                                        '%s'%marker,
                                        color=colors[counter%len(colors)],
                                        markersize=8,
                                        markeredgecolor='k')
                        if isinstance(legend,list):# and compname+'('+facet+')' not in legend:
                                legend.append(compname+'('+facet+')')
            allfacetsdata=np.array(allfacetsdata)
            if not separate_facets and len(allfacetsdata):
                marker='o'
                colors=cm.gist_ncar(np.linspace(0,1,len(scaling_data.keys())+1))
                plt.plot(allfacetsdata[:,0],allfacetsdata[:,1],
                                        '%s'%marker,
                                        color=colors[counter%len(colors)],
                                        markersize=8,
                                        markeredgecolor='k')
                if isinstance(legend,list):# and compname+'('+facet+')' not in legend:
                                legend.append(compname)

            if len(scaling_data[comp].keys()):
                counter+=1
        alldata=np.array(alldata)
    else:
        name = basename+'%s.pdf'%facet
        for j,elem in enumerate(scaling_data):
            if (descriptors[1] in scaling_data[elem] and
               descriptors[0] in scaling_data[elem]):
                alldata.append([scaling_data[elem][descriptors[0]],scaling_data[elem][descriptors[1]]])

                if different_colors:
                    plt.plot(alldata[-1][0],alldata[-1][1],
                            linestyle="None",
                            marker=markersymbol,color=colors[j%len(colors)],
                            markersize=8,markeredgecolor='k')

                if txt:
                    if basename.split('_')[0] != 'dwf':
                        plt.text(alldata[-1][0]+0.05,alldata[-1][1]+0.05,elem,fontweight='bold',fontsize=fontsize,color='k')
                    else:
                        plt.text(alldata[-1][0]+0.02,alldata[-1][1]+0.02,elem,fontweight='bold',fontsize=fontsize)

                elif isinstance(legend,list) and elem+'('+facet+')' not in legend:
                            legend.append(elem+'('+facet+')')
        alldata=np.array(alldata)
        if not different_colors:
            #print(alldata)
            plt.plot(alldata[:,0],alldata[:,1],pointcolor,markeredgecolor='k')

    plt.xlabel(xname,fontsize=15)
    plt.ylabel(yname,fontsize=15)#'$\Delta$E$_{OCCO}$ [eV]')
    if basename.split('_')[0] != 'dwf':
        #plt.ylim(1.8,5.5)
        plt.ylim(-3,1.)
    if yrange:
        plt.ylim(yrange)
    if xrange:
        plt.xlim(xrange)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #print (alldata)
    if len(alldata) > 1 and fit_line:
        coeffs,dummy = curve_fit(lin_fun,alldata[:,0],alldata[:,1])
        plt.plot([alldata[:,0].min()-0.2,alldata[:,0].max()+0.2],
            lin_fun(np.array([alldata[:,0].min()-0.2,alldata[:,0].max()+0.2]),*coeffs),'k')
        if title is not None:
            plt.title(label=title+', E$_y$= %1.2fE$_x$+%1.2f'%(coeffs[0],coeffs[1]))
        else:
            plt.title(label='E$_{%s}$= %1.2fE$_{%s}$+%1.2f'%(descriptors[0],coeffs[0],descriptors[1],coeffs[1]))
    elif fit_line:
        if len(alldata) < 2:
            print('Scaling line '+name+' has not been plotted due to lack of data')
            return False
    else:
        plt.title(label=title,fontsize=20,fontweight='bold')

    plt.tight_layout()
    if return_plt:
        if legend:
            return plt,legend
        else:
            return plt
    else:
        if legend:
            plt.legend(legend,fancybox=True,bbox_to_anchor=(1,1),fontsize=fontsize)
        plt.tight_layout()
        plt.savefig(home+'/results/'+name)
        plt.close()

def plot_alloy_dictionaries(E_bind_min,name,ylabel='CO binding energy [eV]',
                            zeroline=False,CO_gas_line=False,binding_on_Cu=False,
                            #textshift=0.3, Not used anymore
                            fontsize=14,
                            y_range=[],x_range=[],figsize=None,lattices=None,
                            result_path=os.getcwd(),figax=None,
                            show_markers=False,markercolor='k',facet_colors=True,separators=None,markersize=2,
                            visualize=False,PTMs=[],return_plt=False,alpha=1,
                            constant_shift=0):
    if figax is None:
        fig,ax=plt.subplots()
    else:
        fig,ax=figax

    if figsize is not None:
        from pylab import rcParams
        rcParams['figure.figsize'] = figsize

    xindex=0
    plotted_el,plotted_el_dirname = [],[]

    exceptions=['Ga1Ni1121']
    for i,elem in enumerate(E_bind_min.keys()):
        plotted=False
        for facet in E_bind_min[elem].keys():
          if elem+facet not in ['Ge4Pt2101']:
            print(elem,facet,E_bind_min[elem][facet])
            sites=[]
            for site in E_bind_min[elem][facet].keys():
                if site not in PTMs: sites.append(site)
            for site in E_bind_min[elem][facet].keys():
                if site in PTMs: sites.append(site)
            for site in sites:

                in_range=True
                E = E_bind_min[elem][facet][site]+constant_shift # constant shift is e.g. the free energy correction
                if y_range:
                    if (E < y_range[0] or\
                    E > y_range[1]):
                        in_range=False

                if in_range:
                    if facet_colors:
                        if facet in ['111','112','011','101','11-20','221','021'] or elem+facet in exceptions:
                                    color='orange'
                        elif facet in ['100','001','010','2-1-10','1-102']:
                                    color='darkolivegreen'
                        elif facet in ['0001','110','121']:
                                    color='darkslateblue'
                        else:
                                    color='r'
                    else:
                        color='k'

                    if show_markers:
                        if 'OCCO' in site:
                            ax.plot(xindex,E,color+'s')
                            plotted=True
                        elif 'top' in site:
                            ax.plot(xindex,E,color+'.')
                            plotted=True
                        else:
                            if site in ['Cu','Ag','Au','Pt','Pd','Ni','Ir','Rh']:
                                ax.plot(xindex,E,color=color,marker='D',markeredgecolor='k',markersize=markersize,alpha=alpha)
                                if len(PTMs):
                                    ax.plot(xindex,E,color=color,marker='D',markeredgecolor='k',markersize=markersize,alpha=alpha)
                                    for ptm in PTMs:
                                        if ptm in elem:
                                            break
                                    else:
                                    #    ax.plot(xindex,E,color=color,marker='D',markeredgecolor='k',markersize=markersize)
                                        mcolor='w'
                                        if color in ['orange']: mcolor='k'
                                        ax.annotate(site, xy=(xindex,E),ha='center',va='center',color=mcolor,fontsize=0.5*fontsize)
                            else:
                                ax.plot(xindex,E,color=color,marker='o',markeredgecolor='k',markersize=markersize/2,alpha=alpha)
                            plotted=True
                    else:
                        #print(xindex,E)
                        ax.plot(xindex,E,'w.')
                        #plt.text(xindex-textshift,E,site,color=color,fontsize=fontsize)
                        ax.text(xindex,E,site,color=color,fontsize=fontsize,va='center',ha='center')

                        plotted=True
                    if get_reduced_alloy_name(elem) not in plotted_el:
                            plotted_el.append(get_reduced_alloy_name(elem))
                            plotted_el_dirname.append(elem)
        if plotted:
            xindex+=1
    if separators:
     totsep=0
     seplabels=[]
     if isinstance(separators,dict):
         seplabels=list(separators.keys())
         separators=separators.values()
     for i, sep in enumerate(separators):
         totsep+=sep
         if 0:
             ax.axvline(x=totsep-0.5,linestyle='--',color='k')
         else:
             if not sep: continue
             if color == 'k':
                     color = 'w'
             else:
                     color = 'k'
             ax.axvspan(-0.5+totsep-sep, -0.5+totsep, facecolor=color, alpha=0.1, zorder=-10)
#             plt.annotate(seplabels[i], xy=(-0.5+totsep-sep/2, y_range[1]+0.1), ha='center', va='top', fontsize=0.8*fontsize, color='k', zorder=10)
             if len(seplabels):
                 ax.text(-0.5+totsep-sep/2, y_range[1]+(y_range[1]-y_range[0])*0.1, seplabels[i],
                          ha='center', va='top', fontsize=1.2*fontsize, color='k', zorder=10)

    if zeroline:
        #plt.plot([0,len(plotted_el)],[0,0],'k', linestyle='dashed',linewidth=0.5)
        ax.plot([-0.5,len(plotted_el)],[0,0],'k', linestyle='dashed',linewidth=0.5)
    if CO_gas_line:
        ax.plot([0,len(plotted_el)],[2.5316,2.5316],'k', linestyle='dashed',linewidth=0.5)
        #ax.text(-1,2.5816,'E(CO$_g$)',color='k',fontsize=fontsize)

    if isinstance(binding_on_Cu,float):
        ax.axhline(y=binding_on_Cu+constant_shift,color='k', linestyle='dashed',linewidth=0.5)
#        ax.text(len(plotted_el)-0.2,binding_on_Cu,'E$_b$(Cu)',color='k',fontsize=fontsize)
    elif isinstance(binding_on_Cu,dict):
        for element in binding_on_Cu.keys():
            ax.axhline(y=binding_on_Cu[element]+constant_shift,color='k', linestyle='dashed',linewidth=0.5)
            #plt.plot([-0.5,len(plotted_el)-0.3],[binding_on_Cu,binding_on_Cu],'k', linestyle='dashed',linewidth=0.5)
            #ax.text(0.2,binding_on_Cu[element]+0.01,'E$_b$(%s)'%(element),color='k',fontsize=0.5*fontsize,)

    elif binding_on_Cu:
        ax.plot([-0.5,len(plotted_el)-0.3],[-0.67,-0.67],'k', linestyle='dashed',linewidth=0.5)
        dd
        ax.text(len(plotted_el)-0.2,-0.67,'E$_b$(Cu)',color='k',fontsize=fontsize)

    #plt.xticks(np.arange(len(E_bind.keys())),E_bind.keys(),rotation='vertical')
    if not lattices:
        ax.set_xticks(np.arange(len(plotted_el)),plotted_el,rotation='vertical',fontsize=fontsize)
    else:
        ax.set_xticks(np.arange(len(plotted_el)),plotted_el,rotation='vertical',fontsize=fontsize)
        for latcoli,elem in enumerate(plotted_el_dirname):
            print(latcoli,elem)
            if elem not in lattices:
                print(elem, 'is not in lattices!')
                continue
            if lattices[elem] == 'tetragonal':
                ax.get_xticklabels()[latcoli].set_color('blue')
            if lattices[elem] == 'hexagonal':
                ax.get_xticklabels()[latcoli].set_color('goldenrod')
            if lattices[elem] == 'orthorhombic':
                ax.get_xticklabels()[latcoli].set_color('green')
            if lattices[elem] == 'trigonal':
                ax.get_xticklabels()[latcoli].set_color('red')


    #ax.set_yticks(fontsize=fontsize)
    ax.set_xlim(-0.5,len(plotted_el)-0.3)
    if y_range:
        ax.set_ylim(y_range)
    ax.set_ylabel(ylabel,fontsize=fontsize+10)
    plt.tight_layout()
#    if return_plt:
#        return plt

    if visualize:
        plt.show()
    else:
        plt.savefig(result_path+'/'+name)
        print('Results figure has been written at %s'%(result_path+'/'+name))
    if figax is None:
        plt.close()

def plot_large_alloy_dictionaries(E_bind,name,ylabel='CO binding energy [eV]',zeroline=False,CO_gas_line=False):
    xindex=0
    plotted_el = []
    exceptions=['Cu1Pd1110','Au1Cu1101','Au1Cu1011','Ni1Pt1101','Ni1Pt1101','Au4Cu8Pd4110','Au4Cu8Pd4101','Cu3Pt3001','Ag3Pt3001','Ag3Pd3001','Pd3Pt3001','Ag3Pd3012','Pd3Pt3012','Ag3Pt3012','Cu3Pt3012','Ag3Au3001','Ag3Au3012']
    for i,elem in enumerate(E_bind.keys()):
        plotted=False
        for facet in E_bind[elem].keys():
            for termination in E_bind[elem][facet].keys():
                for site in E_bind[elem][facet][termination].keys():
                   #     print(comp,facet,termination,i,E[e]-E_ref)
                        plt.plot(xindex,E_bind[elem][facet][termination][site]['energy'],'w.')
                        #plt.plot(i,E_bind[elem][facet][site],'w.')
                        if facet in ['111','112','011','101'] or elem+facet in exceptions:
                            plt.text(xindex-0.2,E_bind[elem][facet][termination][site]['energy'],E_bind[elem][facet][termination][site]['site'],color='r')
                        elif facet in ['100','001','010','110']:
                            plt.text(xindex-0.2,E_bind[elem][facet][termination][site]['energy'],E_bind[elem][facet][termination][site]['site'],color='k')
                        else:
                            plt.text(xindex-0.2,E_bind[elem][facet][termination][site]['energy'],E_bind[elem][facet][termination][site]['site'],color='b')
                        plotted=True
                        if elem not in plotted_el:
                            plotted_el.append(elem)
        if plotted:
            xindex+=1

    if zeroline:
        plt.plot([0,len(plotted_el)],[0,0],'k', linestyle='dashed',linewidth=0.5)
    if CO_gas_line:
        plt.plot([0,len(plotted_el)],[2.5316,2.5316],'k', linestyle='dashed',linewidth=0.5)
        plt.text(len(plotted_el),2.5316,'E(CO$_g$)',color='k')

    #plt.xticks(np.arange(len(E_bind.keys())),E_bind.keys(),rotation='vertical')
    plt.xticks(np.arange(len(plotted_el)),plotted_el,rotation='vertical')
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(os.getcwd()+'/results/'+name)

def plot_Ebind_vs_q(home,element,facet,ads,Eb,allcoeffs,slab_wf0=None,x_range=[-30,0],linear_fit=True,functiontxt=True,Capacity=20,zero_volt=4.6):
  #    plt.cla()
      fig = plt.figure()
      ax1 = fig.add_subplot(111)
      colors=np.array(['r','k','b','g','m','y'])
      #if ads=='OCCO':
      #    print(Eb)
      #    das
      for i,site in enumerate(allcoeffs.keys()):
          if max(Eb[site][:,0]) > 0.1:
              x_range = [-30,max(Eb[site][:,0])]
          else:
              x_range = [-30,0]
          #print(site)
          coeffs=allcoeffs[site]
          #print(coeffs)
          if linear_fit:
            ax1.plot(x_range,lin_fun(np.array(x_range),*coeffs),colors[i%len(colors)])
            #plt.title('Eb$_{%s}$($\sigma$) = %1.4f*$\sigma$+%1.4f'%(ads,coeffs[0],coeffs[1]),fontsize=20)
            ax1.set_title('%s(%s)'%(element,facet),fontsize=20)
            if functiontxt:
                ax1.text(x_range[0],x_range[1]*coeffs[0]+coeffs[1],'%s: %1.4f*$\sigma$+%1.4f'%(site.split('_')[-1],coeffs[0],coeffs[1]),fontsize=16,color=colors[i%len(colors)])
            ax1.plot(Eb[site][:,0],Eb[site][:,1],'+r',markersize=8,color=colors[i%len(colors)])
          else:
            ax1.plot(x_range,quad_fun(np.array(x_range),*coeffs),'k')
            #print(coeffs)
            ax1.title('Eb$_{%s}$($\sigma$) = %1.4f*$\sigma^2$+%1.4f*$\sigma$+%1.4f'%(ads,coeffs[0],coeffs[1],coeffs[2]),fontsize=14)


      #plt.title(-20,Eb[-1,1],'Eb = %1.4f*$\sigma$+%1.4f'%(coeffs[0],coeffs[1]),fontsize=15)
      ax1.set_xlabel('Surface charge [$\mu$Ccm$^{-2}$]')
      ax1.set_ylabel('Binding energy [eV]')
      #ax2.plot(range(100),range(100))

      if slab_wf0:
        ax2=ax1.twiny()
        ax2.set_xticks(ax1.get_xticks())
        ax2.set_xbound(ax1.get_xbound())
        ax2.set_xticklabels(np.around(ax1.get_xticks()/Capacity + slab_wf0-zero_volt,1))
        ax2.set_xlabel('U$_{SHE}$ [V]')

      if ads+'_binding_vs_q' not in os.listdir(home+'/results'):
            os.mkdir(home+'/results/'+ads+'_binding_vs_q')
      fig.tight_layout()
      fig.savefig(home+'/results/'+ads+'_binding_vs_q/%s_binding_vs_q_%s%s.pdf'%(ads,element,facet))
      plt.close()

def plot_NEB_band(home,element,facet,allimages,system='OC-CO'):
    from ase.neb import NEBTools
    if isinstance(allimages,dict):
        for charge in allimages.keys():
            if system+'_bands' not in os.listdir(home+'/results/'):
                os.mkdir(home+'/results/'+system+'_bands')
            for image in allimages[charge]:
                image.calc.results['energy']/=2.
                image.calc.results['forces']/=2.
            fig = NEBTools(allimages[charge]).plot_band()
            fig.savefig(home+'/results/'+system+'_bands/'+element+'_'+facet+'_q_'+str(charge)+'.pdf')
            plt.close()

        for charge in allimages.keys():
            fig = NEBTools(allimages[charge]).plot_band()
            fig.savefig(home+'/results/'+system+'_bands/%s_%s_allbands.pdf'%(element,facet))

def plot_NEB_band_local(home):
    from ase.neb import NEBTools
    atoms=[]
    count=0
    for dire in os.listdir(home):
        if dire[0] == '0':
            count+=1
    for dire in range(count):
            atoms.append(read(home+'/0'+str(dire)+'/'+'OUTCAR',index='-1'))
    fig = NEBTools(atoms).plot_band()
    fig.tight_layout()
    fig.savefig('NEB_band.pdf')
    plt.close()

#def plot_front_and_back_barrier(basehome,element,facet,system=ads):
#    if isinstance(allimages,dict):
#        for charge in allimages.keys():



def multicolor_ylabel(ax,list_of_strings,list_of_colors,axis='x',anchorpad=0,**kw):
    """this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    # x-axis label
    if axis=='x' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',**kw))
                    for text,color in zip(list_of_strings,list_of_colors) ]
        xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
        anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,frameon=False,bbox_to_anchor=(0.2, -0.09),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_xbox)

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',rotation=90,**kw))
                     for text,color in zip(list_of_strings[::-1],list_of_colors) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(-0.10, 0.2),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_ybox)

def plot_pourbaix(E_bind,adsorbates=['CO','H','OH'],potentials=np.arange(-1.5,1,0.1)):
    comps = [i for i in E_bind['CO'].keys()]
    E_bind_sorted={}
    dashes=['--','-','.','+']
    colors=['k','b','r']
    #Free energy corrections
    #DeltaG_vib = G_vib(*A)-G_vib(A(gas or l))
    CO_free_en=0.022+0.371
    H_free_en=0.063+0.5*0.035
    OH_free_en=0.33-0.5*0.035+0.009
    #Get Cu reference
    Cu_CO_H_line = E_bind['CO']['Cu']['100']['Cu']+CO_free_en - (E_bind['H']['Cu']['100']['Cu']+H_free_en)
    Cu_CO_OH_line = E_bind['OH']['Cu']['100']['Cu']+OH_free_en - (E_bind['CO']['Cu']['100']['Cu']+CO_free_en)
    Cu_CO_line = E_bind['CO']['Cu']['100']['Cu']+CO_free_en
    #das
    for comp in comps:
        E_bind_sorted[comp]={}
        for facet in E_bind['CO'][comp].keys():
            E_bind_sorted[comp][facet]={}
            plotted=False
            for isite,site in enumerate(E_bind['CO'][comp][facet].keys()):
                E_bind_sorted[comp][facet][site]={}
                #for pot in potentials:
                if comp in E_bind['H'].keys():
                    if facet in E_bind['H'][comp].keys():
                        if site in E_bind['H'][comp][facet].keys():
                         E_bind_sorted[comp][facet][site]['H'] = [E_bind['H'][comp][facet][site]+H_free_en+pot for pot in potentials]
                if comp in E_bind['OH'].keys():
                    if facet in E_bind['OH'][comp].keys():
                        if site in E_bind['OH'][comp][facet].keys():
                            E_bind_sorted[comp][facet][site]['OH'] = [E_bind['OH'][comp][facet][site]+OH_free_en-pot for pot in potentials]

                E_bind_sorted[comp][facet][site]['CO'] = [E_bind['CO'][comp][facet][site]+CO_free_en for i in potentials]

                if 'H' in E_bind_sorted[comp][facet][site].keys() or\
                        'OH' in E_bind_sorted[comp][facet][site].keys():
                 for iads,adsorbate in enumerate(adsorbates):
                    if adsorbate in E_bind_sorted[comp][facet][site].keys():
                        if site in ['Pd','Pt','Ni','Cu']:
                            plt.plot(potentials,E_bind_sorted[comp][facet][site][adsorbate],colors[iads]+'-',label=site+adsorbate)
                        else:
                            plt.plot(potentials,E_bind_sorted[comp][facet][site][adsorbate],colors[iads]+'--',label=site+adsorbate)
                 plotted=True
            #try:
            #    print(comp, facet, E_bind['OH'][comp][facet])
            #    print(comp, facet, E_bind['H'][comp][facet])
            #    print(comp, facet, E_bind['CO'][comp][facet])
            #except:
            #    pass

            if plotted:
                plt.legend()
                plt.ylabel('Binding free energy [eV]')
                plt.xlabel('Voltage vs RHE [V]')
                plt.axvline(x=Cu_CO_H_line,color='b',linestyle=':')
                plt.axvline(x=Cu_CO_OH_line,color='r',linestyle=':')
                plt.axhline(y=Cu_CO_line,color='k',linestyle=':')
                plt.xlim(potentials[0],potentials[-1])
                plt.ylim(-1.5,0.5)
                plt.title(get_reduced_alloy_name(comp)+'('+facet+')')
                plt.tight_layout()

                plt.savefig('results/pourbaix/%s_%s.pdf'%(comp,facet))
                plt.close()




def sort_alloys_for_plot(msortlist,data,PTM_all):
   #Sort compositions by the post transition metals
    compsorted={}
    TM='TMs'
    compsorted[TM]=[]
    for comp in data.keys():
        for ptm in PTM_all:
            if  ptm in comp:
                break
        else:
            compsorted[TM].append(comp)

    for tm  in msortlist:
        compsorted[tm]=[i for i in data.keys()
            if (tm in i and i not in compsorted[TM])]
    nperM={tm: len(compsorted[tm]) for tm in compsorted.keys()}
    return compsorted,nperM
