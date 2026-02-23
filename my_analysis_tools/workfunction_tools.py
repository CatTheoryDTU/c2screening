import numpy as np
import sys,os



def get_workfunction(tmplist,facet,system='slab',charges=0.00,site=None,code='VASP',verbose=False,barr_index=None):
    if isinstance(charges,(float,int)):
        charges=[charges]
    if code=='QE':
        if system == 'slab':
                if facet in ['100','111']:
                #if 1:
                    filename = 'slab_%1.2f/log'%charge
                elif facet in ['211']:
                    filename = 'q_%1.2f/log'%charge
        elif system == 'OCCO':
                    filename = '%s/q_%1.2f/log'%(site,charge)

        inlines = open(filename).readlines()
        #print('inlines',filename)
        for line in inlines:
                    if line.split()[1:4] == ['Fermi', 'energy', 'is']:
                        efermi = float(line.split()[-2])
                    elif line.split()[1:4] == ['Fermi', 'energy', 'shift']:
                        shift = float(line.split()[-2])

        slab_wf = -(efermi+shift)
        print( system+' workfunction at q=%1.2f:'%(charge),slab_wf)
        return slab_wf

    elif code=='VASP':
        Phi,Phi_shift={},{}
        for charge in charges:
        #filebase = 'q_%1.2f/OUTCAR'%charge
            Ef_uncorr=None
            Ef_shift=None

            if site:
                filebase=site+'/q_%1.2f/'%charge
            else:
                filebase='q_%1.2f/'%charge

            if barr_index is None:
                try:
                    inlines = open(filebase+'OUTCAR').readlines()
                except:
                    continue
            else:
                try:
                    inlines = open(filebase+'%02d/OUTCAR'%(barr_index)).readlines()
                    #print(os.getcwd())
                    #inlines = open('%02d/OUTCAR'%(barr_index)).readlines()
                except:
                #    print('here',charge)
                #    print(filebase,os.listdir(os.getcwd()),'%02d/OUTCAR'%(barr_index))
                    continue

            for line in inlines:
                try:
                    if line.split()[0] == 'E-fermi':
                        Ef_uncorr = float(line.split()[2])
                except:
                    pass
            if barr_index is None:
                inlines = open(filebase+'vasp.out').readlines()
            else:
                try:
                    inlines = open(filebase+'01/vasp.out').readlines()
                except FileNotFoundError:
                    inlines = open(filebase+'vasp.out').readlines()

            for line in inlines:
                try:
                    if line.split()[0] == 'FERMI_SHIFT' and float(line.split()[2]) > 1e-5:
                        Ef_shift = float(line.split()[2])
                except:
                    pass
            if Ef_uncorr is not None and Ef_shift is not None:
                Phi[np.around(charge,2)]=-(Ef_uncorr+Ef_shift)
                Phi_shift[np.around(charge,2)] = -Ef_shift

            if verbose:
                print(system+' potential at q=%1.2f: %1.5f'%(charge,-(Ef_uncorr+Ef_shift)))
        #return -(Ef_uncorr+Ef_shift),-Ef_shift
        return Phi,Phi_shift

def get_workfunction_GPAW(tmplist,facet,system='slab',charge=0.00,site=None):
    print(os.getcwd())
    if system == 'slab':
                filename = 'q_%1.2f.txt'%charge
    elif system == 'OCCO':
                filename = '%s/q_%1.2f.txt'%(site,charge)

    inlines = open(filename).readlines()
    for line in inlines:
                try:
                    if line[:12] == 'Dipole-layer':#.split('-')[:2] == ['Dipole','layer']:
                        #print (line)
                        wf = float(line.split()[-2])
                except IndexError:
                    pass

    print( system+' workfunction at q=%1.2f:'%(charge),wf)
    return wf
