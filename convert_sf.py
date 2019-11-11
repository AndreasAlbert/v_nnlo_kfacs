#!/usr/bin/env python

from matplotlib import pyplot as plt
import numpy as np
import uproot
import matplotlib
from scipy import signal
font = {'family' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

def read_file(path_to_file):
    state = 'waiting'
    x, y = {}, {}
    with open(path_to_file) as f:
        for l in f.readlines():
            if state == 'waiting':
                if 'BEGIN HISTO' in l:
                    name = l.split()[-1]
                    if 'PDF' in name:
                        continue
                    state = 'histo'
                    x[name] = []
                    y[name] = []
                else:
                    continue
            elif state == 'histo':
                if 'END' in l:
                    state = 'waiting'
                    name = None
                elif l.startswith('#'):
                    continue
                else:
                    parts = l.split()
                    x[name].append((float(parts[0]), float(parts[1])))
                    y[name].append(float(parts[2]))
    for k, v in x.items():
        x[k] = np.array(v)
    for k, v in y.items():
        y[k] = np.array(v)
    return x,y 


def plot_sf(tag):
    x, y = read_file(f'input/{tag}.dat')

    sf_y = y[f'{tag}_pTV_NNLO'] / y[f'{tag}_pTV_NLO']

    ix = 0.5 * (x[f'{tag}_pTV_LO'][:,0] + x[f'{tag}_pTV_LO'][:,1])

    plt.gcf().clf()
    plt.plot(ix, sf_y,'-o',label=tag)
    
    
    
    fit = False
    if fit:
        degree = 6
        pars = np.polyfit(np.log(ix), sf_y, degree)
        fit = np.zeros(len(ix))
        for n, p in enumerate(reversed(pars)):
            print(n,p)
            fit = fit + p * (np.log(ix)**n)
        
        print(fit)
        plt.plot(ix, fit,'o-',label=f'Fit poly {degree}')

    flit = False;
    if flit:
        filt=signal.savgol_filter(sf_y,
                           5, # window size used for filtering
                           2), # order of fitted polyn
        
        plt.plot(ix, filt[0],'-',linewidth=2, label='Filtered')                
    tmp = set(x[f'{tag}_pTV_LO'][:,0])
    tmp = tmp.union(set(x[f'{tag}_pTV_LO'][:,1]))
    
    sf_x = np.array(sorted(list(tmp)))
    outfile[tag] = (sf_y,sf_x)

    plt.gca().set_xlim([80,4e3])
    plt.gca().set_ylim([1,1.15])

    plt.legend()
    plt.xscale('log')
    plt.xlabel('$p_{T}$ (V)')
    plt.ylabel('QCD NLO -> NNLO SF')
    plt.gcf().savefig(f'output/{tag}.pdf')


names = {
    'eej' : 'DY',
    'vvj' : r'Z->$\nu\nu$',
    'evj' : r'W$^{\pm}$',
    'aj' : r'$\gamma$',
}
def overlay_sf():

    for tag in ['eej', 'evj','vvj','aj']:
        x, y = read_file(f'input/{tag}.dat')

        sf_y = y[f'{tag}_pTV_NNLO'] / y[f'{tag}_pTV_NLO']

        ix = 0.5 * (x[f'{tag}_pTV_LO'][:,0] + x[f'{tag}_pTV_LO'][:,1])

        plt.plot(ix, sf_y,'-o',label=names[tag])

    plt.gca().set_xlim([80,4e3])
    plt.gca().set_ylim([1,1.15])

    plt.legend()
    plt.xscale('log')
    plt.xlabel('Boson $p_{T}$ (GeV)')
    plt.ylabel('QCD NLO -> NNLO SF')
    plt.gcf().savefig(f'output/qcd_nnlo_all.pdf')



# outfile = uproot.recreate(f'lindert_qcd_nnlo_sf.root')
# for tag in ['eej', 'evj','vvj','aj']:
#     plot_sf(tag, outfile)

overlay_sf()
