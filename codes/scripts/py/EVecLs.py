def main():
    Tcs = list(range(0,250+1, 50))
    for Tc in Tcs:
        ifn_lbl = "./data/AVE/CCSDT/assets/EL.MTH/T"+str(Tc).zfill(3)+"/full/labels.tsv"
        lbls = get_lbls(ifn_lbl)
        ifn_EVecL = "./data/AVE/CCSDT/assets/HierarchicalLumping/T"+str(Tc).zfill(3)+"/ESys/EVecL.nb.tsv"
        EVecLs = get_EVecLs(ifn_EVecL, lbls)
        fig_EVecLs("./data/AVE/CCSDT/assets/EVecLs/T"+str(Tc).zfill(3)+'.pdf', EVecLs, 4, Tc)
    return


def get_lbls(ifn: str):
    lbls = {}
    print('r: '+ifn)
    ifs = open(ifn, 'r')
    a_line = ifs.readline()
    while a_line:
        if a_line[0] != '#':
            line = a_line.rstrip('\n').split('\t')
            lbls[line[0]] = line[1]
        a_line = ifs.readline()
    ifs.close()
    return(lbls)


def get_EVecLs(ifn: str, lbls: dict):
    EVecLs = []
    print('r: '+ifn)
    ifs = open(ifn, 'r')
    a_line = ifs.readline()
    while a_line:
        if a_line[0] != '#':
            line = a_line.rstrip('\n').split('\t')
            EVecL = {}
            for itr in range(len(line)):
                EVecL[int(lbls[str(itr)])] = float(line[itr])
            EVecLs.append(EVecL)
        a_line = ifs.readline()
    ifs.close()
    return(list(reversed(EVecLs)))


def fig_EVecLs(ofn: str, EVecLs: list, nhat: int, Tc: int):
    import matplotlib.pyplot as plt
    import numpy as np
    cmd = {}
    cmd['rcParams'] = {
        'font.size': 10.0, # 10.0
        'text.usetex': True, # False
        'text.latex.preamble': r"\usepackage{braket}",
        'font.serif': ['Computer Modern Roman', 'Times New Roman'],
        'font.family': ['serif'], # ['sans-serif']
        'lines.linewidth': 1.0, # 1.5
        'figure.dpi': 300.0, # 100.0
        'figure.figsize': [6.0, 3.0], # [6.4, 4.8]
        'axes.grid': False,
        'axes.grid.axis': 'both',
        'axes.grid.which': 'major',
        'grid.color': '0.9'
        }
    plt.rcParams.update(cmd['rcParams'])  # noqa W605
    fig, ax = plt.subplots(1,1)
    X = []
    for itr in range(1, 23+1):
        X.append(np.nan)
        X.append(itr)
        X.append(itr)
    for k in range(1,nhat+1):
        EVecL = []
        EVecL_tld = []
        for itr in range(1, 23+1):
            EVecL.append(np.nan)
            EVecL.append(k*2-1)
            EVecL.append(k*2-1+EVecLs[k][itr])
            #
            EVecL_tld.append(np.nan)
            EVecL_tld.append(k*2-1)
            if EVecLs[k][itr] > 0:
                EVecL_tld.append(k*2-1+EVecLs[k][itr]/max(EVecLs[k].values()))
            else:
                EVecL_tld.append(k*2-1-EVecLs[k][itr]/min(EVecLs[k].values()))
            #
            err = 2e-1
            if EVecLs[k][itr]>0:
                if EVecLs[k][itr]/max(EVecLs[k].values()) < err:
                    rect = plt.Rectangle((itr-.5,k*2-2),1,2,fc='mistyrose')
                else:
                    rect = plt.Rectangle((itr-.5,k*2-2),1,2,fc='salmon')
            else:
                if EVecLs[k][itr]/min(EVecLs[k].values()) < err:
                    rect = plt.Rectangle((itr-.5,k*2-2),1,2,fc='lavender')
                else:
                    rect = plt.Rectangle((itr-.5,k*2-2),1,2,fc='skyblue')
            ax.add_patch(rect)
        ax.plot([.5, 23.5], [k*2-1, k*2-1], color='gray')
        ax.plot(X, np.array(EVecL_tld), lw=10, solid_capstyle='butt',color='k')
        ax.plot(X, np.array(EVecL), lw=2, solid_capstyle='butt',color='darkgray')
    #
    ax.set_xticks(range(1,24))
    ax.set_yticks(sorted(list(range(1,2*nhat,2))+[x+.5 for x in range(2*nhat)]))
    ax.set_yticklabels(['-0.5', '$\\bra{1}$ 0.0', '0.5',
                        '-0.5', '$\\bra{2}$ 0.0', '0.5',
                        '-0.5', '$\\bra{3}$ 0.0', '0.5',
                        '-0.5', '$\\bra{4}$ 0.0', '0.5'])
    ax.set_xlim([.5, 23.5])
    ax.set_ylim([0,8])
    ax.set_ylabel('$'+str(Tc)+'^\circ \mathrm{C}$')
    #
    fig.savefig(ofn, bbox_inches='tight')
    print('w:', ofn)
    return()

main()
