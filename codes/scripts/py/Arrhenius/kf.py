def main():
    from read import read
    import numpy as np
    R = 8.31446261815324
    Tk0 = 273.15
    ifn_pe  = "./data/AVE/CCSDT/assets/RateConstants/kFs/pe.kf.nb.tsv"
    ifn_rl  = "./data/AVE/CCSDT/assets/RateConstants/kFs/rl.kf.nb.tsv"
    ifn_gpe = "./data/AVE/CCSDT/assets/HierarchicalLumping/EVals/id2/gpe.ks.nb.tsv"
    ifn_lmp = "./data/AVE/CCSDT/assets/HierarchicalLumping/EVals/id2/lmp.ks.nb.tsv"
    ifn_exp = "./data/AVE/CCSDT/source/SchulerMurphy.tsv"
    ofn = './data/AVE/CCSDT/assets/Arrhenius/Arrhenius.kf'
    #
    labels = [
            'rate limiting',
            'pre-equilibrium',
            'generalized pre-equilibrium',
            'lumping (this study)',
            'experiment',
            ]
    X, Y = {}, {}
    for label in labels:
        X[label] = []
        Y[label] = []
    #
    # rate limiting
    label = labels[0]
    XY = read.k_fwd(ifn_rl)
    for Tc in XY.keys():
        X[label].append(1 / (R * (Tk0 + float(Tc))))
        Y[label].append(XY[Tc])
    # pre-equilibrium
    label = labels[1]
    XY = read.k_fwd(ifn_pe)
    for Tc in XY.keys():
        X[label].append(1 / (R * (Tk0 + float(Tc))))
        Y[label].append(XY[Tc])
    #
    # generalized pre-equilibrium
    label = labels[2]
    rm = read.rate_matrix(ifn_gpe)
    key = (1,2)
    for Tc in rm.keys():
        X[label].append(1 / (R * (Tk0 + float(Tc))))
        Y[label].append(rm[Tc][key])
    # lumping
    label = labels[3]
    rm = read.rate_matrix(ifn_gpe)
    key = (1,2)
    for Tc in rm.keys():
        X[label].append(1 / (R * (Tk0 + float(Tc))))
        Y[label].append(rm[Tc][key])
    #
    # experiment
    X['experiment'], Y['experiment'] = read.kf_exp(ifn_exp)
    #
    cmd = {}
    cmd['ofn'] = ofn
    cmd['lim'] = [[0.2, 0.45], [1e-13, 1]]
    cmd['colors'] = ['.6\t1\t1', '.8\t1\t1', '1\t1\t1', '1\t1\t0', '1\t0\t.5']
    cmd['alpha'] = ['1', '1', '1', '1', '1']
    cmd['marker'] = ['1', '2', '+', 'x', '.']
    cmd['labels'] = labels
    write_fig(X, Y, cmd)
    return


def write_fig(X: dict, Y: dict, cmd: dict):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    scl = 1e3
    plt.rcParams.update({
        'font.size': 16,
        'text.usetex': True,
        'font.serif': ['Computer Modern Roman'],
        'text.latex.preamble': r"\usepackage{bm}"
        })
    fig = plt.figure(figsize=(6, 6), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('${10}^{-3}/R T~/~[\\mathrm{mol}/\\mathrm{kJ}]$')
    ax.set_ylabel('$k_{\mathrm{P,R}}^{(2)}~/~[1/\\mathrm{s}]$')
    ax.set_yscale('log')
    itr = -1
    ofs = open(cmd['ofn']+'.tsv', 'w')
    ofs.write('\t'.join(['#', 'type', 'A', 'Ea'])+'\n')
    for key in cmd['labels']:
        X[key] = np.array(X[key])
        itr += 1
        s_col = tuple([float(x) for x in cmd['colors'][itr].split('\t')])
        col = colors.hsv_to_rgb(s_col)
        Ea, ST = np.polyfit(scl*X[key], np.log(Y[key]), 1)
        ofs.write('\t'.join([str(x) for x in [key, np.exp(ST), Ea*scl]])+'\n')
        ax.plot(scl*X[key], np.exp(np.multiply(Ea, scl*X[key])+ST),
                color=col, alpha=0.2, zorder=10)
        ax.scatter(scl*X[key], Y[key],
                   label=key, color=col,
                   marker=cmd['marker'][itr],
                   alpha = float(cmd['alpha'][itr]),
                   linewidth=1, s=100, zorder=11)
    ofs.close()
    print('w: '+cmd['ofn']+'.tsv')
    ax.legend(loc='lower left', fontsize='small', title=r'')
    minor_tick = []
    major_tick = []
    for p in range(-12, 0, 2):
        major_tick += [10**p]
    for p in range(-13, 0, 1):
        minor_tick += [x*10**p for x in range(10)]
    ax.set_yticks(minor_tick, minor=True)
    ax.set_yticks(major_tick)
    ax.set_xlim(cmd['lim'][0][0], cmd['lim'][0][1])
    ax.set_ylim(cmd['lim'][1][0], cmd['lim'][1][1])
    ax.grid(True, which="both", ls="-", color='0.9')
    fig.savefig(cmd['ofn']+'.pdf', bbox_inches='tight')
    print('w:', cmd['ofn']+'.pdf')
    return


main()
