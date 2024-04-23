def main():
    import numpy as np
    from read import read
    R = 8.31446261815324
    #
    cmd = dict(ofn="./data/AVE/CCSDT/assets/Arrhenius/Arrhenius.evals.pdf")
    cmd['rcParams'] = {
            'font.size': 16,
            'text.usetex': True,
            'text.latex.preamble': r"\usepackage{bm}"
            }
    cmd['figsize'] = (6, 6)
    cmd['lim'] = [[0.22, 0.45], [5e9, 7e11]]
    cmd['yscale'] = 'log'
    cmd['xscale'] = ''
    cmd['legend'] = {
            'title': '',
            'loc': 'upper left', #'lower left', # 'best',
            'fontsize': 'x-small', # 'medium'
            'bbox_to_anchor': (0, -.15), # (1.1, 1.02),# (.5, -.15),
            'ncol': 3
            }
    cmd['xlabel'] = '${10}^{-3}/R T~/~[\\mathrm{mol}/\\mathrm{kJ}]$'
    cmd['ylabel'] = '$-\\lambda_k^{(\hat{n})}~/~[1/\\mathrm{s}]$'
    # cmd['fit_lrange'] = {keys[0]: [0,5], keys[1]: [1,5]}
    cmd['linestyle'] = '-' # None
    #
    lines = {
            "$-\\lambda^{(3) \\mathrm{ID-EL}}_2$":  dict(line=False, color=".90\t.8\t.8", marker="2", zorder=21, alpha=.5, ifn="./data/AVE/CCSDT/assets/HierarchicalLumping/EVals/id3/lmp.eval2.nb.tsv"),
            "$-\\lambda^{(4) \\mathrm{ID-EL}}_3$":  dict(line=False, color=".56\t.8\t.8", marker="2", zorder=31, alpha=.5, ifn="./data/AVE/CCSDT/assets/HierarchicalLumping/EVals/id4/lmp.eval3.nb.tsv"),
            "$-\\lambda^{(5) \\mathrm{ID-EL}}_4$":  dict(line=False, color=".23\t.8\t.8", marker="2", zorder=41, alpha=.5, ifn="./data/AVE/CCSDT/assets/HierarchicalLumping/EVals/id5/lmp.eval4.nb.tsv"),
            "$-\\lambda^{(3) \\mathrm{PCCA}}_2$":   dict(line=False, color=".10\t.8\t.8", marker="1", zorder=22, alpha=.5, ifn="./data/AVE/CCSDT/assets/HierarchicalLumping/EVals/ss3/gpe.eval2.nb.tsv"),
            "$-\\lambda^{(4) \\mathrm{PCCA}}_3$":   dict(line=False, color=".76\t.8\t.8", marker="1", zorder=32, alpha=.5, ifn="./data/AVE/CCSDT/assets/HierarchicalLumping/EVals/ss4/gpe.eval3.nb.tsv"),
            "$-\\lambda^{(5) \\mathrm{PCCA}}_4$":   dict(line=False, color=".43\t.8\t.8", marker="1", zorder=42, alpha=.5, ifn="./data/AVE/CCSDT/assets/HierarchicalLumping/EVals/ss5/gpe.eval4.nb.tsv"),
            "$-\\lambda_2$":                        dict(line=True,  color="1.0\t.8\t.8", marker="+", zorder=20, alpha=1., ifn="./data/AVE/CCSDT/assets/HierarchicalLumping/EVals/full/eval2.nb.tsv"),
            "$-\\lambda_3$":                        dict(line=True,  color=".66\t.8\t.8", marker="+", zorder=30, alpha=1., ifn="./data/AVE/CCSDT/assets/HierarchicalLumping/EVals/full/eval3.nb.tsv"),
            "$-\\lambda_4$":                        dict(line=True,  color=".33\t.8\t.8", marker="+", zorder=40, alpha=1., ifn="./data/AVE/CCSDT/assets/HierarchicalLumping/EVals/full/eval4.nb.tsv")
            }
    for conf in lines[list(lines.keys())[0]].keys():
        cmd[conf+'s'] = []
    for key in lines.keys():
        for conf in lines[key].keys():
            cmd[conf+'s'].append(lines[key][conf])
    #
    X, Y = {}, {}
    for key in lines.keys():
        XY = read.evals(lines[key]['ifn'])
        keys_XY = list(XY.keys())
        X[key] = +1 / (R * (273.15 + np.array(XY[keys_XY[0]])))
        Y[key] = -1 * np.array(XY[keys_XY[2]])
    #
    write_fig_eval(X, Y, lines.keys(), cmd)
    #
    cmd['insets'] = {
            'A': [[0.250, 0.258], [1.40e11,1.65e11]],
            'B': [[0.368, 0.376], [1.19e11,1.43e11]]
            }
    cmd['yMtick'] = [1.2e11, 1.4e11, 1.6e11, 1.8e11]
    cmd['ymtick'] = [1.3e11, 1.5e11, 1.7e11, 1.9e11]
    cmd['rcParams']['font.size'] = 10
    for key in cmd['insets'].keys():
        write_fig_eval_inset(X, Y, lines.keys(), cmd, cmd['insets'][key], key)
    return


def write_fig_eval(X: dict, Y: dict, keys: list, cmd: dict):
    def flip(items, ncol):
        return itertools.chain(*[items[i::ncol] for i in range(ncol)])

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import math
    import matplotlib.patches as patches
    import itertools
    scl = 1e3
    plt.rcParams.update(cmd['rcParams'])
    fig = plt.figure(figsize=cmd['figsize'], dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    if cmd['yscale'] == 'log':
        ax.set_yscale('log')
        minor_tick = []
        major_tick = []
        yrmin, yrmax = math.floor(np.log10(cmd['lim'][1][0])), math.ceil(np.log10(cmd['lim'][1][1]))
        for p in range(yrmin, yrmax, 1):
            minor_tick += [x*10**p for x in range(10)]
            major_tick += [10**p]
    ax.set_yticks(minor_tick, minor=True)
    ax.set_yticks(major_tick)
    ax.set_xlabel(cmd['xlabel'])
    ax.set_ylabel(cmd['ylabel'])
    ax.set_xlim(cmd['lim'][0][0], cmd['lim'][0][1])
    ax.set_ylim(cmd['lim'][1][0], cmd['lim'][1][1])
    ax.grid(True, which="both", ls="-", color='0.9')
    #
    if 'insets' in cmd:
        for key in cmd['insets'].keys():
            inset = cmd['insets'][key]
            ax.add_patch(patches.Rectangle((inset[0][0], inset[1][0]),
                                            inset[0][1]-inset[0][0],
                                            inset[1][1]-inset[1][0],
                                            ls="-", ec="k", fc="none", zorder=100))
            ax.text(sum(inset[0])/2, sum(inset[1])/2, key, zorder=100, alpha=.5,
                    ha='center', va='center_baseline',
                    fontfamily='sans-serif', fontsize='large')
    #
    itr = -1
    ofn = cmd['ofn']+'.Arrhenius.tsv'
    ofs = open(ofn, 'w')
    ofs.write('\t'.join(['#', 'type', 'A', 'Ea'])+'\n')
    for key in keys:
        X[key] = np.array(X[key])
        itr += 1
        s_col = tuple([float(x) for x in cmd['colors'][itr].split('\t')])
        col = colors.hsv_to_rgb(s_col)
        if 'fit_lrange' in cmd:
            Ea, ST = np.polyfit(scl*X[key][cmd['fit_lrange'][key][0]:cmd['fit_lrange'][key][1]], np.log(Y[key][cmd['fit_lrange'][key][0]:cmd['fit_lrange'][key][1]]), 1)
        else:
            Ea, ST = np.polyfit(scl*X[key], np.log(Y[key]), 1)
        ofs.write('\t'.join([str(x) for x in [key, np.exp(ST), Ea*scl]])+'\n')
        if cmd['lines'][itr]:
            # ax.plot(scl*X[key], np.exp(np.multiply(Ea, scl*X[key])+ST),
            ax.plot(scl*X[key], Y[key],
                    label='',
                    color=col,
                    alpha=0.2,
                    zorder=cmd['zorders'][itr],
                    lw=1.,
                    linestyle=cmd['linestyle'])
        ax.scatter(scl*X[key], Y[key],
                label=key,
                color=col,
                s=150,
                alpha=cmd['alphas'][itr],
                zorder=cmd['zorders'][itr],
                marker=cmd['markers'][itr],
                edgecolor=None
                )
    ofs.close()
    print('w: '+ofn)
    if 'legend' in cmd:
        if 'bbox_to_anchor' in cmd['legend']:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(
                    flip(handles, int(cmd['legend']['ncol'])), flip(labels, int(cmd['legend']['ncol'])),
                    title=cmd['legend']['title'],
                    loc=cmd['legend']['loc'],
                    fontsize=cmd['legend']['fontsize'],
                    ncol=int(cmd['legend']['ncol']),
                    bbox_to_anchor=tuple([float(x) for x in cmd['legend']['bbox_to_anchor']])
                    )
        else:
            ax.legend(
                    title=cmd['legend']['title'],
                    loc=cmd['legend']['loc'],
                    fontsize=cmd['legend']['fontsize']
                    )
    fig.savefig(cmd['ofn'], bbox_inches='tight')
    print(cmd['ofn'])
    return


def write_fig_eval_inset(X: dict, Y: dict, keys: list, cmd: dict, lim: list, title: str):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import math
    import matplotlib.patches as patches
    scl = 1e3
    plt.rcParams.update(cmd['rcParams'])
    fig = plt.figure(figsize=(cmd['figsize'][0]*.25,cmd['figsize'][1]*.25), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    if cmd['yscale'] == 'log':
        ax.set_yscale('log')
    ax.yaxis.tick_right()
    ax.set_yticks([1.2e11, 1.4e11, 1.6e11, 1.8e11])
    ax.set_yticks([1.1e11, 1.3e11, 1.5e11, 1.7e11, 1.9e11], minor=True)
    ax.set_yticklabels(['', '', '', '', ''], minor=True)
    ax.set_xlim(lim[0][0], lim[0][1])
    ax.set_ylim(lim[1][0], lim[1][1])
    ax.grid(True, which="both", ls="-", color='0.9')
    #
    itr = -1
    for key in keys:
        X[key] = np.array(X[key])
        itr += 1
        s_col = tuple([float(x) for x in cmd['colors'][itr].split('\t')])
        col = colors.hsv_to_rgb(s_col)
        if cmd['lines'][itr]:
            ax.plot(scl*X[key], Y[key],
                    label='',
                    color=col,
                    alpha=0.2,
                    zorder=cmd['zorders'][itr],
                    lw=1.,
                    linestyle=cmd['linestyle'])
        ax.scatter(scl*X[key], Y[key],
                label=key,
                color=col,
                s=150,
                alpha=cmd['alphas'][itr],
                zorder=cmd['zorders'][itr],
                marker=cmd['markers'][itr])
    fig.savefig(cmd['ofn']+'.'+title+'.inset.pdf', bbox_inches='tight')
    print(cmd['ofn']+'.'+title+'.inset.pdf')
    return


main()
