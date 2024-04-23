def main():
    def read_configs():
        def read_conf_file(file_name):
            cmd = {}
            ifs = open(file_name, 'r')
            a_line = ifs.readline()
            while a_line:
                if a_line[0] != '#':
                    line = a_line.strip('\n').split('\t')
                    if len(line) == 2:
                        cmd[line[0]] = line[1]
                    else:
                        cmd[line[0]] = line[1:]
                a_line = ifs.readline()
            ifs.close()
            return(cmd)

        import sys
        cmd = {'error': False}
        if len(sys.argv) > 2:
            print('Error:\t too many input: ' + str(len(sys.argv)) + ' (> 2)')
            print('\tread_config in ./src/py/MyPackage')
            print()
            cmd['error'] = True
        elif len(sys.argv) == 1:
            print('Error:\tan input file path is required')
            print('\tread_config in ./src/py/MyPackage')
            print()
            cmd['error'] = True
        else:
            confs = []
            ifs = open(sys.argv[1], 'r')
            # print('\nr: '+sys.argv[1])
            a_line = ifs.readline()
            while a_line:
                if a_line[0] != '#':
                    line = a_line.strip('\n').split('\t')
                    if len(line) == 2:
                        cmd[line[0]] = line[1]
                        if line[1][-5:] == '.conf':
                            confs.append(line[0])
                    else:
                        cmd[line[0]] = line[1:]
                a_line = ifs.readline()
            for conf in confs:
                # print('  '+conf+'\tr: '+cmd[conf])
                cmd[conf] = read_conf_file(cmd[conf])
                # for key in cmd[conf].keys():
                #     print('\t', key, cmd[conf][key])
        return(cmd)
    cmd = read_configs()
    #
    if 'ofn' not in cmd:
        cmd['ofn'] = cmd['ifn']+'.pdf'
    X = read_trj(cmd['ifn'])
    write_fig(X, cmd)
    return




def read_trj(ifn: str):
    ifs = open(ifn, 'r')
    #
    X = []
    a_line = ifs.readline()
    size = len(a_line.split('\t'))
    for itr in range(size):
        X.append([])
    while a_line:
        line = a_line.strip('\n').split('\t')
        for itr in range(len(line)):
            X[itr].append(float(line[itr]))
        a_line = ifs.readline()
    #
    ifs.close()
    print('r: '+ifn)
    return(X)


def write_fig(X: list, cmd: dict):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy.lib.scimath as scimath
    # import scipy as sp
    plt.rcParams.update({
        'font.size': 16,
        'axes.grid': True,
        'text.usetex': True,
        'font.serif': ['Computer Modern Roman'],
        'text.latex.preamble': r"\usepackage{bm}"
        })  # noqa W605
    fig, ax = plt.subplots(1, 1, dpi=300)
    xlim = [float(x) for x in cmd['lim'][0].split(' ')]
    ylim = [float(x) for x in cmd['lim'][1].split(' ')]
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(cmd['xlabel'])
    ax.set_ylabel(cmd['ylabel'])  # noqa W605
    #
    if type(cmd['vline.val']) is str:
        cmd['vline.val'] = [cmd['vline.val']]
    if type(cmd['vline.lbl']) is str:
        cmd['vline.lbl'] = [cmd['vline.lbl']]
    for itr in range(len(cmd['vline.val'])):
        ax.axvline(float(cmd['vline.val'][itr]), color='g', alpha=.3)
        ax.text(float(cmd['vline.val'][itr]), ylim[1], cmd['vline.lbl'][itr],
                ha='center', va='bottom',
                color=colors.hsv_to_rgb((1/3,1,.5)))  # noqa W605
    ax.axhline(scimath.log(1.1), color='g', alpha=.3)
    ax.text(xlim[1], scimath.log(1.1), '$\\ln{(1.1)}$',
            ha='left', va='center',
            rotation=270,
            color=colors.hsv_to_rgb((1/3,1,.5)))  # noqa W605
    for itr in range(len(X[1:])):
        s_col = tuple([float(x) for x in cmd['colors'][itr].split(' ')])
        if len(s_col) == 4:
            col = colors.hsv_to_rgb(s_col[:3])
            ax.loglog(X[0], X[itr+1],
                      label=cmd['labels'][itr], zorder=float(cmd['zorders'][itr]),  # noqa W605
                      c=col, alpha=s_col[3],
                      ls=cmd['linestyles'][itr], lw=cmd['linewidths'][itr])  # noqa W605
        elif len(s_col) == 3:
            col = colors.hsv_to_rgb(s_col)
            ax.loglog(X[0], X[itr+1],
                      label=cmd['labels'][itr], zorder=float(cmd['zorders'][itr]),  # noqa W605
                      c=col,
                      ls=cmd['linestyles'][itr], lw=cmd['linewidths'][itr])  # noqa W605
        else:
            print('error: ', itr ,'th line color-format is inappropriate: ', s_col)  # noqa W605
            quit
    if 'legend' in cmd:
        if 'bbox_to_anchor' in cmd['legend']:
            ax.legend(loc=cmd['legend']['loc'],
                      fontsize=cmd['legend']['fontsize'],
                      title=cmd['legend']['title'],
                      ncol=int(cmd['legend']['ncol']),
                      bbox_to_anchor=tuple([float(x) for x in cmd['legend']['bbox_to_anchor']]))
        else:
            ax.legend(loc=cmd['legend']['loc'],
                      title=cmd['legend']['title'],
                      fontsize=cmd['legend']['fontsize'])
    else:
        ax.legend(fontsize='small', title=cmd['legend.title'], loc=cmd['legend.loc'])  # noqa W605
    ax.grid(True, which="both", ls="-", color='0.9')
    fig.savefig(cmd['ofn'], bbox_inches='tight')
    print('w: '+cmd['ofn'])
    return


main()
