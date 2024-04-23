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
    if 'ofn' not in cmd:
        cmd['ofn'] = cmd['ifn'] + '.pdf'
    #
    cmd['rcParams'] = {'text.usetex': True,
                       'font.family': 'serif',
                       'lines.linewidth': 1,
                       'figure.dpi': 300,
                       'figure.figsize': (2.5, 1.5),
                       'font.size': 8,
                       'axes.spines.right': False,
                       'axes.spines.left': True,
                       'axes.spines.top': False,
                       'axes.spines.bottom': True,
                       }
    #
    XYs = read_XYs(cmd['ifn'])
    fig_XYs(XYs, cmd)
    return


def fig_XYs(XYs: dict, cmd: dict):
    import matplotlib.pyplot as plt
    plt.rcParams.update(cmd['rcParams'])
    #
    cmd['cmap'] = plt.get_cmap(cmd['cmap'])
    #
    fig = plt.figure()
    ax = fig.add_subplot()
    ax_XYs(ax, XYs, cmd)
    plt.savefig(cmd['ofn'], dpi=300, bbox_inches='tight')
    # plt.show()
    print('w: '+cmd['ofn'])
    return


def ax_XYs(ax, XYs: dict, cmd: dict):
    import numpy as np
    #
    cmap = cmd['cmap']
    cmd ['colors'] = []
    for itr in range(len(XYs['Y'])):
        cmd['colors'].append(cmap((itr+1)/(len(XYs['Y']))))
    #
    ax.stackplot(XYs['X'], XYs['Y'], baseline='zero',
                 colors=cmd['colors'], labels=cmd['labels'])
    ax.set_xscale('log')
    #
    # xrng, yrng = [1e-14, 1e-4], [0, 1]
    # xrng, yrng = [1e-14, 1e-6], [0, 1]
    #xrng, yrng = cmd['xrng'], cmd['yrng']
    xrng = [float(x) for x in cmd['xrng']]
    yrng = [float(x) for x in cmd['yrng']]
    #
    ax.set_xlim(xrng)
    Mtick, mtick = auto.tick.log(xrng)
    ax.set_xticks(Mtick)
    ax.set_xticks(mtick, minor=True)
    if 'xticklabels' in cmd:
        ax.set_xticklabels(cmd['xticklabels'])
    if 'yticklabels' in cmd:
        ax.set_xticklabels(cmd['yticklabels'])
    #
    ax.set_ylim(yrng)
    Mtick, mtick = auto.tick.linear(yrng)
    ax.set_yticks(Mtick)
    ax.set_yticks(mtick, minor=True)
    # # major tick
    # ax.tick_params(which='major', labelleft=True, labelbottom=True,
    #                top=False, bottom=False, left=False, right=False)
    # ax.grid(color='lightgray', which='major', linewidth=.5, zorder=1, axis='y')  # noqa
    # # minor tick
    # ax.tick_params(which='minor', labelleft=True, labelbottom=True,
    #                top=False, bottom=False, left=False, right=False)
    # ax.grid(color='lightgray', which='minor', linewidth=.1, zorder=1, axis='y')  # noqa
    #
    #
    # keys = list(XYs.keys())
    # tags = ['$10^{'+str(n)+'}$' for n in [int(np.log10(x)) for x in XYs[keys[0]]]]  # noqa
    # btm = np.array([0 for x in XYs[keys[1]]])
    # for key in keys[1:]:
    #     ax.bar(tags, XYs[key], bottom=btm, label=key, alpha=.3, zorder=2)
    #     btm = btm + np.array(XYs[key])
    if 'bbox_to_anchor' in cmd['legend']:
        ax.legend(loc=cmd['legend']['loc'],
                  fontsize=cmd['legend']['fontsize'],
                  ncol=int(cmd['legend']['ncol']),
                  bbox_to_anchor=tuple([float(x) for x in cmd['legend']['bbox_to_anchor']]))
    else:
        ax.legend(loc=cmd['legend']['loc'],
                  fontsize=cmd['legend']['fontsize'])
    # ax.legend(loc='upper left', fontsize='x-small')
    # ax.legend(loc='lower left', fontsize='x-small')
    # ax.legend(loc='lower right', fontsize='x-small', ncol=1, bbox_to_anchor=(1.4, 0))
    # label
    ax.set_xlabel(cmd['tlabel'])
    ax.set_ylabel(cmd['plabel'])
    return


def read_XYs(ifn: str):
    import numpy as np
    XYs = {'X': []}
    ifs = open(ifn, 'r')
    a_line = ifs.readline()
    line = a_line.rstrip('\n').split('\t')
    #
    Ys, keys = [], []
    for itr in range(1, len(line)):
        keys.append(itr)
        Ys.append([])
    #
    a_line = ifs.readline()
    while a_line:
        line = a_line.rstrip('\n').split('\t')
        XYs['X'].append(float(line[0]))
        for itr in keys:
            Ys[itr-1].append(float(line[itr]))
        a_line = ifs.readline()
    XYs['Y'] = np.vstack(Ys)
    ifs.close()
    print('r: '+ifn)
    return(XYs)


class auto():
    class tick():
        def linear(rng: list):
            import math
            pw = math.floor(math.log10(rng[1] - rng[0]))
            rng_sc = (rng[1] - rng[0]) / 10**pw
            if rng_sc < 3:
                pw = pw - 1
                # pw = math.copysign(abs(pw) - 1, pw)
                if rng_sc < 1.2:
                    i_m = .5 * (10**pw)
                    i_M = 2 * (10**pw)
                else:
                    i_m = 1 * (10**pw)
                    i_M = 5 * (10**pw)
            else:
                if rng_sc < 6:
                    i_m = .2 * (10**pw)
                    i_M = 1 * (10**pw)
                else:
                    i_m = .5 * (10**pw)
                    i_M = 2 * (10**pw)
            Mtick = [(i+int(rng[0]/i_M))*i_M for i in range(math.floor((rng[1] - rng[0]) / i_M)+1)]  # noqa
            mtick = [(i+int(rng[0]/i_m))*i_m for i in range(math.floor((rng[1] - rng[0]) / i_m)+1)]  # noqa
            return(Mtick, mtick)

        def log(rng: list):
            import math
            prng = [math.floor(math.log10(rng[0])), math.ceil(math.log10(rng[1]))]  # noqa
            Mtick, mtick = [], []
            for p in range(prng[0], prng[1], 1):
                Mtick += [10**p]
                mtick += [x*10**p for x in range(2, 10)]
            Mtick += [10**prng[1]]
            return(Mtick, mtick)


main()
