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
        print('\nchecking configuration files...')
        print('===================')
        for conf in confs:
            print(cmd[conf])
            cmd[conf] = read_conf_file(cmd[conf])
        print('===================\n')
    return(cmd)


def check_keys(cmd: dict, keys_cmd: list):
    import os
    print('\nchecking keys...')
    print('===================')
    if len(cmd) == 0:
        return(True)
    #
    error = False
    for key in keys_cmd:
        print(key, end='\t')
        if key not in cmd:
            print('has no input')
            error = True
            continue
        if key[:2] == 'if':
            if not os.path.isfile(cmd[key]):
                print('file does NOT exist', end=':\t')
                error = True
            else:
                print('file exists', end=':\t')
            print(cmd[key])
        elif key == 'win_dir':
            print(cmd[key])
        elif key[:2] == 'of' or key[-3:] == 'dir':
            itr = cmd[key].rfind('/')
            if not os.path.isdir(cmd[key][:itr]):
                print('folder does NOT exist:\t' + cmd[key][:itr])
                error = True
            else:
                print('folder exists:\t' + cmd[key])
        else:
            print(cmd[key])
    print('===================')
    if not error:
        print('keys are all clear')
    else:
        print('input ERROR in check_keys.py')
    return(error)


def main():
    # from MyPackage import read_configs, check_keys
    import EL
    cmd = read_configs()
    check_keys(cmd, ['ofd'])
    print('\neng', end='')
    check_keys(cmd['ifn'], ['EQs', 'TSs', 'lmp'])

    cmd['fig']['rcParams'] = {'text.usetex': True,
                              'font.family': 'serif',
                              'font.serif': ['Computer Modern Roman'],
                              'lines.linewidth': 1,
                              'figure.dpi': 300,
                              'axes.spines.top': False,
                              }
    cmd['fig']['ofd'] = cmd['ofd']
    # read EQ energies
    Es, s_unit = EL.read.energies(cmd['ifn'])
    if s_unit == 'kJ/mol':
        Es['unit'] = 1e3
    #
    mtx_dst, vals = EL.read.LInfLogP(cmd['ifn'], Es['idxs'])
    cmd['lkg'] = {'n_max': len(Es['EQs']),
                  'threshold': 10*max(vals),
                  'undirect': 'min',
                  'linkage': 'complete',
                  }
    lkg = EL.DGs.from_dic(mtx_dst, cmd['lkg'])
    # labeling
    if 'vlabels' in cmd['fig']:
        lbls = {}
        for itr in range(len(Es['idxs'])):
            lbls[itr] = cmd['fig']['vlabels'][itr]
    else:
        lbls = Es['lbls']
    # ordering
    order = EL.DGs.order(cmd['fig'], lbls, Es['EQs'], odr_dir='asc')
    if 'tree_align' in cmd['fig']:
        G = EL.DGs.to_nwk(lkg['mtx'], order, lbls, Es, cmd['fig']['tree_align'])  # noqa
    else:
        G = EL.DGs.to_nwk(lkg['mtx'], order, lbls, Es)  # noqa
    #
    # print('\nEQ indices are labeled\nin\tlabel\torder\tposition')
    # for itr in range(len(lbls)):
    #     print(G.nodes[lbls[itr]]['idx'],
    #           lbls[itr],
    #           G.nodes[lbls[itr]]['odr'],
    #           G.nodes[lbls[itr]]['pos'],
    #           sep='\t')
    # print()
    #
    if 'tag' in cmd['fig']:
        EL.DGs.print_labels(G, lbls, cmd['ofd']+'/labels.'+cmd['fig']['tag']+'.tsv')
        EL.DGs.print_joins(G, cmd['ofd']+'/joins.'+cmd['fig']['tag']+'.tsv')
    else:
        EL.DGs.print_labels(G, lbls, cmd['ofd']+'/labels.tsv')
        EL.DGs.print_joins(G, cmd['ofd']+'/joins.tsv')
    #
    # EL.DGs.fig_gviz(G, cmd['fig'])
    EL.DGs.fig_MTH(G, cmd['fig'])
    #
    if 'nwk' in cmd:
        print('\nplotting network...')
        cmd['nwk']['ofd'] = cmd['ofd']
        cmd['nwk']['labels'] = cmd['fig']['vlabels']
        if 'tag' in cmd['fig']:
            cmd['nwk']['tag'] = cmd['fig']['tag']
        EL.DGs.fig_networks(G, cmd['nwk'], Es)
    #
    return


main()
