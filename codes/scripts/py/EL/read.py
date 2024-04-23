def print_open(ifn: str, open_type: str):
    print(open_type+': '+ifn)
    return(open(ifn, open_type))
# from MyPackage import print_open


class read:
    def energies(cmd):
        def read_EQs(ifn: str, e_unit=1.):
            EQs = {}
            s_unit = []
            lbls = []
            idxs = {}
            #
            itr = -1
            ifs = print_open(ifn, 'r')
            a_line = ifs.readline()
            while a_line:
                line = a_line.rstrip('\n').split('\t')
                if int(line[0]) in idxs:
                    print('error multiple keys: read_EQs with key:'+line[0])
                    exit()
                if line[2] not in s_unit:
                    s_unit.append(line[2])
                # TODO unit convert
                itr += 1
                EQs[itr] = float(line[1])*e_unit
                idxs[int(line[0])] = itr  # label to internal index key
                lbls.append(int(line[0]))  # label to internal index key
                #
                a_line = ifs.readline()
            if len(s_unit) == 1:
                print('EQ Input Energy Unit: '+s_unit[0])
                return(EQs, s_unit[0], idxs, lbls)
            else:
                print('WARNING: EQ Input Energy Unit is inconsistent: ', end='')
                print(s_unit)
                exit()

        def read_TSs(ifn: str, idxs: dict, lbls: list, e_unit=1.):
            TSs = {}
            s_unit = []
            pairs = []
            #
            ifs = print_open(ifn, 'r')
            a_line = ifs.readline()
            while a_line:
                line = a_line.rstrip('\n').split('\t')
                pair = tuple(sorted([int(line[1]), int(line[2])]))
                #
                for key in pair:
                    if key not in lbls:
                        print('terminated by input error')
                        print('TS has a key not in EQ: ', key)
                        exit()
                #
                # pair = tuple(sorted([line[1], line[2]]))
                pairs.append(pair)
                if line[4] not in s_unit:
                    s_unit.append(line[4])
                # TODO unit convert
                key = tuple([idxs[pair[0]], idxs[pair[1]]])
                if key not in TSs:
                    TSs[key] = [float(line[3])*e_unit]
                else:
                    TSs[key].append(float(line[3])*e_unit)
                #
                a_line = ifs.readline()
            if len(s_unit) == 1:
                print('TS Input Energy Unit: '+s_unit[0])
                return(TSs, s_unit[0], pairs)
            else:
                print('WARNING: TS Input Energy Unit is inconsistent: ', end='')
                print(s_unit)
                return(TSs, '')
        print('\n==== r: Energy ====')
        EQs, EQ_unit, idxs, lbls = read_EQs(cmd['EQs'])
        TSs, TS_unit, pairs = read_TSs(cmd['TSs'], idxs, lbls)
        #
        if EQ_unit != TS_unit or EQ_unit == '':
            print('terminated by input error')
            print('unit is not match btw EQ & TS: ', EQ_unit, TS_unit)
            exit()
        #
        print('===================\n')
        return({'EQs': EQs, 'TSs': TSs, 'idxs': idxs, 'lbls': lbls}, EQ_unit)

    def LInfLogP(cmd: dict, EQ_idxs: dict):
        mth_dic = {}
        idx_max = 0
        vals = []
        ifs = print_open(cmd['lmp'], 'r')
        a_line = ifs.readline()
        while a_line:
            line = a_line.rstrip('\n').split('\t')
            # pair = tuple([int(x) - 1 for x in line[:2]])
            pair = tuple([int(x) for x in line[:2]])
            for key in pair:
                if key not in EQ_idxs:
                    print('terminated by input error')
                    print('LInfLogP has a key not in EQ: ', key)
                    exit()
            pair = tuple([EQ_idxs[int(x)] for x in line[:2]])
            mth_dic[pair] = float(line[2])
            mth_dic[(pair[1], pair[0])] = float(line[2])
            vals.append(float(line[2]))
            idx_max = max(idx_max, max(pair))
            a_line = ifs.readline()
        return(mth_dic, sorted(vals))


def read_energy(cmd):
    print('\n==== r: Energy ====')
    EQs, EQ_unit = read_EQs(cmd['ifn_EQs'])
    TSs, TS_unit = read_TSs(cmd['ifn_TSs'])
    if EQ_unit != TS_unit or EQ_unit == '':
        print('terminated by input error')
        exit()
    print('===================\n')
    return({'EQs': EQs, 'TSs': TSs}, EQ_unit)


def read_EQs(ifn: str, e_unit=1.):
    EQs = {}
    s_unit = []
    #
    ifs = print_open(ifn, 'r')
    a_line = ifs.readline()
    while a_line:
        line = a_line.rstrip('\n').split('\t')
        key = int(line[0])
        if line[2] not in s_unit:
            s_unit.append(line[2])
        # TODO unit convert
        if key in EQs:
            print('error multiple keys: read_EQs with key:'+str(key))
        EQs[key] = float(line[1])*e_unit
        #
        a_line = ifs.readline()
    if len(s_unit) == 1:
        print('EQ Input Energy Unit: '+s_unit[0])
        return(EQs, s_unit[0])
    else:
        print('WARNING: EQ Input Energy Unit is inconsistent: ', end='')
        print(s_unit)
        return(EQs, '')


def read_TSs(ifn: str, e_unit=1.):
    TSs = {}
    s_unit = []
    #
    ifs = print_open(ifn, 'r')
    a_line = ifs.readline()
    while a_line:
        line = a_line.rstrip('\n').split('\t')
        key = tuple(sorted([int(line[1]), int(line[2])]))
        if line[4] not in s_unit:
            s_unit.append(line[4])
        # TODO unit convert
        if key not in TSs:
            TSs[key] = [float(line[3])*e_unit]
        else:
            TSs[key].append(float(line[3])*e_unit)
        #
        a_line = ifs.readline()
    if len(s_unit) == 1:
        print('TS Input Energy Unit: '+s_unit[0])
        return(TSs, s_unit[0])
    else:
        print('WARNING: TS Input Energy Unit is inconsistent: ', end='')
        print(s_unit)
        return(TSs, '')


def read_LInfLogP(ifn: str):
    import numpy as np
    mth_dic = {}
    idx_max = 0
    vals = []
    ifs = print_open(ifn, 'r')
    a_line = ifs.readline()
    while a_line:
        line = a_line.rstrip('\n').split('\t')
        key = tuple([int(x) - 1 for x in line[:2]])
        mth_dic[key] = float(line[2])
        vals.append(float(line[2]))
        idx_max = max(idx_max, max(key))
        a_line = ifs.readline()
    #
    mth_mtx = np.zeros((idx_max+1, idx_max+1))
    for key in mth_dic.keys():
        i1, i2 = min(key), max(key)
        mth_mtx[i1, i2] = mth_dic[key]
        mth_mtx[i2, i1] = mth_dic[key]
    return(mth_mtx, sorted(vals))
