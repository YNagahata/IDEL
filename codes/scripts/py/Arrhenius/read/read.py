class read:
    def evals(ifn: str):
        print('r: '+ifn)
        ifs = open(ifn, 'r')
        a_line = ifs.readline()
        XY = {}
        keys = []
        for key in a_line.rstrip('\n').split('\t')[1:]:
            XY[key] = []
            keys.append(key)
        while a_line:
            if a_line[0] != '#':
                line = a_line.rstrip('\n').split('\t')
                XY[keys[0]].append(float(line[0])) # T / [degree of Celsius]
                XY[keys[1]].append(int(line[1])) # k : index of eigenvalue
                XY[keys[2]].append(float(line[2])) # lambda_k / [s] : k th eigenvalue
            a_line = ifs.readline()
        ifs.close()
        return(XY)

    def kf_exp(ifn: str, R = 8.31446261815324):
        ifs = open(ifn, 'r')
        a_line = ifs.readline()
        val = {}
        keys = []
        for key in a_line.rstrip('\n').split('\t')[1:]:
            val[key] = []
            keys.append(key)
        a_line = ifs.readline()
        while a_line:
            if a_line[0] != '#':
                line = a_line.rstrip('\n').split('\t')[1:]
                for itr in range(len(line)):
                    val[keys[itr]].append(float(line[itr]))
            a_line = ifs.readline()
        ifs.close()
        x = [1/(R*z) for z in val['kelvin']]
        y = [z*2.303/60 for z in val['k/2.303 min^-1']]
        return(x, y)

    def k_fwd(ifn: str):
        print('r: '+ifn)
        ifs = open(ifn, 'r')
        a_line = ifs.readline()
        XY = {}
        while a_line:
            if a_line[0] != '#':
                line = a_line.rstrip('\n').split('\t')
                XY[line[0]] = float(line[1])
            a_line = ifs.readline()
        ifs.close()
        return(XY)

    def rate_matrix(ifn: str):
        print('r: '+ifn)
        ifs = open(ifn, 'r')
        a_line = ifs.readline()
        rm = {}
        while a_line:
            if a_line[0] != '#':
                line = a_line.rstrip('\n').split('\t')
                if line[0] not in rm:
                    rm[line[0]] = {}
                rm[line[0]][(int(line[1]), int(line[2]))] = float(line[3])
            a_line = ifs.readline()
        ifs.close()
        return(rm)

    def ks(ifn: str):
        print('r: '+ifn)
        ifs = open(ifn, 'r')
        a_line = ifs.readline()
        XY = {}
        keys = []
        while a_line:
            if a_line[0] != '#':
                line = a_line.rstrip('\n').split('\t')
                if line[0] not in XY:
                    XY[line[0]] = {} # T / [degree of Celsius]
                XY[line[0]][(int(line[1]), int(line[2]))] = float(line[3]) # k_{i,j} / [s] : rate constant
            a_line = ifs.readline()
        ifs.close()
        return(XY)
