import networkx as nx


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
            MtickL = [(i+int(rng[0]/i_M))*i_M for i in range(math.floor((rng[1] - rng[0]) / i_M)+1)]  # noqa
            return(Mtick, mtick, MtickL)

        def log(rng: list):
            import math
            prng = [math.floor(math.log10(rng[0])), math.ceil(math.log10(rng[1]))]  # noqa
            #
            for u in [1,2,5,10,20]:
                if (prng[1]-prng[0])/u > 6:
                    continue
                u_M = u
                break
            #
            Mtick, mtick, MtickL = [], [], []
            for p in range(prng[0], prng[1]+1, 1):
                Mtick += [10**p]
                mtick += [x*10**p for x in range(2, 10)]
                if p%u_M == 0:
                    MtickL += ["$10^{"+str(p)+"}$"]
                else:
                    MtickL += [""]
            return(Mtick, mtick, MtickL)


class DGs:
    def from_dic(dst_dic, cmd):
        import numpy as np
        from scipy.cluster.hierarchy import linkage
        from scipy.spatial.distance import squareform
        key_miss = []
        for key in ['n_max', 'threshold', 'undirect', 'linkage']:
            if key not in cmd:
                key_miss.append(key)
        if len(key_miss) > 0:
            print('missing keys: '+', '.join([str(x) for x in key_miss]))
            quit()
        #
        n_max = cmd['n_max']
        #
        dst_out = {}
        mtx = np.full((n_max, n_max), cmd['threshold'])
        if cmd['undirect'] == 'min':
            for key in dst_dic.keys():
                i1, i2 = key[0], key[1]
                if max(i1, i2) > n_max - 1:
                    dst_out[(min(key), max(key))] = dst_dic[key]
                    continue
                elif i1 != i2:
                    val = min([dst_dic[i1, i2], dst_dic[i2, i1]])
                    mtx[i1, i2] = val
                    mtx[i2, i1] = val
        elif cmd['undirect'] == 'max':
            for key in dst_dic.keys():
                i1, i2 = key[0], key[1]
                if max(i1, i2) > n_max - 1:
                    dst_out[(min(key), max(key))] = dst_dic[key]
                    continue
                elif i1 != i2:
                    val = max([dst_dic[i1, i2], dst_dic[i2, i1]])
                    mtx[i1, i2] = val
                    mtx[i2, i1] = val
        for itr in range(n_max):
            mtx[itr, itr] = 0
        mtx_dst = squareform(mtx)
        lkg_out = {'mtx': linkage(mtx_dst, cmd['linkage']),
                   'DCs': dst_out}
        return(lkg_out)

    def get_weights(Es: dict):
        weights = {}
        DCs = {}
        E_MAX = 0
        for key in Es['TSs'].keys():
            for E_TS in Es['TSs'][key]:
                if E_MAX < E_TS:
                    E_MAX = E_TS
        for key in Es['TSs'].keys():
            if len(key) != 2:
                print('ERROR: TST_SI: len(key) != 2 : ', key)
                exit()
            idxs = []
            for idx in key:
                if idx in Es['EQs']:
                    idxs.append(idx)
            if len(idxs) == 2:
                weights[(key[1], key[0])] = 0
                weights[(key[0], key[1])] = 0
                for E_TS in Es['TSs'][key]:
                    dG = (E_TS - Es['EQs'][key[0]])  # * E_unit
                    if weights[(key[1], key[0])] < dG:
                        weights[(key[1], key[0])] = dG
                    dG = (E_TS - Es['EQs'][key[1]])  # * E_unit
                    if weights[(key[0], key[1])] < dG:
                        weights[(key[0], key[1])] = dG
            elif len(idxs) == 1:
                dG = (Es['TSs'][key][0] - Es['EQs'][idxs[0]])  # * E_unit
                DCs[idxs[0]] = dG
        #
        keys, itrs, itr = {}, {}, -1
        for key in Es['EQs'].keys():
            itr += 1
            keys[itr] = key
            itrs[key] = itr
        mtx = {}
        for key in weights.keys():
            mtx[(itrs[key[0]], itrs[key[1]])] = weights[key]
        #
        rng = {'min': min(weights.values()), 'max': max(weights.values())}
        edge_out = {'rng': rng,
                    'weights': weights,
                    'mtx': mtx,
                    'keys': keys,
                    'itrs': itrs,
                    'wDCs': DCs,
                    }
        return(edge_out)

    def fig_lkg(lkg_mtx: list, cmd: dict):
        from scipy.cluster.hierarchy import dendrogram
        from scipy.cluster.hierarchy import set_link_color_palette
        import numpy as np
        import matplotlib.pyplot as plt
        #
        key_miss = []
        for key in ['rcParams', 'vrange', 'vlabels', 'orientation',
                    'color_threshold', 'vlabel', 'ofn_fig']:
            if key not in cmd:
                key_miss.append(key)
        if len(key_miss) > 0:
            print('missing keys in lkg.fig_ddgm: '
                  + ', '.join([str(x) for x in key_miss]))
            quit()
        #
        # import math
        plt.rcParams.update(cmd['rcParams'])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if 'title' in cmd:
            plt.title(cmd['title'])
        if 'text' in cmd:
            ax.text(float(cmd['text'][0]), float(cmd['text'][1]),
                    cmd['text'][2])
        vrng = [float(x) for x in cmd['vrange']]
        #
        if 'scale' in cmd:
            if cmd['scale'] == 'log':
                Mtick, mtick, MtickL = auto.tick.log(vrng)
            else:
                Mtick, mtick, MtickL = auto.tick.linear(vrng)
        else:
            Mtick, mtick, MtickL = auto.tick.linear(vrng)
        ax.tick_params(axis='both', which='minor', labelsize=0)
        #
        set_link_color_palette(['k'])
        ret = dendrogram(lkg_mtx, ax=ax,
                         labels=cmd['vlabels'],
                         orientation=cmd['orientation'],
                         color_threshold=float(cmd['color_threshold']),
                         above_threshold_color='lightgray')
        x = np.mean(ret['icoord'][-1][1:3])
        y1 = ret['dcoord'][-1][1]
        if cmd['orientation'] == 'top' or cmd['orientation'] == 'bottom':
            if cmd['scale'] == 'log':
                ax.set_yscale('log')
            ax.set_ylabel(cmd['vlabel'])
            plt.ylim(vrng[0], vrng[1])
            ax.set_yticks(Mtick)
            ax.set_yticklabels(MtickL)
            ax.set_yticks(mtick, minor=True)
            # get lumped position
            for lst in ret['dcoord']:
                ax.axhline(lst[1], c='gray', lw=.1)
            #
            x = np.mean(ret['icoord'][-1][1:3])
            ax.plot([x, x], [y1, vrng[1]], color=ret['color_list'][-1])
        elif cmd['orientation'] == 'right' or cmd['orientation'] == 'left':
            if cmd['scale'] == 'log':
                ax.set_xscale('log')
            ax.set_xlabel(cmd['vlabel'])
            plt.xlim(vrng[0], vrng[1])
            ax.set_xticks(Mtick)
            ax.set_xticklabels(MtickL)
            ax.set_xticks(mtick, minor=True)
            # get lumped position
            for lst in ret['dcoord']:
                ax.axvline(lst[1], c='gray', lw=.1)
            #
            x = np.mean(ret['icoord'][-1][1:3])
            ax.plot([y1, vrng[1]], [x, x], color=ret['color_list'][-1])
        #
        plt.savefig(cmd['ofn_fig']+'.lkg.pdf', dpi=300, bbox_inches='tight')
        print('w: '+cmd['ofn_fig']+'.lkg.pdf')
        return

    def fig_gviz(G_in: nx.DiGraph, cmd: dict):
        Graph = G_in.copy()
        # Change the shape of leaves
        for key in Graph.nodes:
            if type(key) is int:
                Graph.nodes[key]['shape'] = 'box'
        # Draw
        gviz = nx.nx_agraph.to_agraph(Graph)
        gviz.draw(cmd['ofn_fig']+'.pdf', prog='dot')
        print('gviz: '+cmd['ofn_fig']+'.gviz.pdf')
        return

    def order(cmd: dict, lbls: dict, vals={}, odr_dir='asc'):
        odr_out = {}
        if 'order' in cmd:
            for itr in range(len(cmd['order'])):
                odr_out[cmd['order'][itr]] = itr
        elif len(vals) != 0:
            if odr_dir == 'asc':
                itr = -1
                for pair in sorted(vals.items(), key=lambda x: x[1]):
                    itr += 1
                    odr_out[lbls[pair[0]]] = (itr)
            elif odr_dir == 'dsc':
                itr = -1
                for pair in sorted(vals.items(), key=lambda x: x[1], reverse=True):  # noqa 501
                    itr += 1
                    odr_out[lbls[pair[0]]] = (itr)
            else:
                print('EL.DGs.order: odr_dir: type error')
        else:
            print('EL.DGs.order: no order input in cmd or vals')
        return(odr_out)

    def to_nwk(mtx_lkg: list, order: list, lbls: dict, Es: dict, pred_pos='mid'):  # noqa E501
        def set_odr(G_in: nx.DiGraph, order: list, lbls: dict):
            # set 'odr' of leaves
            for key in lbls.keys():
                G_in.nodes[lbls[key]]['odr'] = order[lbls[key]]
            # sort joins by length of the key: from smaller to larger
            joins = sorted(list(G_in.nodes['root']['joins']), key=len)
            # set 'odr' of joins
            for join in joins:
                odrs = [G_in.nodes[key]['odr'] for key in G_in.succ[join].keys()]  # noqa 501
                G_in.nodes[join]['odr'] = min(odrs)
            return(G_in)

        def set_pos(G_in: nx.DiGraph, pred_pos='mid'):
            poss = {}  # position according to odr key
            # get value from the root node
            lst = list(G_in.succ['root'].keys())
            odr_pred = G_in.nodes[lst[0]]['odr']
            G_in.nodes[lst[0]]['pos'] = 1
            poss[odr_pred] = 1
            #
            if pred_pos == 'mid':
                while len(lst) > 0:
                    key = lst.pop(0)
                    keys = list(G_in.succ[key].keys())
                    if len(keys) < 2:
                        G_in.nodes[key]['pos'] = poss[G_in.nodes[key]['odr']]
                        continue
                    # key length
                    key_len = []
                    for k in keys:
                        if type(k) is int:
                            key_len.append(1)
                        else:
                            key_len.append(len(k))
                    # 'pos'
                    odr_pred = G_in.nodes[key]['odr']
                    if G_in.nodes[keys[0]]['odr'] == odr_pred:
                        margin = poss[odr_pred] + key_len[0]
                        G_in.nodes[key]['pos'] = margin - .5
                        poss[G_in.nodes[keys[1]]['odr']] = margin
                    else:
                        margin = poss[odr_pred] + key_len[1]
                        G_in.nodes[key]['pos'] = margin - .5
                        poss[G_in.nodes[keys[0]]['odr']] = margin
                    lst += [keys[0], keys[1]]
            elif pred_pos == 'center':
                while len(lst) > 0:
                    key = lst.pop(0)
                    keys = list(G_in.succ[key].keys())
                    if len(keys) < 2:
                        G_in.nodes[key]['pos'] = poss[G_in.nodes[key]['odr']]
                        continue
                    # key length
                    key_len = []
                    for k in keys:
                        if type(k) is int:
                            key_len.append(1)
                        else:
                            key_len.append(len(k))
                    # 'pos'
                    odr_pred = G_in.nodes[key]['odr']
                    G_in.nodes[key]['pos'] =\
                        poss[odr_pred] + (key_len[0] + key_len[1] - 1)/2
                    if G_in.nodes[keys[0]]['odr'] == odr_pred:
                        poss[G_in.nodes[keys[1]]['odr']] =\
                            poss[odr_pred] + key_len[0]
                    else:
                        poss[G_in.nodes[keys[0]]['odr']] =\
                            poss[odr_pred] + key_len[1]
                    lst += [keys[0], keys[1]]
            elif pred_pos == 'align':
                while len(lst) > 0:
                    key = lst.pop(0)
                    keys = list(G_in.succ[key].keys())
                    if len(keys) < 2:
                        continue
                    # key length
                    key_len = []
                    for k in keys:
                        if type(k) is str:
                            key_len.append(1)
                        else:
                            key_len.append(len(k))
                    # 'pos'
                    odr_pred = G_in.nodes[key]['odr']
                    mar = poss[odr_pred]
                    if G_in.nodes[keys[0]]['odr'] == odr_pred:
                        G_in.nodes[keys[0]]['pos'] = mar
                        G_in.nodes[keys[1]]['pos'] = mar + key_len[0]
                        poss[G_in.nodes[keys[1]]['odr']] = mar + key_len[0]
                    else:
                        G_in.nodes[keys[1]]['pos'] = mar
                        G_in.nodes[keys[0]]['pos'] = mar + key_len[1]
                        poss[G_in.nodes[keys[0]]['odr']] = mar + key_len[1]
                    lst += [keys[0], keys[1]]
            elif pred_pos != 'branch':
                print("ERROR lkg.set_pos: ", end='')  # noqa 501
                print("pred_pos is not 'mid', 'center', 'align' nor 'branch'")  # noqa 501
                exit
            return(G_in)

        Graph = nx.DiGraph()
        Graph.nodes(data=True)
        clstr_key, clstr = [], []
        joins = []
        # for itr in range(len(mtx_lkg) + 1):
        for itr in range(len(lbls)):
            clstr_key.append(lbls[itr])
            clstr.append(tuple([lbls[itr]]))
        # leaves = clstr_key
        for ele in mtx_lkg:
            clstrs = [clstr[int(ele[0])], clstr[int(ele[1])]]
            keys = [clstr_key[int(ele[0])], clstr_key[int(ele[1])]]
            key_new = tuple(sorted(list(set(clstrs[0]) | set(clstrs[1]))))
            clstr_key.append(key_new)
            clstr.append((key_new))
            # construct new node with edges
            Graph.add_node(key_new, weight=ele[2])
            Graph.add_edge(key_new, keys[0])
            Graph.add_edge(key_new, keys[1])
            # keep the info of the construction
            joins.append(key_new)
        # construct the root node
        Graph.add_edge('root', key_new)
        Graph.nodes['root']['joins'] = joins
        leaves = [x for x in Graph.nodes() if Graph.out_degree(x) == 0 and Graph.in_degree(x) == 1]  # noqa E501
        Graph.nodes['root']['leaves'] = leaves
        Graph.nodes['root']['lbls'] = lbls
        for itr in range(len(lbls)):
            Graph.nodes[lbls[itr]]['En'] = Es['EQs'][itr]
            Graph.nodes[lbls[itr]]['idx'] = itr # Es['idxs'][itr]
        #
        Graph = set_odr(Graph, order, lbls)
        Graph = set_pos(Graph, pred_pos=pred_pos)
        #
        return(Graph)

    def fig_networks(G_in: nx.DiGraph, cmd: dict, Es: dict):
        import matplotlib.pyplot as plt
        import numpy as np
        import matplotlib.patches as patches
        #
        #  DEFAULT  #
        s = 0.1  # size of the nodes
        rcParams = {'font.family': 'suns-serif',
                    'xtick.top': False,
                    'xtick.bottom': False,
                    'ytick.left': False,
                    'ytick.right': False,
                    'xtick.labeltop': False,
                    'xtick.labelbottom': False,
                    'ytick.labelleft': False,
                    'ytick.labelright': False,
                    }
        plt.rcParams.update(rcParams)
        #  INPUTS  #
        #
        belonging = {}
        i_cls = -1
        clusters = []
        Rs = []
        for clstr in cmd['cluster']:
            i_cls += 1
            cluster = {}
            for key in [x for x in clstr.split(' ')]:
                belonging[key] = {'cls': i_cls,
                                  'pos': G_in.nodes[key]['pos']}
                cluster[G_in.nodes[key]['pos']] = {'idx': G_in.nodes[key]['idx'],
                                                   'lbl': key}
            clusters.append({'ele': cluster,
                             'inf': {'R': 1.5*s/np.sin(np.pi/len(cluster))}})
            Rs.append(s/np.sin(np.pi/len(cluster)))
        Rs = sorted(Rs)
        R_cls = 2*(Rs[-1] + Rs[-2])
        #
        # input graph allocations
        itr0 = -1
        for itr in range(len(clusters)):
            itr0 += 1
            clusters[itr]['inf']['alc'] = [R_cls*np.cos(itr0*2*np.pi/len(clusters)),
                                           R_cls*np.sin(itr0*2*np.pi/len(clusters))]
            n = len(list(clusters[itr]['ele'].keys()))
            R = clusters[itr]['inf']['R']
            itr1 = -1
            for pos in sorted(list(clusters[itr]['ele'].keys())):
                itr1 += 1
                clusters[itr]['ele'][pos]['alc'] = [R*np.cos(itr1*2*np.pi/n),
                                                    R*np.sin(itr1*2*np.pi/n)]
        #
        #  FIGURE DRAWING  #
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        ax.set_xlim(-2*R_cls, 2*R_cls)
        ax.set_ylim(-2*R_cls, 2*R_cls)
        # regular polygon
        X, Y = [], []
        for cluster in clusters:
            lst = sorted(list(cluster['ele'].keys()))
            for key in lst:
                X.append(cluster['ele'][key]['alc'][0]+cluster['inf']['alc'][0])
                Y.append(cluster['ele'][key]['alc'][1]+cluster['inf']['alc'][1])
            for key in lst[:1]:
                X.append(cluster['ele'][key]['alc'][0]+cluster['inf']['alc'][0])
                Y.append(cluster['ele'][key]['alc'][1]+cluster['inf']['alc'][1])
            X.append(np.nan)
            Y.append(np.nan)
            ax.plot(X, Y)
        # lines from center
        X, Y = [], []
        for cluster in clusters:
            for key in cluster['ele'].keys():
                X.append(cluster['ele'][key]['alc'][0]+cluster['inf']['alc'][0])
                Y.append(cluster['ele'][key]['alc'][1]+cluster['inf']['alc'][1])
                X.append(1.5*cluster['ele'][key]['alc'][0]+cluster['inf']['alc'][0])
                Y.append(1.5*cluster['ele'][key]['alc'][1]+cluster['inf']['alc'][1])
                X.append(np.nan)
                Y.append(np.nan)
            ax.plot(X, Y)
        # node
        X, Y = [], []
        for cluster in clusters:
            for key in cluster['ele'].keys():
                x = cluster['ele'][key]['alc'][0]+cluster['inf']['alc'][0]
                y = cluster['ele'][key]['alc'][1]+cluster['inf']['alc'][1]
                ax.add_patch(patches.Circle(xy=(x, y), radius=s,
                                            fc='w', ec='k', zorder=10))
                ax.annotate(cluster['ele'][key]['lbl'], xy=(x, y),
                            fontsize=10, va="center", ha="center", zorder=11)
                # tag = input id
                x = 1.5*cluster['ele'][key]['alc'][0]+cluster['inf']['alc'][0]
                y = 1.5*cluster['ele'][key]['alc'][1]+cluster['inf']['alc'][1]
                ax.annotate(cluster['ele'][key]['idx'], xy=(x, y),
                            fontsize=10, va="center", ha="center", zorder=11)
        # edge
        # TSs = read_TSs(cmd['ifn_TSs'])
        X, Y = [], []
        # for key in TSs.keys():
        for key in Es['TSs'].keys():
            key0 = G_in.nodes['root']['lbls'][key[0]]
            key1 = G_in.nodes['root']['lbls'][key[1]]
            cluster = clusters[belonging[key0]['cls']]
            pos = belonging[key0]['pos']
            X.append(cluster['ele'][pos]['alc'][0]+cluster['inf']['alc'][0])
            Y.append(cluster['ele'][pos]['alc'][1]+cluster['inf']['alc'][1])
            cluster = clusters[belonging[key1]['cls']]
            pos = belonging[key1]['pos']
            X.append(cluster['ele'][pos]['alc'][0]+cluster['inf']['alc'][0])
            Y.append(cluster['ele'][pos]['alc'][1]+cluster['inf']['alc'][1])
            X.append(np.nan)
            Y.append(np.nan)
        ax.plot(X, Y, color='gray')
        #
        if 'ofn_fig' not in cmd:
            if 'tag' in cmd:
                cmd['ofn_fig'] = cmd['ofd']+'/EL.DGs.fig.MTH.nwk.'+cmd['tag']+'.pdf'
            else:
                cmd['ofn_fig'] = cmd['ofd']+'/EL.DGs.fig.MTH.nwk.pdf'
        plt.savefig(cmd['ofn_fig'], dpi=300, bbox_inches='tight')
        print('w: '+cmd['ofn_fig'])
        return

    def fig_MTH(G_in: nx.DiGraph, cmd: dict):
        import matplotlib.pyplot as plt
        from matplotlib import gridspec

        def set_rcParams(cmd: dict, n_leaves: int):
            if 'font.size' in cmd:
                cmd['rcParams']['font.size'] = float(cmd['font.size'])
            # fig.size
            if type(cmd['figure.figsize']) == list and len(cmd['figure.figsize']) == 2:  # noqa E501
                cmd['rcParams']['figure.figsize'] = [float(cmd['figure.figsize'][0]),  # noqa E501
                                                     float(cmd['figure.figsize'][1])]  # noqa E501
            elif cmd['orientation'] == 'top' or cmd['orientation'] == 'bottom':
                cmd['rcParams']['figure.figsize'] = [n_leaves/9,
                                                     float(cmd['figure.figsize'])]  # noqa E501
            elif cmd['orientation'] == 'right' or cmd['orientation'] == 'left':
                cmd['rcParams']['figure.figsize'] = [float(cmd['figure.figsize']),  # noqa E501
                                                     n_leaves/9]
            # print
            print('\n==== rcParams =====')
            for key in cmd['rcParams']:
                print(key, cmd['rcParams'][key])
            print('===================\n')
            return(cmd)

        def XY_MTH(G_in: nx.DiGraph, cmd: dict):
            import numpy as np
            # from matplotlib.patches import Polygon
            # import math
            leaves = G_in.nodes['root']['leaves']
            joins = G_in.nodes['root']['joins']
            vrng = [float(x) for x in cmd['vrange']]
            # plot from leaves
            d_lvs = {'X': [], 'Y': []}
            for key in leaves:
                pred = list(G_in.pred[key].keys())[0]
                X = [G_in.nodes[key]['pos'], G_in.nodes[key]['pos']]
                Y = [vrng[0], G_in.nodes[pred]['weight']]
                if G_in.nodes[pred]['weight'] == float(cmd['color_threshold']):
                    continue
                if G_in.nodes[key]['pos'] != G_in.nodes[pred]['pos']:
                    X.append(G_in.nodes[pred]['pos'])
                    Y.append(G_in.nodes[pred]['weight'])
                d_lvs['X'] += X + [np.nan]
                d_lvs['Y'] += Y + [np.nan]
            #
            # plot from joins
            d_jns = []
            for join in joins:
                # dendrogram
                pred = list(G_in.pred[join].keys())[0]
                if pred != 'root':
                    if G_in.nodes[pred]['weight'] == float(cmd['color_threshold']):  # noqa E501
                        continue
                    X_fwd = [G_in.nodes[join]['pos'], G_in.nodes[join]['pos']]
                    Y_fwd = [G_in.nodes[join]['weight'], G_in.nodes[pred]['weight']]  # noqa E501
                    if G_in.nodes[join]['pos'] != G_in.nodes[pred]['pos']:
                        X_fwd.append(G_in.nodes[pred]['pos'])
                        Y_fwd.append(G_in.nodes[pred]['weight'])
                else:
                    if G_in.nodes[join]['weight'] < vrng[1]:  # noqa E501
                        X_fwd = [G_in.nodes[join]['pos'], G_in.nodes[join]['pos']]  # noqa E501
                        Y_fwd = [G_in.nodes[join]['weight'], vrng[1]]
                    else:
                        continue
                d_jns.append({'X': X_fwd, 'Y': Y_fwd, 'len': len(join)})
            #
            XYs = {'leaves': d_lvs,
                   'joins': d_jns,
                   }
            return(XYs)

        def ax_MTH(ax, G_in: nx.DiGraph, cmd: dict, show_basin):
            import matplotlib.pyplot as plt
            import numpy as np
            ax.spines['right'].set_visible(False)
            #
            if 'title' in cmd:
                ax.set_title(cmd['title'])
            vrng = [float(x) for x in cmd['vrange']]
            #
            # tick
            if 'scale' in cmd:
                if cmd['scale'] == 'log':
                    Mtick, mtick, MtickL = auto.tick.log(vrng)
                else:
                    Mtick, mtick, MtickL = auto.tick.linear(vrng)
            else:
                Mtick, mtick, MtickL = auto.tick.linear(vrng)
            #
            # label
            ktick, ktlbl = [], []
            if 'labels' in cmd: # TODO
                for key in G_in.nodes['root']['leaves']:
                    ktlbl.append(key)
                    ktick.append(G_in.nodes[key]['pos'])
            else:
                for key in G_in.nodes['root']['leaves']:
                    ktlbl.append(str(key))
                    ktick.append(G_in.nodes[key]['pos'])
            #
            XYs = XY_MTH(G_in, cmd)
            # text
            if 'text' in cmd:
                text = {}
                #
                if type(cmd['text']) != list:
                    text['text'] = cmd['text']
                else:
                    text['text'] = cmd['text'][-1]
                    if len(cmd['text']) == 3:
                        text['x'] = float(cmd['text'][0])
                        text['y'] = float(cmd['text'][1])
                #
                if len(text.keys()) != 3:
                    text['x'] = max(XYs['leaves']['X'])-.5
                    text['y'] = vrng[1]*.9
                    if 'scale' in cmd:
                        if cmd['scale'] == 'log':
                            text['y'] = pow(10, np.log10(vrng[1]/vrng[0])*.9)*vrng[0]  # noqa
            #
            lw_dgm = 2
            if 'lw_dgm' in cmd:
                lw_dgm = float(cmd['lw_dgm'])
            lw_tir = 1
            if 'lw_tir' in cmd:
                lw_tir = float(cmd['lw_tir'])
            #
            if cmd['orientation'] == 'top' or cmd['orientation'] == 'bottom':
                if 'scale' in cmd:
                    if cmd['scale'] == 'log':
                        ax.set_yscale('log')
                ax.set_ylabel(cmd['vlabel'])
                ax.set_yticks(Mtick)
                ax.set_yticklabels(MtickL)
                ax.set_yticks(mtick, minor=True)
                ax.set_yticklabels(["" for x in mtick], minor=True)
                ax.set_xticks(ktick)
                if show_basin:
                    ax.set_xticklabels([])  # noqa 501
                else:
                    ax.set_xticklabels(ktlbl, fontsize='small', rotation=90, va='top', ha='center')  # noqa 501
                plt.ylim(vrng[0], vrng[1])
                plt.xlim(.5, max(XYs['leaves']['X'])+.5)
                ax.text(text['x'], text['y'], text['text'])
                #
                # dendrogram
                ax.plot(XYs['leaves']['X'], XYs['leaves']['Y'],
                        lw=lw_dgm, color='k', solid_capstyle='butt', zorder=3)
                for XY in XYs['joins']:
                    ax.plot(XY['X'], XY['Y'], color='k', solid_capstyle='butt',
                            zorder=3, lw=lw_dgm)
                    ax.plot(XY['X'], XY['Y'], color='lightgray', solid_capstyle='butt',
                            zorder=2, lw=(XY['len']-1)*lw_tir+lw_dgm)
            elif cmd['orientation'] == 'right' or cmd['orientation'] == 'left':
                if cmd['scale'] == 'log':
                    ax.set_xscale('log')
                ax.set_xlabel(cmd['vlabel'])
                ax.set_xticks(Mtick)
                ax.set_xticklabels(MtickL)
                ax.set_xticks(mtick, minor=True)
                ax.set_xticklabels(["" for x in mtick], minor=True)
                ax.set_yticks(ktick)
                ax.set_yticklabels(ktlbl, fontsize='small', va='center', ha='right')  # noqa 501
                plt.xlim(vrng[0], vrng[1])
                plt.ylim(.5, max(XYs['leaves']['X'])+.5)
                ax.text(text['y'], text['x'], text['text'])
                #
                # dendrogram
                ax.plot(XYs['leaves']['Y'], XYs['leaves']['X'],
                        lw=lw_dgm, color='k', solid_capstyle='butt', zorder=3)
                for XY in XYs['joins']:
                    ax.plot(XY['Y'], XY['X'], color='k', solid_capstyle='butt',
                            zorder=3, lw=lw_dgm)
                    ax.plot(XY['Y'], XY['X'], color='lightgray', solid_capstyle='butt',
                            zorder=2, lw=(XY['len']-1)*lw_tir+lw_dgm)
            return

        def ax_basins(ax, G_in: nx.DiGraph, cmd: dict):
            import numpy as np
            #
            ax.spines['right'].set_visible(True)
            # EQ bars
            d_EQs, EQs = [], []
            for key in G_in.nodes['root']['leaves']:
                EQs.append(G_in.nodes[key]['En'])
            for key in G_in.nodes['root']['leaves']:
                if G_in.nodes[key]['En'] != min(EQs):
                    d_EQ = {'col': 'gray'}
                    d_EQ['X'] = [G_in.nodes[key]['pos']-.3,
                                 G_in.nodes[key]['pos']+.3,
                                 G_in.nodes[key]['pos']+.3,
                                 G_in.nodes[key]['pos']-.3]
                    d_EQ['Y'] = [G_in.nodes[key]['En'], G_in.nodes[key]['En'], min(EQs), min(EQs)]  # noqa
                    d_EQs.append(d_EQ)
            # plot
            if 'erange' in cmd:
                if len(cmd['erange']) == 2:
                    erng = [float(x) for x in cmd['erange']]
                elif len(cmd['erange']) == 1:
                    erng = [float(cmd['erange'][0]), max(EQs)]
                else:
                    erng = [min(EQs), max(EQs)]
            else:
                erng = [min(EQs), max(EQs)]
            Mtick, mtick, MtickL = auto.tick.linear(erng)
            # EQ lines
            l_EQs = {'X': [], 'Y': [], 'col': 'black'}
            for key in G_in.nodes['root']['leaves']:
                l_EQs['X'] += [G_in.nodes[key]['pos'], G_in.nodes[key]['pos'], np.nan]  # noqa
                l_EQs['Y'] += [erng[0], erng[1], np.nan]  # noqa
            # label
            ktick, ktlbl = [], []
            if 'labels' in cmd:
                for key in G_in.nodes['root']['leaves']:
                    ktick.append(G_in.nodes[key]['pos'])
                    if G_in.nodes[key]['En'] != min(EQs):
                        ktlbl.append(cmd['labels'][key])
                    else:
                        ktlbl.append('GM '+cmd['labels'][key])
            else:
                for key in G_in.nodes['root']['leaves']:
                    ktick.append(G_in.nodes[key]['pos'])
                    if G_in.nodes[key]['En'] != min(EQs):
                        ktlbl.append(str(key))
                    else:
                        ktlbl.append('GM '+str(key))
            #
            lw_EQs = 5
            if 'lw_EQs' in cmd:
                lw_EQs = float(cmd['lw_EQs'])
            #
            if cmd['orientation'] == 'top' or cmd['orientation'] == 'bottom':
                plt.ylim(erng[0], erng[1])
                plt.xlim(.5, len(G_in.nodes['root']['leaves'])+.5)
                ax.set_ylabel(cmd['elabel'])
                ax.set_yticks(Mtick)
                ax.set_yticklabels(MtickL)
                ax.set_yticks(mtick, minor=True)
                ax.set_xticks(ktick)
                ax.set_xticklabels(ktlbl, fontsize='small', rotation=90, va='top', ha='center')  # noqa 501
                for d_EQ in d_EQs:
                    ax.fill(d_EQ['X'], d_EQ['Y'], color=d_EQ['col'], zorder=1)
                # ax.text(d_EQ_X_gm, min(EQs), 'Global Minimum', color='gray')
                ax.plot(l_EQs['X'], l_EQs['Y'], color=l_EQs['col'], zorder=2,
                        lw=lw_EQs, solid_capstyle='butt')
            elif cmd['orientation'] == 'right' or cmd['orientation'] == 'left':
                plt.xlim(erng[0], erng[1])
                plt.ylim(.5, len(G_in.nodes['root']['leaves'])+.5)
                ax.set_xlabel(cmd['elabel'])
                ax.set_xticks(Mtick)
                ax.set_xticklabels(MtickL)
                ax.set_xticks(mtick, minor=True)
                ax.set_yticks(ktick)
                ax.set_yticklabels(ktlbl, fontsize='small', va='center', ha='right')  # noqa 501
                #
                for d_EQ in d_EQs:
                    ax.fill(d_EQ['Y'], d_EQ['X'], color=d_EQ['col'], zorder=1)
                #
                ax.plot(l_EQs['Y'], l_EQs['X'], color=l_EQs['col'], zorder=2,
                        lw=lw_EQs, solid_capstyle='butt')
            return
        #
        show_basin = False
        if 'show_basin' in cmd:
            if cmd['show_basin'] == 'True':
                show_basin = True
            elif cmd['show_basin'] != 'False':
                print('fig.conf: show_basin is inappropriate')
        # rcParams
        cmd = set_rcParams(cmd, len(G_in.nodes['root']['leaves']))
        plt.rcParams.update(cmd['rcParams'])
        #
        fig = plt.figure()
        if cmd['orientation'] == 'top' or cmd['orientation'] == 'bottom':
            spec = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[11, 3], hspace=.0)  # noqa
            ax = fig.add_subplot(spec[0])
            ax_MTH(ax, G_in, cmd, show_basin)
            if show_basin:
                ax_basins(fig.add_subplot(spec[1]), G_in, cmd)
        elif cmd['orientation'] == 'right' or cmd['orientation'] == 'left':
            spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[3, 11], wspace=.0)  # noqa
            ax = fig.add_subplot(spec[1])
            ax_MTH(ax, G_in, cmd, show_basin)
            if show_basin:
                ax_basins(fig.add_subplot(spec[0]), G_in, cmd)
        #
        if 'ofn_fig' not in cmd:
            tr_al = ''
            if 'tree_align' in cmd:
                tr_al = cmd['tree_align']
            if 'tag' in cmd:
                tr_al += '.'+cmd['tag']
            cmd['ofn_fig'] = cmd['ofd']+'/EL.DGs.fig.MTH.'+tr_al+'.pdf'
        plt.savefig(cmd['ofn_fig'], bbox_inches='tight')
        print('\nw: '+cmd['ofn_fig'])
        # plt.show()
        return

    def print_labels(G_in: nx.DiGraph, lbls: dict, ofn: str):
        ofs = open(ofn, 'w')
        ofs.write('label\torder\tposition\n')
        for itr in range(len(lbls)):
            out = [G_in.nodes[lbls[itr]]['idx'],
                   lbls[itr],
                   G_in.nodes[lbls[itr]]['odr'],
                   G_in.nodes[lbls[itr]]['pos']]
            ofs.write('\t'.join([str(x) for x in out])+'\n')
        print('w: '+ofn)
        ofs.close()
        return

    def print_joins(G_in: nx.DiGraph, ofn: str):
        joins = []
        for key in G_in.nodes['root']['joins']:
            join = [G_in.nodes[key]['weight']]
            for node in list(G_in.succ[key].keys()):
                if type(node) is tuple:
                    join.append(" ".join([str(x) for x in node]))
                else:
                    join.append(node)
            joins.append(join)
        ofs = open(ofn, "w")
        ofs.write('${\Delta t}_{I,J}$\t$I$\t$J$\n')
        for val in sorted(joins, key=lambda x: x[0]):
            ofs.write('\t'.join([str(x) for x in val])+'\n')
        print('w: '+ofn)
        ofs.close()
        return
