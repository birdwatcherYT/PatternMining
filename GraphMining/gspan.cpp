#include "gspan.h"

using namespace std;

namespace GSPAN {

    gSpan::gSpan(const string& filename, unsigned int minsup, unsigned int maxpat, bool directed) {
        this->minsup = minsup;
        this->maxpat = maxpat;
        this->directed = directed;
        ifstream fp(filename);
        read(fp);
        // for single vertex pattern
        for (unsigned int id = 0; id < TRANS.size(); ++id) {
            for (unsigned int nid = 0 ; nid < TRANS[id].size() ; ++nid) {
                if (singleVertex[id][TRANS[id][nid].label] == 0) {
                    // number of graphs it appears in
                    singleVertexLabel[TRANS[id][nid].label] += 1;// support
                }
                singleVertex[id][TRANS[id][nid].label] += 1;//count
            }
        }
        for (map<unsigned int, unsigned int>::iterator it = singleVertexLabel.begin () ; it != singleVertexLabel.end () ; ++it) {
            if (it->second < minsup)
                continue;

            vector<bool> x(TRANS.size(),0);
            for (map<unsigned int, map<unsigned int, unsigned int> >::iterator cur = singleVertex.begin(); cur != singleVertex.end(); ++cur)
                x[cur->first] = (cur->second[it->first]>0);

            DFS_CODE.createSingleVertex(it->first);
            singleTree.push_back(Save(x, it->second, DFS_CODE, NULL));
            DFS_CODE.clear();
        }
        //-----
        cout<<"filename = "<<filename<<endl;
    }

    std::istream &gSpan::read(std::istream &is) {
        Graph g(directed);
        while (true) {
            g.read(is);
            if (g.empty()) break;
            TRANS.push_back(g);
        }
        return is;
    }

    std::map<unsigned int, unsigned int> gSpan::support_counts(Projected &projected) {
        std::map<unsigned int, unsigned int> counts;

        for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur)
            counts[cur->id] += 1;

        return counts;
    }

    unsigned int gSpan::support(Projected &projected) {
        unsigned int oid = 0xffffffff, size = 0;

        for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) {
            if (oid != cur->id)
                ++size;
            oid = cur->id;
        }
        return size;
    }

    void gSpan::mining(vecVec& Xt, PatternSupportList& patternSupportList) {
        Xt.clear();
        patternSupportList.clear();

        // for single vertex pattern
        for (map<unsigned int, unsigned int>::iterator it = singleVertexLabel.begin () ; it != singleVertexLabel.end () ; ++it) {
            if (it->second < minsup)
                continue;

            DFS_CODE.createSingleVertex(it->first);
            patternSupportList.push_back(pair<DFSCode, int>(DFS_CODE, it->second));
            DFS_CODE.clear();

            vector<bool> x(TRANS.size());
            for (map<unsigned int, map<unsigned int, unsigned int> >::iterator cur = singleVertex.begin(); cur != singleVertex.end(); ++cur)
                x[cur->first] = (cur->second[it->first]>0);
            Xt.push_back(x);
        }
        //-----

        EdgeList edges;
        Projected_map3 root;

        for (unsigned int id = 0; id < TRANS.size(); ++id) {
            Graph &g = TRANS[id];
            for (unsigned int from = 0; from < g.size(); ++from) {
                if (get_forward_root(g, g[from], edges)) {
                    for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                        root[g[from].label][(*it)->elabel][g[(*it)->to].label].push(id, *it, 0);
                }
            }
        }

        for (Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end(); ++fromlabel) {
            for (Projected_iterator2 elabel = fromlabel->second.begin(); elabel != fromlabel->second.end(); ++elabel) {
                for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                    /* Build the initial two-node graph.  It will be grown
                     * recursively within project.
                     */
                    DFS_CODE.push(0, 1, fromlabel->first, elabel->first, tolabel->first);
                    __mining(tolabel->second, Xt, patternSupportList);
                    DFS_CODE.pop();
                }
            }
        }
    }

    /* Recursive subgraph mining function (similar to subprocedure 1
     * Subgraph_Mining in [Yan2002]).
     */
    void gSpan::__mining(Projected &projected, vecVec& Xt, PatternSupportList& patternSupportList) {
        /* Check if the pattern is frequent enough.
         */
        unsigned int sup = support(projected);
        if (sup < minsup || DFS_CODE.nodeCount() >= maxpat)
            return;
        /* The minimal DFS code check is more expensive than the support check,
         * hence it is done now, after checking the support.
         */
        if (!is_min())
            return;

        patternSupportList.push_back(pair<DFSCode, int>(DFS_CODE, sup));

        vector<bool> x(TRANS.size());
        for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur)
            x[cur->id] = 1;
        Xt.push_back(x);

        /* We just outputted a frequent subgraph.  As it is frequent enough, so
         * might be its (n+1)-extension-graphs, hence we enumerate them all.
         */
        const RMPath &rmpath = DFS_CODE.buildRMPath();
        int minlabel = DFS_CODE[0].fromlabel;
        int maxtoc = DFS_CODE[rmpath[0]].to;

        Projected_map3 new_fwd_root;
        Projected_map2 new_bck_root;
        EdgeList edges;

        /* Enumerate all possible one edge extensions of the current substructure.
         */
        for (unsigned int n = 0; n < projected.size(); ++n) {

            unsigned int id = projected[n].id;
            PDFS *cur = &projected[n];
            History history(TRANS[id], cur);

            // XXX: do we have to change something here for directed edges?

            // backward
            for (int i = (int) rmpath.size() - 1; i >= 1; --i) {
                Edge *e = get_backward(TRANS[id], history[rmpath[i]], history[rmpath[0]], history);
                if (e)
                    new_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push(id, e, cur);
            }

            // pure forward
            // FIXME: here we pass a too large e->to (== history[rmpath[0]]->to
            // into get_forward_pure, such that the assertion fails.
            //
            // The problem is:
            // history[rmpath[0]]->to > TRANS[id].size()
            if (get_forward_pure(TRANS[id], history[rmpath[0]], minlabel, history, edges))
                for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                    new_fwd_root[maxtoc][(*it)->elabel][TRANS[id][(*it)->to].label].push(id, *it, cur);

            // backtracked forward
            for (int i = 0; i < (int) rmpath.size(); ++i)
                if (get_forward_rmpath(TRANS[id], history[rmpath[i]], minlabel, history, edges))
                    for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                        new_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][TRANS[id][(*it)->to].label].push(id, *it, cur);
        }

        /* Test all extended substructures.
         */
        // backward
        for (Projected_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
            for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
                DFS_CODE.push(maxtoc, to->first, -1, elabel->first, -1);
                __mining(elabel->second, Xt, patternSupportList);
                DFS_CODE.pop();
            }
        }

        // forward
        for (Projected_riterator3 from = new_fwd_root.rbegin(); from != new_fwd_root.rend(); ++from) {
            for (Projected_iterator2 elabel = from->second.begin(); elabel != from->second.end(); ++elabel) {
                for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                    DFS_CODE.push(from->first, maxtoc + 1, -1, elabel->first, tolabel->first);
                    __mining(tolabel->second, Xt, patternSupportList);
                    DFS_CODE.pop();
                }
            }
        }
        return;
    }

    tree<gSpan::Save>::iterator gSpan::createRoot(){
        if (Tree.empty()){
            DFS_CODE.clear();
            EdgeList edges;
            tree<Save>::iterator root=Tree.insert(Tree.begin(), Save());
            
            Projected_map3 &Root=root->new_fwd_root;
            for (unsigned int id = 0; id < TRANS.size(); ++id) {
                Graph &g = TRANS[id];
                for (unsigned int from = 0; from < g.size(); ++from) {
                    if (get_forward_root(g, g[from], edges)) {
                        for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                            Root[g[from].label][(*it)->elabel][g[(*it)->to].label].push(id, *it, 0);
                    }
                }
            }
            root->nextCheck=true;

            for (Projected_iterator3 fromlabel = Root.begin(); fromlabel != Root.end(); ++fromlabel) {
                for (Projected_iterator2 elabel = fromlabel->second.begin(); elabel != fromlabel->second.end(); ++elabel) {
                    for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                        DFS_CODE.push(0, 1, fromlabel->first, elabel->first, tolabel->first);

                        Projected *projected=&(tolabel->second);
                        vector<bool> x(TRANS.size(), 0);
                        unsigned int sup = 0;
                        for (Projected::iterator cur = projected->begin(); cur != projected->end(); ++cur)
                            if (!x[cur->id])
                                ++sup, x[cur->id] = 1;
                        if (sup >= minsup && DFS_CODE.nodeCount() <= maxpat && is_min())
                            Tree.append_child(root, Save(x,sup,DFS_CODE,projected));

                        DFS_CODE.pop();
                    }
                }
            }
        }
        return Tree.begin();
    }
    
    void gSpan::createChildren(const tree<gSpan::Save>::iterator &node){
        if (!node->nextCheck){
            DFS_CODE=node->dfscode;
            
            const RMPath &rmpath = DFS_CODE.buildRMPath();
            int minlabel = DFS_CODE[0].fromlabel;
            int maxtoc = DFS_CODE[rmpath[0]].to;
            
            EdgeList edges;
            for (unsigned int n = 0; n < node->projected->size(); ++n) {
                unsigned int id = (*node->projected)[n].id;
                PDFS *cur = &(*node->projected)[n];
                History history(TRANS[id], cur);

                for (int i = (int) rmpath.size() - 1; i >= 1; --i) {
                    Edge *e = get_backward(TRANS[id], history[rmpath[i]], history[rmpath[0]], history);
                    if (e)
                        node->new_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push(id, e, cur);
                }

                if (get_forward_pure(TRANS[id], history[rmpath[0]], minlabel, history, edges))
                    for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                        node->new_fwd_root[maxtoc][(*it)->elabel][TRANS[id][(*it)->to].label].push(id, *it, cur);

                for (int i = 0; i < (int) rmpath.size(); ++i)
                    if (get_forward_rmpath(TRANS[id], history[rmpath[i]], minlabel, history, edges))
                        for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                            node->new_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][TRANS[id][(*it)->to].label].push(id, *it, cur);
            }
            
            node->nextCheck=true;
            for (Projected_iterator2 to = node->new_bck_root.begin(); to != node->new_bck_root.end(); ++to) {
                for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
                    DFS_CODE.push(maxtoc, to->first, -1, elabel->first, -1);

                    Projected *projected=&(elabel->second);
                    vector<bool> x(TRANS.size(), 0);
                    unsigned int sup = 0;
                    for (Projected::iterator cur = projected->begin(); cur != projected->end(); ++cur)
                        if (!x[cur->id])
                            ++sup, x[cur->id] = 1;
                    if (sup >= minsup && DFS_CODE.nodeCount() <= maxpat && is_min())
                        Tree.append_child(node, Save(x,sup,DFS_CODE,projected));

                    DFS_CODE.pop();
                }
            }

            for (Projected_riterator3 from = node->new_fwd_root.rbegin(); from != node->new_fwd_root.rend(); ++from) {
                for (Projected_iterator2 elabel = from->second.begin(); elabel != from->second.end(); ++elabel) {
                    for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
                        DFS_CODE.push(from->first, maxtoc + 1, -1, elabel->first, tolabel->first);

                        Projected *projected=&(tolabel->second);
                        vector<bool> x(TRANS.size(), 0);
                        unsigned int sup = 0;
                        for (Projected::iterator cur = projected->begin(); cur != projected->end(); ++cur)
                            if (!x[cur->id])
                                ++sup, x[cur->id] = 1;
                        if (sup >= minsup && DFS_CODE.nodeCount() <= maxpat && is_min())
                            Tree.append_child(node, Save(x,sup,DFS_CODE,projected));

                        DFS_CODE.pop();
                    }
                }
            }
        }
    }

    double gSpan::getMaxValue(const std::vector<double>& v) {
        visit=0;
        double maxval = 0;
        // for single vertex pattern
        for (const Save &node : singleTree) {
            ++visit;

            double dot = 0;
            for (int i = 0; i < TRANS.size(); ++i) 
                dot+=v[i]*node.x[i];
            maxval = fmax(maxval, dot);
        }
        //-----

        tree<Save>::iterator root=createRoot();

        for (tree<Save>::sibling_iterator it=Tree.begin(root);it!=Tree.end(root);++it)
            __getMaxValue(it, maxval, v);

        cout<<"getMaxValue:visit = "<<visit<<endl;
        return maxval;
    }

    void gSpan::__getMaxValue(const tree<Save>::iterator &node, double& maxval, const std::vector<double>& v) {
        ++visit;
        
        double p = 0, m = 0;
        for (size_t i=0;i<TRANS.size();++i) 
            (v[i] > 0) ? p += v[i]*node->x[i] : m += v[i]*node->x[i];
        if (fmax(p, -m) < maxval)
            return;
        double dot = fabs(p + m);
        maxval = fmax(maxval, dot);

        createChildren(node);

        for (tree<Save>::sibling_iterator it=Tree.begin(node);it!=Tree.end(node);++it)
            __getMaxValue(it, maxval, v);
    }

    void gSpan::safePatternPruning(double r, const std::vector<double> &theta, vecVec& Xt, PatternSupportList& patternSupportList) {
        Xt.clear();
        patternSupportList.clear();
        visit=0;

        // for single vertex pattern
        for (const Save &node : singleTree) {
            ++visit;

            double dot = 0;
            for (size_t i=0;i<TRANS.size();++i)
                dot += theta[i]*node.x[i];
            double n = TRANS.size();
            int sup=node.support;
            double ub = fabs(dot) + r * sqrt(sup - sup * sup / n);
            if (ub >= 1) {
                Xt.push_back(node.x);
                patternSupportList.push_back(pair<DFSCode, int>(node.dfscode, sup));
            }
        }
        //-----

        tree<Save>::iterator root=createRoot();

        for (tree<Save>::sibling_iterator it=Tree.begin(root);it!=Tree.end(root);++it)
            __safePatternPruning(it, r, theta, Xt, patternSupportList);
        cout<<"safePatternPruning:visit = "<<visit<<endl;
    }

    void gSpan::__safePatternPruning(const tree<Save>::iterator &node, double r, const std::vector<double> &theta, vecVec& Xt, PatternSupportList& patternSupportList) {
        ++visit;

        double p = 0, m = 0;
        for (size_t i=0;i<TRANS.size();++i)
            (theta[i] > 0) ? p += theta[i]*node->x[i] : m += theta[i]*node->x[i];
        double n = TRANS.size();
        int sup=node->support;
        double sppc = fmax(p, -m) + r * sqrt(sup);
        if (sppc < 1)
            return;
        double ub = fabs(p + m) + r * sqrt(sup - sup * sup / n);
        if (ub >= 1) {
            Xt.push_back(node->x);
            patternSupportList.push_back(pair<DFSCode, int>(node->dfscode, sup));
        }
        
        createChildren(node);
        
        for (tree<Save>::sibling_iterator it=Tree.begin(node);it!=Tree.end(node);++it)
            __safePatternPruning(it, r, theta, Xt, patternSupportList);
    }

    void gSpan::regularizationPath(int splitNum, double R, int maxloop, double eps, int freq) {
        int n = TRANS.size();
        vector<double> y(n);
        for (int i = 0; i < n; ++i)
            y[i] = TRANS[i].y;
        vector<double> w, error(y);
        
        // b=mean(y)
        double b = 0;
        for (const double& yi : y)
            b += yi;
        b /= n;
        
        // error = y - b
        for (double& e : error)
            e -= b;
        
        // lambda_max
        double lam = getMaxValue(error);
        cout<<"lambda_max = "<<lam<<endl;
        
        // dual variable
        vector<double> theta(error);
        for (double& t : theta)
            t /= lam;
        
        // pattern -> weight
        map<string, double> wSaveDict;
        double bSave = b;
        
        // Xt: patternNum x sampleSize, 0/1 matrix
        vecVec Xt;
        PatternSupportList pat;
        
        bool first = true;
        double wl1norm = 0;
        int dim = 0;
        for (int index = 1; index < splitNum; ++index) {
            lam *= R;
            cout << "-------------------" << endl << "[" << index << "] lam = " << lam << endl;
            clock_t start = clock();
            for (int loop = 0; loop < maxloop; ++loop) {
                if (!first) {
                    double maxval = 0;
                    for (int i = 0; i < dim; ++i) {
                        double tmp = fabs(error * Xt[i]);
                        if (tmp > maxval)
                            maxval = tmp;
                    }
                    for (int i = 0; i < n; ++i)
                        theta[i] = error[i] / maxval;
                }
                first = false;
                double primal = 0.5 * (error * error) + lam*wl1norm;
                double dual = -0.5 * lam * lam * (theta * theta) + lam * (y * theta);
                double gap = primal - dual;
                cout << "primal = "<<primal<<", dual = "<<dual<<", gap = " << gap << endl;
                if (gap <= eps * primal) {// convergence
                    wSaveDict.clear();
                    for (int i = 0; i < dim; ++i)
                        wSaveDict[pat[i].first.toString()] = w[i];
                    bSave = b;
                    cout << "w = " << w << "b = " << b << endl << "pat = " << pat;
                    break;
                }
                if (loop % freq == 0) {
                    double r = sqrt(2 * gap) / lam;
                    if (loop == 0) {
                        safePatternPruning(r, theta, Xt, pat);
                        w.clear();
                        for (const pair< DFSCode, int> &it : pat)
                            w.push_back(wSaveDict[it.first.toString()]);
                        b = bSave;
                    } else {// dynamic screening
                        PatternSupportList newpat;
                        vecVec newXt;
                        vector<double> new_w;
                        for (int i = 0; i < dim; ++i) {
                            double ub = fabs((Xt[i] * theta)) + r * sqrt(pat[i].second - pat[i].second * pat[i].second / (double) n);
                            if (ub >= 1) {
                                newXt.push_back(Xt[i]);
                                newpat.push_back(pat[i]);
                                new_w.push_back(w[i]);
                            }
                        }
                        Xt = newXt;
                        pat = newpat;
                        w = new_w;
                    }
                    dim = Xt.size();
                    cout << "dim = " << dim << endl;
                }
                wl1norm = 0;
                vector<size_t> randomIdx=randomIndex(dim);
                // update w by shooting algorithm
                for (int rand_idx = 0; rand_idx < dim; ++rand_idx) {
                    size_t k = randomIdx[rand_idx];
                    double yk = 0;
                    for (int i = 0; i < n; ++i) {
                        if (Xt[k][i] == 0)
                            continue;
                        double tmp = y[i] - b;
                        for (int j = 0; j < dim; ++j) {
                            if (j == k)
                                continue;
                            tmp -= Xt[j][i] * w[j];
                        }
                        yk += tmp * Xt[k][i];
                    }
                    if (lam < yk)
                        w[k] = (yk - lam) / pat[k].second;
                    else if (yk < -lam)
                        w[k] = (yk + lam) / pat[k].second;
                    else
                        w[k] = 0;
                    wl1norm += fabs(w[k]);
                }
                // update b: b=mean(y-Xw)
                b = 0;
                for (int i = 0; i < n; ++i) {
                    b += y[i];
                    for (int j = 0; j < dim; ++j)
                        b -= Xt[j][i] * w[j];
                }
                b /= n;
                // error = y - (Xw+b)
                for (int i = 0; i < n; ++i) {
                    error[i] = y[i] - b;
                    for (int j = 0; j < dim; ++j)
                        error[i] -= Xt[j][i] * w[j];
                }
            }
            cout << "time = " << (double) (clock() - start) / CLOCKS_PER_SEC << endl;
        }
    }

    std::ostream& operator<<(std::ostream& os, const PatternSupportList & patternSupportList) {
        os << "[";
        for (const std::pair<DFSCode, int>& p : patternSupportList)
            os << "([ " << p.first << " ], " << p.second << "), ";
        os << "]" << std::endl;
        return os;
    }

}
