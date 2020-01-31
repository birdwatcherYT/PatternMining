#include "ItemsetMining.hpp"

using namespace std;

std::ostream& operator<<(std::ostream& os, const PatternSupportList& patternSupportList) {
    os << "[";
    for (const pair<Pattern, int>& p : patternSupportList) {
        os << "([ ";
        for (const int& element : p.first)
            os << element << " ";
        os << "], " << p.second << "), ";
    }
    os << "]" << endl;
    return os;
}

ItemsetMining::ItemsetMining(const string &filename, int minsup, int maxpat) {
    this->minsup = minsup;
    this->maxpat = maxpat;

    ifstream ifs(filename);
    if (ifs.fail()) {
        cerr << "ファイルオープンに失敗" << endl;
        exit(1);
    }
    //行ごと
    string line;
    double value;
    while (getline(ifs, line)) {
        stringstream ss(line);
        vector<int> items;
        ss >> value;
        y.push_back(value);
        int item;
        while (ss >> item)
            items.push_back(item);
        transaction.push_back(items);
    }
    n = y.size();
}

void ItemsetMining::mining(vecVec& Xt, PatternSupportList& patternSupportList) {
    ProjectDB pdb;
    for (int i = 0; i < n; ++i)
        pdb.push_back(pair<int, int>(i, -1));
    Xt.clear();
    patternSupportList.clear();

    vector<int> pattern;
    vector<int> supportList;
    __mining(pdb, Xt, pattern, supportList, patternSupportList);
}

void ItemsetMining::__mining(const ProjectDB& pdb, vecVec& Xt, std::vector<int>& pattern, std::vector<int> &supportList, PatternSupportList& patternSupportList) {
    if (pattern.size() >= maxpat)
        return;
    map<int, ProjectDB> counter;
    for (const pair<int, int>& p : pdb) {
        int id = p.first;
        for (int j = p.second + 1; j < transaction[id].size(); ++j)
            counter[transaction[id][j]].push_back(pair<int, int>(id, j));
    }
    for (const auto& it : counter) {
        vector<bool> x(n, 0);
        int support = 0, oldID = -1;
        for (const pair<int, int>& p : it.second) {
            int id = p.first;
            if (oldID != id) {
                ++support;
                x[id] = 1;
                oldID = id;
            }
        }
        if (support < minsup)
            continue;
        pattern.push_back(it.first);
        supportList.push_back(support);
        Xt.push_back(x);
        patternSupportList.push_back(pair<vector<int>, int>(pattern, support));

        __mining(it.second, Xt, pattern, supportList, patternSupportList);

        pattern.pop_back();
        supportList.pop_back();
    }
}

tree<ItemsetMining::Save>::iterator ItemsetMining::createRoot(){
    if (Tree.empty()){
        ProjectDB pdb;
        for (int i = 0; i < n; ++i)
            pdb.push_back(pair<int, int>(i, -1));
        tree<Save>::iterator root=Tree.insert(Tree.begin(), Save(vector<bool>(), 0, Pattern(), pdb));
        createChildren(root);
    }
    return Tree.begin();
}

void ItemsetMining::createChildren(const tree<Save>::iterator &node){
    if (!node->nextCheck){
        node->nextCheck=true;
        if (node->pattern.size() >= maxpat)
            return;
        map<int, ProjectDB> counter;
        for (const pair<int, int>& p : node->pdb) {
            int id = p.first;
            for (int j = p.second + 1; j < transaction[id].size(); ++j)
                counter[transaction[id][j]].push_back(pair<int, int>(id, j));
        }
        Pattern pattern(node->pattern);
        for (const auto& it : counter) {
            vector<bool> x(n, 0);
            int support = 0, oldID = -1;
            for (const pair<int, int>& p : it.second) {
                int id = p.first;
                if (oldID != id) {
                    ++support;
                    x[id] = 1;
                    oldID = id;
                }
            }
            if (support < minsup)
                continue;
            pattern.push_back(it.first);
            Tree.append_child(node, Save(x, support, pattern, it.second));
            pattern.pop_back();
        }
    }
}

double ItemsetMining::getMaxValue(const std::vector<double>& v) {
    visit=0;
    double maxval = 0;
    tree<Save>::iterator root=createRoot();
    for (tree<Save>::sibling_iterator it=Tree.begin(root);it!=Tree.end(root);++it)
        __getMaxValue(it, maxval, v);
    cout<<"getMaxValue:visit = "<<visit<<endl;
    return maxval;
}

void ItemsetMining::__getMaxValue(const tree<Save>::iterator &node, double &maxval, const std::vector<double>& v) {
    ++visit;

    double p = 0, m = 0;
    for (size_t i=0;i<n;++i) 
        (v[i] > 0) ? p += v[i]*node->x[i] : m += v[i]*node->x[i];
    if (fmax(p, -m) < maxval)
        return ;
    double dot = fabs(p + m);
    maxval = fmax(maxval, dot);

    createChildren(node);

    for (tree<Save>::sibling_iterator it=Tree.begin(node);it!=Tree.end(node);++it)
        __getMaxValue(it, maxval, v);
}

void ItemsetMining::safePatternPruning(double r, const std::vector<double> &theta, vecVec& Xt, PatternSupportList& patternSupportList) {
    Xt.clear();
    patternSupportList.clear();
    visit=0;

    tree<Save>::iterator root=createRoot();

    for (tree<Save>::sibling_iterator it=Tree.begin(root);it!=Tree.end(root);++it)
        __safePatternPruning(it, r, theta, Xt, patternSupportList);
    cout<<"safePatternPruning:visit = "<<visit<<endl;
}

void ItemsetMining::__safePatternPruning(const tree<Save>::iterator &node, double r, const std::vector<double> &theta, vecVec& Xt, PatternSupportList& patternSupportList) {
    ++visit;

    double p = 0, m = 0;
    for (size_t i=0;i<n;++i)
        (theta[i] > 0) ? p += theta[i]*node->x[i] : m += theta[i]*node->x[i];
    int sup=node->support;
    double sppc = fmax(p, -m) + r * sqrt(sup);
    if (sppc < 1)
        return;
    double ub = fabs(p + m) + r * sqrt(sup - sup * sup /(double) n);
    if (ub >= 1) {
        Xt.push_back(node->x);
        patternSupportList.push_back(pair<Pattern, int>(node->pattern, sup));
    }
    
    createChildren(node);
    
    for (tree<Save>::sibling_iterator it=Tree.begin(node);it!=Tree.end(node);++it)
        __safePatternPruning(it, r, theta, Xt, patternSupportList);
}

void ItemsetMining::regularizationPath(int splitNum, double R, int maxloop, double eps, int freq) {
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
    map<Pattern, double> wSaveDict;
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
                    wSaveDict[pat[i].first] = w[i];
                bSave = b;
                cout << "w = " << w << "b = " << b << endl << "pat = " << pat;
                break;
            }
            if (loop % freq == 0) {
                double r = sqrt(2 * gap) / lam;
                if (loop == 0) {
                    safePatternPruning(r, theta, Xt, pat);
                    w.clear();
                    for (const pair< Pattern, int> &it : pat)
                        w.push_back(wSaveDict[it.first]);
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

std::vector<size_t> randomIndex(size_t n){
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::random_shuffle(idx.begin(), idx.end());
    return idx;
}
