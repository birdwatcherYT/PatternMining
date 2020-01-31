#ifndef SEQUENCEHUBER_HPP
#define SEQUENCEHUBER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <map>
#include <vector>
#include <cmath>
#include <ctime>
#include <numeric>
#include <climits>
#include <cfloat>
#include "tree.hh"

typedef std::vector< std::pair<int, int> > ProjectDB;
typedef std::vector< std::vector<bool> > vecVec;
typedef std::vector<int> Pattern;
typedef std::vector< std::pair< Pattern, int> > PatternSupportList;

std::ostream& operator<<(std::ostream& os, const PatternSupportList& patternSupportList);

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector< std::vector<T> >& mat) {
    os << "[" << std::endl;
    for (const std::vector<T>& v : mat)
        os << v;
    os << "]" << std::endl;
    return os;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> vec) {
    os << "[ ";
    for (const T &e : vec)
        os << e << " ";
    os << "]" << std::endl;
    return os;
}

template<class T1, class T2>
double operator*(const std::vector<T1>& v1, const std::vector<T2> v2) {
    double dot = 0;
    for (int i = 0, size = v1.size(); i < size; ++i)
        dot += v1[i] * v2[i];
    return dot;
}

std::vector<size_t> randomIndex(size_t n) ;

class SequenceMining {
private:
    std::vector< Pattern > transaction;
    int n;
    int minsup;
    int maxpat;
    std::vector<double> y;
    
    class Save {
    public:
        std::vector<bool> x;
        unsigned int support;
        Pattern pattern;
        ProjectDB pdb;
        bool nextCheck=false;
        Save(){}
        Save(const std::vector<bool> &x, unsigned int support, const Pattern &pattern, const ProjectDB &pdb){
            this->x=x, this->support=support, this->pattern=pattern, this->pdb=pdb;
        }
    };
    tree<Save> Tree;
    int visit;

    void __mining(const ProjectDB& pdb, vecVec& Xt, std::vector<int>& pattern, std::vector<int> &supportList, PatternSupportList& patternSupportList);
    void __getMaxValue(const tree<Save>::iterator &node, double &maxval, const std::vector<double>& v);
    void __safePatternPruning(const tree<Save>::iterator &node, double r, const std::vector<double> &theta, vecVec& Xt, PatternSupportList& patternSupportList);
    tree<Save>::iterator createRoot();
    void createChildren(const tree<Save>::iterator &node);
public:
    SequenceMining(const std::string& filename, int minsup, int maxpat);
    void mining(vecVec& Xt, PatternSupportList& patternSupportList);
    double getMaxValue(const std::vector<double>& v);
    void safePatternPruning(double r, const std::vector<double>& theta, vecVec& Xt, PatternSupportList& patternSupportList);
    void regularizationPath(int splitNum, double R, int maxloop, double eps, int freq);
};

#endif /* SEQUENCEHUBER_HPP */

