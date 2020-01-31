// LASSO for itemset data
#include "ItemsetMining.hpp"
using namespace std;

const int splitNum = 100;
const int maxloop = INT_MAX;
const double eps = 1e-6;
const int freq = 5;
const double R=pow(0.01, 1.0/(splitNum-1));
const int minsup=1;
const int maxpat=10;
string filepath="./data/";

int main(int argc, char **argv) {
    if (argc!=2){
        cerr<<"usage: ./run filename"<<endl;
        return 0;
    }else{
        filepath+=string(argv[1]);
    }
    // LASSO solve
    ItemsetMining itemset(filepath,minsup, maxpat);
    itemset.regularizationPath(splitNum, R, maxloop, eps, freq);
    
    // all pattern output
    // ItemsetMining itemset(filepath,minsup, maxpat);
    // vecVec Xt;
    // PatternSupportList patternSupportList;
    // itemset.mining(Xt, patternSupportList);
    // cout<<patternSupportList<<endl;
    return 0;
}
