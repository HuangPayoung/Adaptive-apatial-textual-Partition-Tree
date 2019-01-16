#include <iostream>
#include "aptree.h"
#include <assert.h>
using namespace std;

int main()
{
    functionTest();
    timeTest();

//    aptree<int> ex(10, 100, 5, 3, 3);
////buildIndex
//    vector<typename aptree<int>::query> Q;
//    typename aptree<int>::query Qtmp;
//    ifstream fin("query1000-10000.txt");
//    assert(fin);
//    while(fin>>Qtmp)
//        {Q.push_back(Qtmp);}
//    ex.buildIndex(Q);
//    fin.close();
//
//    vector<typename map<int, typename aptree<int>::query>::iterator> R;
////objectmatching
//    fin.open("object.txt");
//    typename aptree<int>::object oex;
//    while(fin>>oex)
//    {
//        cout << ex.ObjectMatching(oex, R) << endl;
//        for(int i=0; i<int(R.size()); ++i) cout << R[i]->second<<endl;
//        R.clear();
//        cout<<endl;
//    }
//    fin.close();
    return 0;
}
