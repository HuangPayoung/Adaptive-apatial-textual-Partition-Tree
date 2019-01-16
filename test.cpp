#include "aptree.h"
#include <windows.h>
#include <iomanip>
#include <assert.h>

using namespace std;

void functionTest()
{
    aptree<int> ex(10, 100, 5, 3, 3);
//buildIndex
    vector<aptree<int>::query> Q;
    aptree<int>::query Qtmp;
    ifstream fin("query1000.txt");
    assert(fin);
    while(fin>>Qtmp)
        {Q.push_back(Qtmp);}
    ex.buildIndex(Q);
    fin.close();
    vector<map<int, aptree<int>::query>::iterator> R;
//objectmatching
    fin.open("object.txt");
    aptree<int>::object oex;
    for(int i=0; i<10; ++i)
    {
        fin>>oex;
        cout << ex.ObjectMatching(oex, R) << endl;
        for(int i=0; i<int(R.size()); ++i) cout << R[i]->second<<endl;
        R.clear();
        cout<<endl;
    }
    fin.close();

//regis与deregis
    fin.open("query.txt");
    assert(fin);
    aptree<int>::query ex1;
    fin>>ex1;
    fin.close();

    cout<<"before registering:"<<endl;
    fin.open("object1.txt");
    assert(fin);
    aptree<int>::object oex1;
    fin>>oex1;
    fin.close();
    cout<<ex.ObjectMatching(oex1, R)<<"个结果"<<endl;
    for(int i=0; i<int(R.size()); ++i) cout << R[i]->second<<endl;
    R.clear();

    cout<<"registering ex1, which is:"<<endl;
    cout<<ex1<<endl;
    ex.regis(ex1);
    cout<<"registeration finished"<<endl;

    cout<<"checking, predicted result:"<<endl;
    cout<<ex1<<endl;

    R.clear();
    ex.ObjectMatching(oex1, R);
    cout<<"the result is:"<<endl;
    for(int i=0; i<int(R.size()); ++i) cout << R[i]->second<<endl;
    cout<<endl;
    R.clear();

    cout<<"before deregistering"<<endl;
    fin.open("object1.txt");
    assert(fin);
    fin>>oex1;
    fin.close();
    cout<<ex.ObjectMatching(oex1, R);
    R.clear();
    cout<<endl;

    cout<<"deregistering ex1..."<<endl;
    ex.deregis(1000);   //该用户id是1000
    cout<<"checking, predicted result:"<<endl;
    cout<<"0"<<endl;

    cout<<"the result is:"<<endl;
    cout<<ex.ObjectMatching(oex1, R);
    R.clear();
    cout<<endl;

//多次regis与derigis
    cout<<"将query1000-10000.txt中的前100个query依次注册注销,然后进行匹配"<<endl;
    fin.open("query1000-10000.txt");
    assert(fin);
    for(int i=0; i<100; ++i)
    {
        fin>>Qtmp;
        ex.regis(Qtmp);
    }
    fin.close();
//ex.print("aptree1RE.txt");
    cout<<"after registering"<<endl;
    fin.open("object.txt");
    for(int i=0; i<10; ++i)
    {
        fin>>oex;
        cout << ex.ObjectMatching(oex, R) << endl;
        for(int i=0; i<int(R.size()); ++i) cout << R[i]->second<<endl;
        R.clear();
        cout<<endl;
    }
    fin.close();
    for(int id=1001; id<=1100; ++id)
    {
        ex.deregis(id);
    }

    cout<<"after deregistering"<<endl;
    fin.open("object.txt");
    for(int i=0; i<10; ++i)
    {
        fin>>oex;
        cout << ex.ObjectMatching(oex, R) << endl;
        for(int i=0; i<int(R.size()); ++i) cout << R[i]->second<<endl;
        R.clear();
        cout<<endl;
    }
    fin.close();
}

void build_match(int Size)
{
    LARGE_INTEGER beginTime;
    LARGE_INTEGER endTime;
    LARGE_INTEGER freq;
    double executeTime;

    QueryPerformanceFrequency(&freq);

    cout<<left;
    aptree<int> ex(10, 10, 5, 3, 3);
    vector<aptree<int>::query> Q;
    aptree<int>::query Qtmp;
    ifstream fin("query1000-10000.txt");
    for(int i=0; i<Size; ++i)
    {
        fin>>Qtmp;
        Q.push_back(Qtmp);
    }
    cout<<setw(6)<<Q.size();
    QueryPerformanceCounter(&beginTime);
    ex.buildIndex(Q);
    QueryPerformanceCounter(&endTime);
    executeTime=(double)(endTime.QuadPart-beginTime.QuadPart)/freq.QuadPart*1000;//单位:毫秒
    cout<<setw(16)<<executeTime;

    vector<map<int, aptree<int>::query>::iterator> R;
    fin.close();
    fin.open("object10000.txt");
    aptree<int>::object oex;
    QueryPerformanceCounter(&beginTime);
//objectmatching
    while(fin>>oex)
    {
        ex.ObjectMatching(oex, R);
        R.clear();
    }
    QueryPerformanceCounter(&endTime);
    executeTime=(double)(endTime.QuadPart-beginTime.QuadPart)/freq.QuadPart*1000;//单位:毫秒
    cout<<executeTime<<endl;

    fin.close();
}

void regis_deregis(aptree<int> &ex)
{
    LARGE_INTEGER beginTime;
    LARGE_INTEGER endTime;
    LARGE_INTEGER freq;
    double executeTime;

    QueryPerformanceFrequency(&freq);

    aptree<int>::query Qtmp;
    ifstream fin;

    cout<<left;
    int id=1000;
    fin.open("query1000-10000.txt");
    for(int i=100; i<=1000; i+=100)
    {
        cout<<setw(6)<<i;

        QueryPerformanceCounter(&beginTime);
        for(int j=0; j<i; ++j)
        {
            fin>>Qtmp;
            ex.regis(Qtmp);
        }
        QueryPerformanceCounter(&endTime);
        id+=i;
        executeTime=(double)(endTime.QuadPart-beginTime.QuadPart)/freq.QuadPart*1000;//单位:毫秒
        cout<<setw(14)<<executeTime;

        QueryPerformanceCounter(&beginTime);
        for(int j=0; j<i; ++j)
        {
            ex.deregis(id-i+j);
        }
        QueryPerformanceCounter(&endTime);
        executeTime=(double)(endTime.QuadPart-beginTime.QuadPart)/freq.QuadPart*1000;//单位:毫秒
        cout<<executeTime<<endl;
    }
    fin.close();
}

void match(const aptree<int> &ex)
{
    LARGE_INTEGER beginTime;
    LARGE_INTEGER endTime;
    LARGE_INTEGER freq;
    double executeTime;

    QueryPerformanceFrequency(&freq);
    cout<<"query:1000"<<'\t';

    vector<map<int, aptree<int>::query>::iterator> R;
    ifstream fin("object10000.txt");
    aptree<int>::object oex;
    QueryPerformanceCounter(&beginTime);
//objectmatching
    while(fin>>oex)
    {
        ex.ObjectMatching(oex, R);
        R.clear();
    }
    QueryPerformanceCounter(&endTime);
    executeTime=(double)(endTime.QuadPart-beginTime.QuadPart)/freq.QuadPart*1000;//单位:毫秒
    cout<<"time(10000 times)(ms):"<<executeTime<<endl;

    fin.close();
}

void timeTest()
{
    cout<<left;
    cout<<"buildIndex objectMatching time:"<<endl;
    cout<<setw(6)<<"size"<<setw(16)<<"buildIndex(ms)"<<"objectMatching(10000 times)(ms)"<<endl;
    for(int i=100; i<=1000; i+=100)
        build_match(i);

    cout<<endl<<"register deregister time:"<<endl;
    cout<<setw(6)<<"size"<<setw(14)<<"register(ms)"<<"deregister(ms)"<<endl;
    aptree<int> ex(10, 5000, 5, 3, 3);
    vector<aptree<int>::query> Q;
    aptree<int>::query Qtmp;
    ifstream fin("query1000.txt");
    while(fin>>Qtmp)
        {Q.push_back(Qtmp);}
    ex.buildIndex(Q);
    fin.close();

    regis_deregis(ex);

    cout<<endl<<"after register and derigister, objectmatching:"<<endl;
    match(ex);
}
