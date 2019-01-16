#ifndef APTREE_H_INCLUDED
#define APTREE_H_INCLUDED

#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <assert.h>

using namespace std;

void functionTest();
void timeTest();

template<class kType>
class aptree
{
friend void functionTest();
public:
    struct location
    {
        double x, y;   //横纵坐标
    };
    struct space
    {
        double x1, x2;  //横坐标范围x1<x2
        double y1, y2;  //纵坐标范围y1<y2
        friend ostream &operator<<(ostream &os, const space &o){os<<o.x1<<'\t'<<o.x2<<'\t'<<o.y1<<'\t'<<o.y2; return os;}
        double area()const{return (x2-x1)*(y2-y1);}
    };
    struct object
    {
        vector<kType> key;
        location loc;
        friend ifstream &operator>>(ifstream &fin, object &o)
        {
            o.key.clear();
            fin >> o.loc.x >> o.loc.y;
            kType kt; int kn;
            fin>>kn;
            for(int i=0; i<kn; ++i){fin>>kt; o.key.push_back(kt);}
            return fin;
        }
        friend ostream &operator<<(ostream &os, const object &o)
        {
            os <<"location:"<<'\t'<< o.loc.x << '\t' << o.loc.y << endl;
            os<<"keys:"<<'\t';
            for(int i=0; i<int(o.key.size()); ++i) os<<'\t'<<o.key[i]<<'\t';
            os<<endl;
            return os;
        }
    };
private:
    static const int K=1, S=2, Q=3; //用来表示结点的类
    static const int Inf=2147483647;

    bool LinS(const location &l, const space &s)const   //判断点l是否落在区域s
    {
        if(l.x>=s.x1&&l.x<=s.x2 && l.y<=s.y2&&l.y>=s.y1) return true;
        else return false;
    }
    bool overlap(const space &s1, const space &s2)const //判断s1与s2是否有交集
    {
        if(s1.x1>=s2.x2 || s1.x2<=s2.x1) return false;
        if(s1.y1>=s2.y2 || s1.y2<=s2.y1) return false;
        return true;
    }
    bool cover(const space &s1,const space &s2)const    //判断s1是否完全被s2覆盖
    {
        if (s1.x1>=s2.x1&&s1.x2<=s2.x2&&s1.y1>=s2.y1&&s1.y2<=s2.y2) return true;
        else return false;
    }
    struct node
    {
        int property;   //取值为K,S,Q,分别表示k-node,s-node,q-node
    };
    struct kNode:public node
    {
        int offset;                 //该结点匹配第几个关键字
        vector<vector<kType> > cut;
        //每个元素对应一个cut
        //每个元素为一个指向其对应所有关键字组织成的线性表
        vector<node*> N;   //指向各个子节点,关于dummy cut的处理同sNode
    };
    struct sNode:public node
    {
        space s;                    //当前分配的空间
        vector<space> cell;         //每个元素对应一个空间区域
        vector<node*> N;
        //指向各个子节点
        //dummy cell放在最后
        //用N.size()>cell.size()判断是否存在dummy cell
        //用N.back()调用dummy cell
    };
public:
    struct query
    {
        friend ifstream &operator>>(ifstream &fin, query &o)
        {
            o.key.clear();
            fin >> o.reg.x1 >> o.reg.x2 >> o.reg.y1 >> o.reg.y2;
            kType kt; int kn;
            fin>>kn;
            for(int i=0; i<kn; ++i)
            {
                fin>>kt;
                o.key.push_back(kt);
            }
            return fin;
        }
        friend ostream &operator<<(ostream &os, const query &q)
        {
            os <<"space:"<<'\t'<< q.reg.x1 << '\t' << q.reg.x2 << '\t' << q.reg.y1 << '\t' << q.reg.y2 << endl;
            os <<"keys:"<<'\t';
            for(int i=0; i<int(q.key.size()); ++i) os<<q.key[i]<<'\t';
            os<<endl;
            return os;
        }
        friend class aptree;
        vector<kType> key;  //存放所有关键字
        space reg;      //对应的空间范围
    private:
        vector<node*> N; //所有包含该 query 的 q-node
    };
private:
    struct qNode:public node
    {
        vector<typename map<int, query>::iterator> q;           //指向包含的各个结点,元素为迭代器
    };

    int thre;   //每个q-node最多包含的query数
    double KL;  //query负载的KL-divergence上限,实现为注册注销数量上限
    int kf;     //每个k-node对应的分支(cut)数
    int sfm, sfn;     //每个s-node将对应区域划分为 sfm*sfn的矩阵区域
    int currentId;
    int timesChange;    //注册+注销数目

//    int total;          //关键字总频数
    space SpaceOfTree;
    vector<kType> KeyOfTree;
    map<kType, int> frequency;   //各个关键字的频数

    node* root; //根结点
    map<int, query> Qmap;  //把所有的query组织在STL中的动态查找表Qmap中，int对应的成员可以认为是用户的id

    void ObjectMatching(const object &o, int ita, node *n, vector<typename map<int, query>::iterator> &R)const;
    //n: 从结点n开始
    //ita: 当前匹配第ita个keyword
    int cFind(const vector<vector<kType> > &key, const kType &k, int s, int e)const;
    //在kN中二分查找与k相匹配的cut，返回对应的下标
    //起始下标s，终止下标e，均包括
    //未找到返回-1
    int sFind(const vector<space> &reg, const location &s)const;    //外包函数
    int sFind(const vector<space> &reg, const location &s, int st, int e)const;
    int kFind(const vector<kType> &key, const kType &k, int s, int e)const;
    void reconstruct();
    void buildIndex(node *&N, const vector<typename map<int, query>::iterator> Qset,
                     int l, const space &s, bool kP, bool sP);
    //N:当前结点,Q:query集合
    //l:当前匹配的是第几个关键字，kP:可以继续通过keyword划分，sP:可以继续通过space划分
    double Kpart(const vector<typename map<int, query>::iterator> &Qset,
               int l, const vector<kType> &V, vector<vector<kType> > &re);
    //返回按照关键字划分方案对应的时间消耗单位
    //传入参数Qset: 要被划分的query集合
    //传入参数re: 划分结果存入re, re的每个元素代表一组关键字
    //每个kNode对应的各个分支关键字有序（递增）
    //l: 要对query的第l个关键字划分
    double Spart(const vector<typename map<int, query>::iterator> &Qset,
               const space &reg, vector<space> &re);
    //返回按照空间划分方案对应的时间消耗单位
    //传入参数Qset: 要被划分的query集合
    //传入参数re: 划分结果存入re, re的每个元素代表一个区域
    //划分形状:sfm*sfn的矩阵
    //每个sNode对应的各个分支有序（横坐标递增，每个横坐标对应的纵坐标递增）
    //如分成2*2区域：[(1,2)(1,2)] [(1,2)(2,3)] [(2,3)(1,2)] [(2,3)(2,3)]
    double calCost(const vector<vector<kType> > &re,
                   const vector<typename map<int, query>::iterator> &Qset,
                   int total, int l)const;
    void Clear()
    {
        frequency.clear();
        KeyOfTree.clear();
        for(typename map<int, query>::iterator pos = Qmap.begin(); pos!=Qmap.end(); ++pos)
        {
            pos->second.N.clear();
        }
        Clear(root);
    }
    void Clear(node *n);
    void regis(query &q,node *N);
//    void deregis(query &q,node *N)
public:
    int ObjectMatching(const object &o, vector<typename map<int, query>::iterator> &R)const
    {ObjectMatching(o, 0, root, R); return int(R.size());}
    //返回匹配到的query个数
    //匹配到的query的地址储存到传入的R参数
    //o: 带匹配的object
    aptree(int thre0, double KL0, int kf0, int sfm0, int sfn0)
        :thre(thre0), KL(KL0), kf(kf0), sfm(sfm0), sfn(sfn0), currentId(0), timesChange(0)
    {
        assert(kf0>1);
        assert(thre0>0);
        assert(KL0>0);
    }
    ~aptree();
    void buildIndex(const vector<query> &Q);
    int regis(const query &q);     //注册用户q
    void deregis(int ID);           //注销id为ID的用户q
    void print(const char name[])const;    //用于debug
};

template<class kType>
aptree<kType>::~aptree()
{
    Clear(reinterpret_cast<node*>(root));
}

#define BUILDINDEX
#define OBJECTMATCHING
#define PRINT
#define TOOL
#define PART
#define REGIS
#include "buildIndex.cpp"
#include "objectMatching.cpp"
#include "print.cpp"
#include "tool.cpp"
#include "part.cpp"
#include "regis.cpp"

#endif // APTREE_H_INCLUDED
