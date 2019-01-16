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
        double x, y;   //��������
    };
    struct space
    {
        double x1, x2;  //�����귶Χx1<x2
        double y1, y2;  //�����귶Χy1<y2
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
    static const int K=1, S=2, Q=3; //������ʾ������
    static const int Inf=2147483647;

    bool LinS(const location &l, const space &s)const   //�жϵ�l�Ƿ���������s
    {
        if(l.x>=s.x1&&l.x<=s.x2 && l.y<=s.y2&&l.y>=s.y1) return true;
        else return false;
    }
    bool overlap(const space &s1, const space &s2)const //�ж�s1��s2�Ƿ��н���
    {
        if(s1.x1>=s2.x2 || s1.x2<=s2.x1) return false;
        if(s1.y1>=s2.y2 || s1.y2<=s2.y1) return false;
        return true;
    }
    bool cover(const space &s1,const space &s2)const    //�ж�s1�Ƿ���ȫ��s2����
    {
        if (s1.x1>=s2.x1&&s1.x2<=s2.x2&&s1.y1>=s2.y1&&s1.y2<=s2.y2) return true;
        else return false;
    }
    struct node
    {
        int property;   //ȡֵΪK,S,Q,�ֱ��ʾk-node,s-node,q-node
    };
    struct kNode:public node
    {
        int offset;                 //�ý��ƥ��ڼ����ؼ���
        vector<vector<kType> > cut;
        //ÿ��Ԫ�ض�Ӧһ��cut
        //ÿ��Ԫ��Ϊһ��ָ�����Ӧ���йؼ�����֯�ɵ����Ա�
        vector<node*> N;   //ָ������ӽڵ�,����dummy cut�Ĵ���ͬsNode
    };
    struct sNode:public node
    {
        space s;                    //��ǰ����Ŀռ�
        vector<space> cell;         //ÿ��Ԫ�ض�Ӧһ���ռ�����
        vector<node*> N;
        //ָ������ӽڵ�
        //dummy cell�������
        //��N.size()>cell.size()�ж��Ƿ����dummy cell
        //��N.back()����dummy cell
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
        vector<kType> key;  //������йؼ���
        space reg;      //��Ӧ�Ŀռ䷶Χ
    private:
        vector<node*> N; //���а����� query �� q-node
    };
private:
    struct qNode:public node
    {
        vector<typename map<int, query>::iterator> q;           //ָ������ĸ������,Ԫ��Ϊ������
    };

    int thre;   //ÿ��q-node��������query��
    double KL;  //query���ص�KL-divergence����,ʵ��Ϊע��ע����������
    int kf;     //ÿ��k-node��Ӧ�ķ�֧(cut)��
    int sfm, sfn;     //ÿ��s-node����Ӧ���򻮷�Ϊ sfm*sfn�ľ�������
    int currentId;
    int timesChange;    //ע��+ע����Ŀ

//    int total;          //�ؼ�����Ƶ��
    space SpaceOfTree;
    vector<kType> KeyOfTree;
    map<kType, int> frequency;   //�����ؼ��ֵ�Ƶ��

    node* root; //�����
    map<int, query> Qmap;  //�����е�query��֯��STL�еĶ�̬���ұ�Qmap�У�int��Ӧ�ĳ�Ա������Ϊ���û���id

    void ObjectMatching(const object &o, int ita, node *n, vector<typename map<int, query>::iterator> &R)const;
    //n: �ӽ��n��ʼ
    //ita: ��ǰƥ���ita��keyword
    int cFind(const vector<vector<kType> > &key, const kType &k, int s, int e)const;
    //��kN�ж��ֲ�����k��ƥ���cut�����ض�Ӧ���±�
    //��ʼ�±�s����ֹ�±�e��������
    //δ�ҵ�����-1
    int sFind(const vector<space> &reg, const location &s)const;    //�������
    int sFind(const vector<space> &reg, const location &s, int st, int e)const;
    int kFind(const vector<kType> &key, const kType &k, int s, int e)const;
    void reconstruct();
    void buildIndex(node *&N, const vector<typename map<int, query>::iterator> Qset,
                     int l, const space &s, bool kP, bool sP);
    //N:��ǰ���,Q:query����
    //l:��ǰƥ����ǵڼ����ؼ��֣�kP:���Լ���ͨ��keyword���֣�sP:���Լ���ͨ��space����
    double Kpart(const vector<typename map<int, query>::iterator> &Qset,
               int l, const vector<kType> &V, vector<vector<kType> > &re);
    //���ذ��չؼ��ֻ��ַ�����Ӧ��ʱ�����ĵ�λ
    //�������Qset: Ҫ�����ֵ�query����
    //�������re: ���ֽ������re, re��ÿ��Ԫ�ش���һ��ؼ���
    //ÿ��kNode��Ӧ�ĸ�����֧�ؼ������򣨵�����
    //l: Ҫ��query�ĵ�l���ؼ��ֻ���
    double Spart(const vector<typename map<int, query>::iterator> &Qset,
               const space &reg, vector<space> &re);
    //���ذ��տռ仮�ַ�����Ӧ��ʱ�����ĵ�λ
    //�������Qset: Ҫ�����ֵ�query����
    //�������re: ���ֽ������re, re��ÿ��Ԫ�ش���һ������
    //������״:sfm*sfn�ľ���
    //ÿ��sNode��Ӧ�ĸ�����֧���򣨺����������ÿ���������Ӧ�������������
    //��ֳ�2*2����[(1,2)(1,2)] [(1,2)(2,3)] [(2,3)(1,2)] [(2,3)(2,3)]
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
    //����ƥ�䵽��query����
    //ƥ�䵽��query�ĵ�ַ���浽�����R����
    //o: ��ƥ���object
    aptree(int thre0, double KL0, int kf0, int sfm0, int sfn0)
        :thre(thre0), KL(KL0), kf(kf0), sfm(sfm0), sfn(sfn0), currentId(0), timesChange(0)
    {
        assert(kf0>1);
        assert(thre0>0);
        assert(KL0>0);
    }
    ~aptree();
    void buildIndex(const vector<query> &Q);
    int regis(const query &q);     //ע���û�q
    void deregis(int ID);           //ע��idΪID���û�q
    void print(const char name[])const;    //����debug
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
