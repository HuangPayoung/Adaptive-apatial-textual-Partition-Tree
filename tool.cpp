#ifdef TOOL

template<class kType>
int aptree<kType>::cFind(const vector<vector<kType> > &key, const kType &k, int s, int e)const
{
    if(k<key[s].front() || k>key[e].back()) return -1;
    if(s==e) return s;
    int mid=(s+e)/2;
    if(k < key[mid].front()) return cFind(key, k, s, mid-1);
    else if(k > key[mid].back()) return cFind(key, k, mid+1, e);
    else return mid;
}

template<class kType>
int aptree<kType>::sFind(const vector<space> &reg, const location &s)const
{
    if(s.x<reg[0].x1 || s.x>reg[reg.size()-1].x2 || s.y<reg[0].y1 || s.y>reg[reg.size()-1].y2) return -1;   //该语句不能放入内函数
    else return sFind(reg, s, 0, reg.size()-1);
}

template<class kType>
int aptree<kType>::sFind(const vector<space> &reg, const location &s, int st, int e)const
{
    if(st==e) return st;
    int mid=(st+e)/2;
    if(s.x<reg[mid].x1 || (s.x<=reg[mid].x2&&s.y<reg[mid].y1))
        return sFind(reg, s, st, mid-1);
    else if(s.x>reg[mid].x2 || (s.x>=reg[mid].x1&&s.y>reg[mid].y2))
        return sFind(reg, s, mid+1, e);
    else return mid;
}

template<class kType>
int aptree<kType>::kFind(const vector<kType> &key, const kType &k, int s, int e)const
{
    if(k<key[s] || k>key[e]) return -1;
    if(s==e) return s;
    int mid=(s+e)/2;
    if(k < key[mid]) return kFind(key, k, s, mid-1);
    else if(k > key[mid]) return kFind(key, k, mid+1, e);
    else return mid;
}

template<class kType>
void aptree<kType>::reconstruct()
{
    Clear();
    typename map<int, query>::iterator p;
    vector<typename map<int, query>::iterator> Q;
    for(p = Qmap.begin(); p!=Qmap.end(); ++p) Q.push_back(p);
    double x1min=Q[0]->second.reg.x1, x2max=Q[0]->second.reg.x2, y1min=Q[0]->second.reg.y1, y2max=Q[0]->second.reg.y2;
    for(int i=0; i<int(Q.size()); ++i)
    {
        if(Q[i]->second.reg.x1<x1min) x1min = Q[i]->second.reg.x1;
        if(Q[i]->second.reg.x2>x2max) x2max = Q[i]->second.reg.x2;
        if(Q[i]->second.reg.y1<y1min) y1min = Q[i]->second.reg.y1;
        if(Q[i]->second.reg.y2>y2max) y2max = Q[i]->second.reg.y2;
    }
    space sp = {x1min, x2max, y1min, y2max};
    SpaceOfTree = sp;
    vector<kType> ktmp;
    typename map<kType, int>::iterator keyPos;
    for (int i=0;i<int(Q.size());i++)
    {
        for (int j=0;j<int(Q[i]->second.key.size());j++)
        {
            keyPos = frequency.find(Q[i]->second.key[j]);
            if(keyPos == frequency.end())//没找到
            {
                frequency.insert(pair<kType, int>(Q[i]->second.key[j],1 ));
                ktmp.push_back(Q[i]->second.key[j]);
            }
            else
            {
                ++keyPos->second;
            }
        }
    }
    sort(ktmp.begin(),ktmp.end());
    for (int i=0;i<int(ktmp.size());i++)
    {
        KeyOfTree.push_back(ktmp[i]);
    }
    buildIndex(root, Q, 0, sp, true, true);
}

template<class kType>
void aptree<kType>::Clear(node *n)
{
    if(Qmap.size()==0) return;
    if(n->property == Q)
    {
        qNode *qn = reinterpret_cast<qNode*>(n);
        delete qn;
    }
    else if(n->property == K)
    {
        kNode *kn = reinterpret_cast<kNode*>(n);
        for(int j=0; j<int(kn->N.size()); ++j)
        {
            Clear(kn->N[j]);
        }
        delete kn;
    }
    else
    {
        sNode *sn = reinterpret_cast<sNode*>(n);
        for(int j=0; j<int(sn->N.size()); ++j)
        {
            Clear(sn->N[j]);
        }
        delete sn;
    }
}

template<class kType>
double aptree<kType>::calCost(const vector<vector<kType> > &re,
        const vector<typename map<int, query>::iterator> &Qset,int total, int l)const
{
    vector<double> cutP;    //每一个分支落上的概率
    for(int i=0; i<int(re.size()); ++i)
    {
        cutP.push_back(0);
        for(int j=0; j<int(re[i].size()); ++j)
        {
            cutP[i] += frequency.find(re[i][j])->second/double(total);
        }
    }
    double cost = 0;
    int index;
    for(int i=0; i<int(Qset.size()); ++i)
    {
        index = cFind(re, Qset[i]->second.key[l], 0, re.size()-1);
        if(index!=-1)
        {
            cost+=cutP[index];
        }
    }
    return cost;
}
#endif
