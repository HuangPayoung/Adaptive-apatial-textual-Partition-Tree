#ifdef BUILDINDEX

template<class kType>
class aptree;

template<class kType>
void aptree<kType>::buildIndex(const vector<query> &Q)
{
    vector<typename map<int, query>::iterator> Qset;
    for(int i=0; i<int(Q.size()); ++i)
    {
        Qmap.insert(pair<int,query>(currentId, Q[i]));
        Qset.push_back(Qmap.find(i));
        ++currentId;    //Qmap中id与Q中下标一一对应
    }
    double x1min=Q[0].reg.x1, x2max=Q[0].reg.x2, y1min=Q[0].reg.y1, y2max=Q[0].reg.y2;
    for(int i=0; i<int(Q.size()); ++i)
    {
        if(Q[i].reg.x1<x1min) x1min = Q[i].reg.x1;
        if(Q[i].reg.x2>x2max) x2max = Q[i].reg.x2;
        if(Q[i].reg.y1<y1min) y1min = Q[i].reg.y1;
        if(Q[i].reg.y2>y2max) y2max = Q[i].reg.y2;
    }
    space sp = {x1min, x2max, y1min, y2max};
    SpaceOfTree = sp;
    vector<kType> ktmp;
    typename map<kType, int>::iterator keyPos;
    for (int i=0;i<int(Qset.size());i++)
    {
        for (int j=0;j<int(Qset[i]->second.key.size());j++)
        {
            keyPos = frequency.find(Qset[i]->second.key[j]);
            if(keyPos == frequency.end())//没找到
            {
                frequency.insert(pair<kType, int>(Qset[i]->second.key[j],0 ));
                ktmp.push_back(Qset[i]->second.key[j]);
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
    buildIndex(root, Qset, 0, sp, true, true);
}

template<class kType>
void aptree<kType>::buildIndex
(node *&N, const vector<typename map<int, query>::iterator> Qset,
  int l, const space &s, bool kP, bool sP)
{
    if((!kP&&!sP)||int(Qset.size())<thre)
    {
        N = new qNode;
        qNode *qN=reinterpret_cast<qNode*>(N);
        qN->property = Q;
        for(int i=0; i<int(Qset.size()); ++i)
        {
            qN->q.push_back(Qset[i]);
            Qset[i]->second.N.push_back(qN);
        }
        return;
    }
    double Ck=Inf, Cs=Inf;
    vector<vector<kType> > reK;
    vector<space> reS;
    vector<kType> Kscope, kScopeT;
    if(kP)
    {
        for(int i=0; i<int(Qset.size()); ++i)
        {
            if(int(Qset[i]->second.key.size()) > l)
                kScopeT.push_back(Qset[i]->second.key[l]);
        }
        sort(kScopeT.begin(),kScopeT.end());
        if(kScopeT.size()>1)
        {
            Kscope.push_back(kScopeT[0]);
            for(int j=1; j<int(kScopeT.size()); ++j)
            {
                if(kScopeT[j]!=kScopeT[j-1]) Kscope.push_back(kScopeT[j]);
            }
        }
        else Kscope = kScopeT;
        Ck = Kpart(Qset, l, Kscope, reK);
    }
    if(sP)
    {
        Cs = Spart(Qset, s, reS);
    }
    if(Ck<Cs)
    {
        N = new kNode;
        kNode *kN = reinterpret_cast<kNode*>(N);
        kN->property = K;
        kN->offset = l;
        kN->cut = reK;

        for(int i=0; i<int(reK.size()); ++i){kN->N.push_back(NULL);}

        vector<typename map<int, query>::iterator> QsetT;
        bool flag=false;
        for(int i=0; i<int(Qset.size()); ++i)   //检查dummy cut
        {
            if(int(Qset[i]->second.key.size()) <= l)
            {QsetT.push_back(Qset[i]); flag=true;}
        }
        if(flag)
        {
            kN->N.push_back(NULL);
            buildIndex(kN->N.back(), QsetT, l+1, s, false, sP);
        }

        vector<typename map<int, query>::iterator> *QsetTset;
        QsetTset = new vector<typename map<int, query>::iterator>[reK.size()];
        int cIndex;
        if(flag)
        {
            for(int i=0; i<int(Qset.size()); ++i)
            {
                if(int(Qset[i]->second.key.size()) > l)
                {
                    cIndex = cFind(reK, Qset[i]->second.key[l], 0, reK.size()-1);
                    QsetTset[cIndex].push_back(Qset[i]);
                }
            }
        }
        else
        {
            for(int i=0; i<int(Qset.size()); ++i)
            {
                cIndex = cFind(reK, Qset[i]->second.key[l],  0, reK.size()-1);
                QsetTset[cIndex].push_back(Qset[i]);
            }
        }
        for(int i=0; i<int(reK.size()); ++i)
            buildIndex(kN->N[i], QsetTset[i], l+1, s, kP, sP);
        delete[] QsetTset;
    }
    else//s-node
    {
        N = new sNode;
        sNode *sN = reinterpret_cast<sNode*>(N);
        sN->property = S;
        sN->cell = reS;
        sN->s = s;
        for(int i=0; i<int(reS.size()); ++i){sN->N.push_back(NULL);}

        vector<typename map<int, query>::iterator> QsetT;
        bool flag=false;
        for(int i=0; i<int(Qset.size()); ++i)   //检查dummy cut
        {
            if(cover(sN->s,Qset[i]->second.reg))
            {QsetT.push_back(Qset[i]); flag=true;}
        }
        if(flag)
        {
            sN->N.push_back(NULL);
            buildIndex(sN->N.back(), QsetT, l, s, kP, false);
        }

        vector<typename map<int, query>::iterator> *QsetTset;
        QsetTset = new vector<typename map<int, query>::iterator>[reS.size()];
        for(int i=0; i<int(reS.size()); ++i)
        {
            if(flag)
            {
                for(int j=0; j<int(Qset.size()); ++j)
                {
                    if(overlap(Qset[j]->second.reg, reS[i])&&
                        !cover(sN->s,Qset[j]->second.reg))
                        QsetTset[i].push_back(Qset[j]);
                }
            }
            else
            {
                for(int j=0; j<int(Qset.size()); ++j)
                {
                    if(overlap(Qset[j]->second.reg, reS[i]))
                    {
                        QsetTset[i].push_back(Qset[j]);
                    }
                }
            }
        }
        for(int i=0; i<int(reS.size()); ++i)
            buildIndex(sN->N[i], QsetTset[i], l, reS[i], kP, sP);
        delete[] QsetTset;
    }
}

#endif
