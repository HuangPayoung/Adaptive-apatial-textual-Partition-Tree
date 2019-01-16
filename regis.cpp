#ifdef REGIS

template <class kType>
int aptree<kType>::regis(const query &q)
{
    Qmap.insert(pair<int,query>(currentId, q));
    currentId++;
    bool flag=true;
    for (int i=0;i<int(q.key.size());i++)
    {
        if (kFind(KeyOfTree,q.key[i],0,KeyOfTree.size()-1)<0) {flag=false;break;}
    }
    if (!flag||!overlap(q.reg, SpaceOfTree)) {reconstruct();/*throw 1;*/}
    regis(Qmap.find(currentId-1)->second, root);
    if(++timesChange>KL){reconstruct(); timesChange=0;}
    return currentId-1;
}

template <class kType>
void aptree<kType>::regis(query &q,node *N)
{
    if (N->property==Q)
    {
        qNode *qN=reinterpret_cast<qNode*>(N);
        typename map<int,query>::iterator qtmp=Qmap.find(currentId-1);
        qN->q.push_back(qtmp);
        q.N.push_back(qN);
    }
    else if(N->property==K)
    {
        kNode *kN=reinterpret_cast<kNode*>(N);
        if (kN->offset>int(q.key.size()))              //dummy cut
        {
            if (kN->N.size()>kN->cut.size())      //存在dummy cut
            {
                regis(q, kN->N.back());
            }
            else                                  //不存在dummy cut
            {
                node *qN=new qNode;
                kN->N.push_back(qN);
                qN->property=Q;
                regis(q,qN);
            }
        }
        else
        {
            int index = cFind(kN->cut, q.key[kN->offset], 0, kN->cut.size()-1);
            if(index==-1)
            {
                vector<kType> cutT;
                cutT.push_back(q.key[kN->offset]);
                kN->cut.push_back(cutT);
                node *newC = new qNode;
                newC->property = Q;
                kN->N.push_back(newC);
                regis(q, newC);
            }
            else regis(q, kN->N[index]);
        }
    }
    else
    {
        sNode *sN=reinterpret_cast<sNode*>(N);
        if (cover(sN->s,q.reg))                    //dummy cell
        {
            if (sN->N.size()>sN->cell.size())      //存在dummy cell
            {
                regis(q,reinterpret_cast<node*>(sN->N.back()));
            }
            else                                   //不存在dummy cell
            {
                node *qN=new qNode;
                sN->N.push_back(qN);
                qN->property=Q;
                regis(q,qN);
            }
        }
        else
        {
            for (int i=0;i<int(sN->cell.size());i++)
            {
                if (overlap(sN->cell[i],q.reg))
                {
                    regis(q,sN->N[i]);
                }
            }
        }
    }

}

template <class kType>
void aptree<kType>::deregis(int ID)
{
    typename map<int,query>::iterator qtmp=Qmap.find(ID);
    if(qtmp==Qmap.end()) cerr<<endl<<"找不到ID"<<ID<<endl;
    query q=qtmp->second;
    qNode *qN;
    for(int i=0; i<int(q.N.size()); ++i)
    {
        qN = reinterpret_cast<qNode*>(q.N[i]);
        typename vector<typename map<int,query>::iterator>::iterator pos;
        bool flag=true;//检查错误
        for(pos=qN->q.begin(); pos!=qN->q.end(); ++pos)
        {
            if(*pos==qtmp)
            {
                qN->q.erase(pos);
                flag = false;
                break;
            }
        }
        if(flag)
        {
            cerr<<endl<<"从ID"<<ID<<"的第"<<i<<"个结点未找到应有的query"<<endl;
        }
    }
    Qmap.erase(ID);
    if(++timesChange>KL){reconstruct(); timesChange=0;}
    return;
}

#endif
