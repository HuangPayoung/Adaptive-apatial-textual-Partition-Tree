#ifdef PRINT

template<class kType>
void aptree<kType>::print(const char name[])const
{
    cout<<"参数设定:"<<endl;
    cout<<"thre:"<<thre<<' '<<"KL:"<<KL<<' '<<"f:"<<kf<<' '<<"m:"<<sfm<<' '<<"n:"<<sfn<<endl;
    queue<pair<int, node*> > que;
    que.push(pair<int, node*>(0,root));
    pair<int,node*> tmp;
    ofstream fout(name);
    int last=0;
    int Count=0;
    cout<<"printing file aptree "<<name<<endl;
    while(!que.empty())
    {
        tmp = que.front();
        que.pop();
        if(tmp.first>last) Count = 0;
        ++Count;
        last = tmp.first;
        fout<<tmp.first<<' '<<Count<<endl;
        if(tmp.second->property==K)
        {
            kNode *kn=reinterpret_cast<kNode*>(tmp.second);
            fout<<"kNode\n";
            fout<<"offset:"<<kn->offset<<endl;
            if(kn->cut.size()<kn->N.size())fout<<"dummy"; else fout<<"ordinary";
            fout<<endl;
            for(int i=0; i<int(kn->cut.size()); ++i)
            {
                for(int j=0; j<int(kn->cut[i].size()); ++j)
                    fout<<kn->cut[i][j]<<'\t';
                fout<<endl;
            }
            for(int i=0; i<int(kn->N.size()); ++i)
            {
                que.push(pair<int, node*>(tmp.first+1, kn->N[i]));
            }
            fout<<endl;
        }
        if(tmp.second->property==S)
        {
            fout<<"sNode\n";
            sNode *sn=reinterpret_cast<sNode*>(tmp.second);
            fout<<sn->s<<endl;
            if(sn->cell.size()<sn->N.size())fout<<"dummy"; else fout<<"ordinary";
            fout<<endl;
            for(int i=0; i<int(sn->cell.size()); ++i)
            {
                fout<<sn->cell[i]<<endl;
            }
            for(int i=0; i<int(sn->N.size()); ++i)
            {
                que.push(pair<int, node*>(tmp.first+1, sn->N[i]));
            }
            fout<<endl;
        }
        if(tmp.second->property==Q)
        {
            fout<<"qNode\n";
            qNode *qn=reinterpret_cast<qNode*>(tmp.second);
            for(int i=0; i<int(qn->q.size()); ++i)
                fout << qn->q[i]->first << ' ';
            fout<<endl<<endl;
        }
    }
    fout.close();
    cout<<"printed"<<endl<<endl;
}

#endif // PRINT
