#ifndef Edge_H
#define Edge_H

#include "Node.h"

template <class G>
class Edge
{
    public:
        typedef typename G::E E;
        typedef typename G::node node;

        node* nodes[2];

    public:
        E Data;
        bool dir;

    public:
        Edge(E d,node* i,node* f,bool di):Data(d),dir(di){
            nodes[0]=i;
            nodes[1]=f;
        }

        void kill(){
            (nodes[0]->edges).remove(this); nodes[0]=NULL; delete nodes[0];
            (nodes[1]->edges).remove(this); nodes[1]=NULL; delete nodes[1];
            delete this;
        }

        E getData(){return Data;}

        bool& getDir(){return dir;}
};

#endif
