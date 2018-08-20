#ifndef Grafo_H
#define Grafo_H

#include "includes.h"
#include "Node.h"
#include "Edge.h"

using namespace std;


struct GraphTraits{
    typedef char * N;
    typedef int E;
};


template <class Tr>
class Graph
{
    public:
        typedef Graph<Tr> self;
        typedef Node<self> node;
        typedef Edge<self> edge;
        typedef vector<node*> NodeSeq;
        typedef list<edge*> EdgeSeq;
        typedef typename Tr::N N;
        typedef typename Tr::E E;
        typedef typename NodeSeq::iterator NodeIte;
        typedef typename EdgeSeq::iterator EdgeIte;

    public:
        NodeSeq nodes;
        NodeIte ni;
        EdgeIte ei;

    public:
        Graph(){};

        node* searchNode(N n,int& pos){
            pos=-1;
            for(ni=nodes.begin();ni!=nodes.end();++ni){
                pos++;
                if((*ni)->getData()==n)
                    return *ni;
            }
            return NULL;
        }

        bool searchEdge(node* i, node* f){
            for(ei=(i->edges).begin();ei!=(i->edges).end();++ei){
                if((*ei)->nodes[1]==f)
                    return true;
            }
            return false;
        }

        node* insertNode(N n,double x,double y){
            node* new_node = new node(n,x,y);
            nodes.push_back(new_node);
            return new_node;
        }

        bool deleteNode(N n){
            int temp;
            node* tmp=searchNode(n,temp);
            if(tmp==NULL) return false;

            while(!(tmp->edges).empty())
                ((tmp->edges).front())->kill();

            swap(nodes.back(), nodes[temp]);
            nodes.pop_back();
            return true;
        }

        bool insertEdge(E e, node* i, node* f, bool dir){
            int temp;
            if(!f || !i) return false;

            if(searchEdge(i, f)) return false;
            edge* nuevo=new edge(e, i, f, dir);
            (i->edges).push_back(nuevo);

            if(dir) return true;

            if(searchEdge(f, i)) return false;
            nuevo=new edge(e, f, i, dir);
            (f->edges).push_back(nuevo);
            return true;
        }

        void insertEdge(E e, N ini, N fin, bool dir){
            int temp;
            node* i = searchNode(ini, temp);
            if(!i) return;
            node* f = searchNode(fin, temp);
            if(!f) return;

            if(searchEdge(i, f) || searchEdge(f, i)) return;

            edge* nuevo=new edge(e, i, f, dir);

            (i->edges).push_back(nuevo);
            (f->edges).push_back(nuevo);
        }

        bool deleteEdge(N ini, N fin){
            int temp;
            EdgeIte ei2;
            node* i = searchNode(ini, temp);
            if(!i) return false;
            node* f = searchNode(fin, temp);
            if(!f) return false;
            for(ei=(i->edges).begin(); ei!=(i->edges).end(); ++ei){
                if((*ei)->nodes[0]==i && (*ei)->nodes[1]==f &&
                        !(*ei)->getDir()){
                    (*ei)->nodes[0]=f;
                    (*ei)->nodes[1]=i;
                    (*ei)->getDir()=1;
                    return true;
                }

                if((*ei)->nodes[0]==i && (*ei)->nodes[1]==f &&
                        (*ei)->getDir()){
                    (*ei)->kill();
                    return true;
                }
            }
            return false;
        }

        float Density(){
            float n_edges=0;
            for(ni=nodes.begin();ni!=nodes.end();++ni)
                n_edges+=(*ni)->edges.size();

            return float((n_edges)/((nodes.size())*(nodes.size()-1)));
        }

        void clear(){
            for(ni=nodes.begin();ni!=nodes.end();++ni)
                (*ni)->edges.clear();
            nodes.clear();
        }

        bool empty(){
            if(nodes.empty()) return true;
            return false;
        }

        ~Graph(){};
};

#endif
