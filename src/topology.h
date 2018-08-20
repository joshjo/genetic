#ifndef TOPOLOGY_H
#define TOPOLOGY_H
using namespace std;

#include "Graph/graph.h"
#include "chromosome.h"
// #include "Genocop.h"


class TopologyTraits{
public:
    typedef Chromosome ch;
    // typedef ChromoGeno<double> ch;
    typedef ch * N;
    typedef int E;
};

template <class Tr>
class Topology{
public:
    typedef typename Tr::N N;
    typedef typename Tr::ch ch;
    typedef typename Tr::E E;

    typedef Graph<Tr> GraphT;
    typedef Node<GraphT> node;
    typedef Edge<GraphT> edge;
    typedef vector<node*> NodeSeq;
    typedef list<edge*> EdgeSeq;

    typedef typename EdgeSeq::iterator EdgeIte;
    typedef typename NodeSeq::iterator NodeIte;

private:
    bool old_build_edges(bool toroidal) {
        int matrix_size = matrix.size();
        int weight = 1;
        bool dir = false;
        for (int i=0; i < matrix_size; ++i) {
            int matrix_size_i = matrix[i].size();
            for (int j=0; j < matrix_size_i; j++) {
                int p_i = i+1;
                if(p_i < matrix_size){
                    G->insertEdge(weight, matrix[i][j], matrix[p_i][j], dir);
                } else if(toroidal){
                    G->insertEdge(weight, matrix[i][j], matrix[0][j], dir);
                }
                p_i = i-1;
                if(p_i >= 0){
                    G->insertEdge(weight, matrix[i][j], matrix[p_i][j], dir);
                } else if(toroidal){
                    G->insertEdge(weight, matrix[i][j], matrix[matrix_size_i-1][j], dir);
                }
                int p_j = j+1;
                if(p_j < matrix_size_i){
                    G->insertEdge(weight, matrix[i][j], matrix[i][p_j], dir);
                } else if(toroidal){
                    G->insertEdge(weight, matrix[i][j], matrix[i][0], dir);
                }
                p_j = j-1;
                if(p_j >= 0){
                    G->insertEdge(weight, matrix[i][j], matrix[i][p_j], dir);
                } else if(toroidal){
                    G->insertEdge(weight, matrix[i][j], matrix[i][matrix_size_i-1], dir);
                }
            }
        }
        return true;
    }

    bool build_edges(bool toroidal, int radius, char mode='*'){
        /**
         * Init function. Used in the constructor.
         * @param toroidal: Boolean to decide if the matrix is toroidal
         * @param radius: Radius scope of each node
         * @param mode: How to connect nodes: + in cross mode, * in asterisk mode.
         * **/
        int matrix_size = matrix.size();
        bool dir = false;
        int scaling_factor = 2;
        for (int f=0; f<matrix_size; ++f){
            int matrix_size_i = matrix[f].size();
            for (int c=0; c<matrix_size_i; c++){
                for (int i=f-radius; i <= f + radius; ++i){
                    for (int j=c-radius; j <= c + radius; ++j){
                        if(i == f && j == c) continue;
                        if(mode == '+' && !(i==f || j==c)) continue;
                        if(mode == '*' && (abs(i-f) > radius-abs(j-c))) continue;
                        int w = sqrt(pow(abs(i) - f, 2) + pow(abs(j) - c, 2));
                        int p_i = i;
                        int p_j = j;
                        bool modified = false;
                        if(p_i < 0){
                            p_i = matrix_size+i;
                            modified = true;
                        } else if(p_i >= matrix_size){
                            p_i = p_i % matrix_size;
                            modified = true;
                        }
                        if(p_j < 0){
                            p_j = matrix_size_i+j;
                            modified = true;
                        } else if(p_j >= matrix_size_i){
                            p_j = p_j % matrix_size_i;
                            modified = true;
                        }
                        if(toroidal || !modified){
                            G->insertEdge(
                                w, matrix[f][c], matrix[p_i][p_j], dir);
                        }
                    }
                }
            }
        }
        return true;
    }

    bool init(  int dimension, NumberParameters & number_params,
                bool toroidal, int c, int r, int radius, char mode){
        /**
         * Init function. Used in the constructor.
         * @param dimension:
         * @param number_params: A class
         * @param toroidal: Boolean to decide if the matrix is toroidal
         * @param c: Number of columns
         * @param r: Number of rows
         * @param radius: Radius scope of each node
         * @param mode: How to connect nodes: + in cross mode, * in asterisk mode.
         * **/
        G = new GraphT();
        int counter = 0;
        for(int i=0; i<c; i++){
            vector <node *> vnodes;
            for(int j = 0; j < r; j++){
                N tmp = new ch(dimension, number_params);
                node* n = G->insertNode(tmp, i, j);
                vnodes.push_back(n);
                counter++;
            }
            matrix.push_back(vnodes);
        }
        return build_edges(toroidal, radius, mode);
    }

public:
    GraphT * G;
    vector <vector <node *> > matrix;
    int cols;
    int rows;

    Topology(int dimension, NumberParameters & number_params,
            bool toroidal=true, int c=10, int r=10, int ratio=1,
            char mode='+'){
        cols = c;
        rows = r;
        init(dimension, number_params, toroidal, c, r, ratio, mode);
    };

    void print_status(int verbosity=1){
        for(int i=0; i<matrix.size(); ++i){
            for(int j=0;j<matrix[i].size(); ++j){
                cout << i << " - " << j << endl;
                matrix[i][j]->Data->print_value();
                if(verbosity > 1){
                    cout << "Edges: " << endl;
                    EdgeSeq es = matrix[i][j]->edges;
                    for(EdgeIte it=es.begin(); it!=es.end(); ++it){
                        node* n = (*it)->nodes[1];
                        cout << "x: " << n->x << ", y: " << n->y << endl;
                        n->Data->print_value();
                    }
                }
                cout << endl;
            }
        }
    }
};

#endif // TOPOLOGY_H
