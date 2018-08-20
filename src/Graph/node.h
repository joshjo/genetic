#ifndef Node_H
#define Node_H

template <class G>
class Node
{
    public:
        typedef typename G::N N;
        typedef typename G::E E;
        typedef typename G::edge edge;
        typedef typename G::EdgeSeq EdgeSeq;

        EdgeSeq edges;

    public:
        N Data;
        float x;
        float y;

    public:
        Node(N d, float p_x, float p_y):Data(d),x(p_x),y(p_y){};

        N getData() {return Data;}

        float getX() {return x;}

        float getY() {return y;}
};

#endif
