#ifndef NUMBERPARAMETERS_H
#define NUMBERPARAMETERS_H


struct NumberParameters{
    double left;
    double right;
    int precission;
    vector<double>* lefts;
    vector<double>* rights;
    NumberParameters(double l, double r, int p):left(l), right(r), precission(p){}
    NumberParameters(vector<double>* l, vector<double>* r, int p = 0){
        lefts = l;
        rights = r;
        precission = p;
    }

};

#endif // CHROMOSOME_H
