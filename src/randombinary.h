#ifndef RANDOM_BINARY_H
#define RANDOM_BINARY_H

#include "includes.h"

class RandomBinary{
public:
    bit_vector * array;
    int length;
    pair<double, double> * limit;

public:
    /*
        get_l comes from functions.h
    */
    RandomBinary(RandomBinary & other){
        length = other.length;
        int array_size = other.get_array()->size();
        array = new bit_vector(array_size);
        for(int i=0; i<array_size; ++i)
            (*array)[i] = other.get_array()->at(i);
        limit = new pair <double, double> (*other.get_limit());
    }

    RandomBinary(double left, double right, int precission){
        length = get_l(left, right, precission);
        limit = new pair <double, double> (left, right);
        array = 0;
        get_array();
    }

    RandomBinary(NumberParameters& params){
        length = get_l(params.left, params.right, params.precission);
        limit = new pair <double, double> (params.left, params.right);
        array = 0;
        get_array();
    }

    bit_vector * get_array(){
        if(!array){
            array = new bit_vector();
            for (int i=0; i<length; i++)
                array->push_back(rand()%2);
        }
        return array;
    }

    pair <double, double> * get_limit(){
        return limit;
    }

    int get_length(){
        return length;
    }

    double get_real(){
        int sum = 0;
        int left = limit->first;
        int right = limit->second;
        if (!array) get_array();
        for(int i=0; i<length; ++i)
            sum += (*array)[i]*pow(2,length-i-1);
        return left+sum*(double(right-left)/(pow(2, length)-1));
    }

    virtual ~RandomBinary() {}
};

#endif // RANDOM_BINARY_H
