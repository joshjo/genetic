#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "randombinary.h"

class Chromosome{
public:
    vector<RandomBinary*> * value;

public:
    Chromosome(Chromosome & other){
        int array_size = other.get_value()->size();
        value = new vector<RandomBinary*> (array_size);
        for(int i=0; i<array_size; ++i)
            (*value)[i] = new RandomBinary(*other.get_value()->at(i));
    }

    Chromosome(int dimension, NumberParameters & params){
        value = new vector<RandomBinary*>;
        for(int i=0; i<dimension; ++i){
            // value->push_back(new RandomBinary(params));
            value->push_back(new RandomBinary(
                (*params.lefts)[i], (*params.rights)[i], params.precission));
        }
    }

    void crossing(Chromosome & other){
        vector <int> points = get_random_points(1, length());
        sort(points.begin(), points.end());
        points.erase(unique(points.begin(), points.end()), points.end());
        int chromosome_length = length();
        auto it_this = value->begin();
        auto it_other = other.value->begin();
        auto current_this = (*it_this)->get_array();
        auto current_other = (*it_other)->get_array();
        int diff_index = 0;
        while(!points.empty()){
            if (points.size() == 1) {
                points.push_back(chromosome_length-1);
            }
            int point_a = points[0];
            int point_b = points[1];
            if(point_a >= current_this->size()){
                diff_index += current_this->size();
                point_a -= diff_index;
                point_b -= diff_index;
                it_this++;
                it_other++;
                if (it_this == value->end())
                    break;
                current_this = (*it_this)->get_array();
                current_other = (*it_other)->get_array();
            }
            else if(point_b >= current_this->size()){
                points.insert(points.begin()+1, current_this->size());
                points.insert(points.begin()+2, current_this->size());
                point_b = points[1];
            }
            swap_vector(*current_this, *current_other, point_a, point_b);
            points.erase(points.begin(), points.begin()+2);
        }
    }

    vector <RandomBinary *> * get_value(){
        return value;
    }

    vector <double> get_real(){
        vector <double> resp;
        for(auto i: *value) resp.push_back(i->get_real());
        return resp;
    }

    int length(){
        int sum = 0;
        for(auto i: *value) sum += i->get_length();
        return sum;
    }

    void mutation(){
        bit_vector * pbv = value->at(
            reindint(0, value->size()-1))->get_array();
        int atr = reindint(0, pbv->size()-1);
        (*pbv)[atr] = ((*pbv)[atr] == 0)? 1 : 0;
    }

    void print_value(int verbosity=1){
        cout << "Chromosome: " << this << endl;
        if(verbosity > 2){
            int size = value->size();
            for(int i=0; i<size; ++i)
                print_vector(*((*value)[i]->get_array()));
        }
        vector <double> d = get_real();
        print_vector(d);
    }

    virtual ~Chromosome() {

    }
};

#endif // CHROMOSOME_H
