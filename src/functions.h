#ifndef FUCNTIONS_H
#define FUCNTIONS_H
#include <random>

double reindom(double from=0, double to=1){
    if(to > from){
        swap(to, from);
    }
    uniform_real_distribution<double> unif(from, to);
    random_device rand_dev;
    mt19937 rand_engine(rand_dev());
    return unif(rand_engine);
}

double reindint(int from, int to){
    // Redone with random_device as reindom
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    if(to < from){
        swap(to, from);
    }
    mt19937 r(seed);
    return (r()%(to-from+1) + from);
}

double reindint(int to){
    return reindint(0, to);
}


int get_l(double left, double right, int precission){
    return ceil(log2(right - left) + precission * log2(10));
}


template<typename T>
void print_vector(vector <T> array, bool print_endl=false){
    cout << "vector [";
    for (auto it = array.begin(); it != array.end(); ++it)
        cout << *it << ", ";
    cout << "] ";
    if(print_endl) cout << endl;
}


template<typename T>
void swap_vector(vector <T> & array, vector <T> & other, int ini, int end){
    if(end >= array.size())
        end = array.size()-1;
    if(ini > array.size() || end > array.size() || array.size() != other.size())
        return;
    for(int i=ini; i < end; ++i)
        swap(array[i], other[i]);
}


vector <int> get_random_points(int number_of_points, int max_limit){
    vector <int> resp;
    int r = reindint(max_limit-1);
    for(int i=0; i<number_of_points; ++i)
        resp.push_back(r);
    return resp;
}

#endif  // FUCNTIONS_H
