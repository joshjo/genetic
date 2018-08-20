#ifndef GENETIC_H
#define GENETIC_H

#include <omp.h>


template <class C>
class Genetic{
public:
    typedef vector <C *> chromosome_vector;

    opt_function eval_function;
    chromosome_vector population;
    function <bool(double, double)> cmp;
    int len_worsts;
    string type_opt_str;
    vector<double> bests;
    vector<double> averages;
    vector<double> offline;
    vector<double> online;
    double optimal;
    double len_worsts_inc;
    int iteration_convergence;
    bool has_converge;
    bool optimization;
    bool type_opt;
    double mode_worst;
    map<C*, double> pfitnesses;



    Genetic(
        int population_length,
        int dimension,
        opt_function fun,
        NumberParameters & number_params,
        bool type_opt,
        double o=0
    ){
        for(int i=0; i<population_length; ++i)
            population.push_back(new C(dimension, number_params));
        optimal = o;
        optimization = type_opt;
        cmp = less<double>();
        has_converge = false;
        type_opt_str = "Min: ";
        if(type_opt == MAXIMIZATION){
            cmp = greater<double>();
            type_opt_str = "Max: ";
        }
        iteration_convergence = -1;
        eval_function = fun;
        mode_worst = 0;
        len_worsts_inc = 1;
    }

    struct doCompare{
        doCompare( const Genetic& info ) : m_info(info) { }
        const Genetic& m_info;
        bool operator()(C * c1, C * c2){
            return m_info.cmp(m_info.eval_function(c2->get_real()), m_info.eval_function(c1->get_real()));
        }
    };

    chromosome_vector tournament_selection(){
        chromosome_vector new_population;
        int size = population.size();
        if(mode_worst <= 0){
            len_worsts_inc *= 0.98;
        } else{
            len_worsts_inc = mode_worst;
        }
        len_worsts = int(len_worsts_inc*size);
        for(int i=0; i<len_worsts; ++i){
            C * first_item = population[reindint(size-1)];
            C * second_item = population[reindint(size-1)];
            if(cmp(pfitnesses[first_item], pfitnesses[second_item]) ) {
                new_population.push_back(new C(*first_item));
            } else{
                new_population.push_back(new C(*second_item));
            }
        }
        return new_population;
    }

    void replace_worsts(chromosome_vector& new_population){
        int i;
        sort(population.begin(), population.end(), doCompare(*this));

        # pragma omp parallel for num_threads(4) private(i) \
                shared(population, new_population)
        for (i = 0; i < len_worsts; i += 1){
            pfitnesses.erase(population[i]);
            population[i] = new_population[i];
            pfitnesses[new_population[i]] = eval_function(new_population[i]->get_real());
        }

    }

    void precalculate_fitnesses(){
        int i, size;
        size = population.size();
        # pragma omp parallel for num_threads(4) private(i) \
                shared(new_population)
        for (i = 0; i < size; i += 1) {
            pfitnesses[population[i]] = eval_function(population[i]->get_real());
        }
    }

    void status_iteration(){
        int t, i, population_size;
        t = bests.size();
        population_size = population.size();
        C * best = population[0];
        double avg = 0;
        for (auto i: pfitnesses) {
            cout << i.first << " - " << i.second << endl;
        }
        for(i = 0; i < population_size; i += 1) {
            cout << "-" << population[i] << endl;
        }
        for(auto i:population){
            // cout << "i: " << i << endl;
            // cout << pfitnesses[i] << " - " << i->get_real() << endl;
            double value = eval_function(i->get_real());
            avg += value;
            if(cmp(value, eval_function(best->get_real()))){
                best = i;
            }
        }
        double besty = eval_function(best->get_real());
        if(besty == optimal && !has_converge){
            has_converge = true;
            iteration_convergence = averages.size();
        }
        avg = avg/population.size();
        averages.push_back(avg);
        bests.push_back(besty);
        double sum_avgs = 0;
        double sum_best = 0;
        for (double n : averages) sum_avgs += n;
        for (double n : bests) sum_best += n;
        double delta = float(1)/(t+1);
        online.push_back(delta*sum_avgs);
        offline.push_back(delta*sum_best);
        // cout << type_opt_str << eval_function(best->get_real()) << " -> ";
        // best->print_value();
        // cout << endl <<endl;
    }

    void plot_bests(ostream & os){
        int size = bests.size();
        for(int i=0; i<size; ++i){
            os << i << " " << bests[i] << " " << averages[i] << endl;
        }
    }

    vector<vector<double> > population_points(){
        vector<vector<double> > resp;
        for(auto it: population){
            resp.push_back(it->get_real());
        }
        return resp;
    }

    void plot_population_points(ostream & os){
        for(auto it: population){
            for(auto data: it->get_real()){
                os << data << " ";
            }
            os << endl;
        }
    }
};


#endif // GENETIC_H
