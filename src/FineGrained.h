#ifndef FINEGRAINED_H
#define FINEGRAINED_H
#include "Topology.h"

template <class Tr, class c_cmp>
class FineGrained{
private:
public:
    typedef typename Tr::N N;
    typedef typename Tr::ch ch;
    typedef typename Tr::E E;
    Topology<Tr>* T;
    typedef Graph<Tr> GraphT;
    typedef Node<GraphT> node;
    typedef Edge<GraphT> edge;
    typedef vector<node*> NodeSeq;
    typedef list<edge*> EdgeSeq;

public:
    double prob_crossing;
    double prob_mutation;
    opt_function eval_function;
    function <bool(double, double)> cmp;
    map<N, double> mfitness;
    map<N, node*> node_chromosome;
    map<node*, bool> worsts;
    map<node*, N> temp_gen;
    typedef typename EdgeSeq::iterator EdgeIte;
    typedef typename NodeSeq::iterator NodeIte;
    vector<double> bests;
    vector<double> averages;
    vector<double> offline;
    vector<double> online;
    N best_elem;
    double len_worsts_inc;
    int num_nodes;
    bool has_converge;
    bool optimization;
    double optimal;
    double t_mutation;
    int iteration_convergence;

    FineGrained(int dimension, opt_function fun,
            NumberParameters & number_params, double o, bool toroidal=false,
            int c=10, int r=10, int ratio=1, char mode='+'){
        cmp = c_cmp();
        if(cmp(1, 2)){
            optimization = MINIMIZATION;
        } else {
            optimization = MAXIMIZATION;
        }
        num_nodes = c * r;
        T = new Topology<Tr>(
            dimension, number_params, toroidal, c, r, ratio, mode);
        eval_function = fun;
        prob_crossing = .75;
        prob_mutation = .1;
        len_worsts_inc = .99;
        has_converge = false;
        optimal = o;
        iteration_convergence = -1;
        t_mutation = 1;
        thread_count = 4;
    }

    // OpenMP
    //

    void evolution(
            int num_iterations=100,
            bool de=false
        ){

        // printf("===> %d", thread_count);

        // # pragma omp parallel for num_threads(thread_count)
        for(int k = 0; k < num_iterations; ++k){
            precalculate_fitnesses();
            NodeSeq nodes = T->G->nodes;
            if (de) {
                // for(auto it: nodes) node_evolution_de(it);
            } else {
                for(auto it: nodes) node_evolution(it);
            }
            for (auto it: nodes) {
                if(temp_gen[it]){
                    it->Data = temp_gen[it];
                }
            }
            t_mutation *= 0.99;
            // cout << "Iteration: " << k+1 << endl;
            status_iteration();
            // cout << endl;
        }
    }

    void parallel_evolution(int num_iterations, bool de=false){
        NodeSeq nodes = T->G->nodes;
        # pragma omp parallel num_threads(num_nodes)
        for (int k = 0; k < num_iterations; ++k) {
            #pragma omp sections
            {
                precalculate_fitnesses();
            }
            # pragma omp for
            for(int i=0; i<nodes.size(); ++i){
                if(de) node_evolution_de_thread();
                else node_evolution_thread();
            }
            #pragma omp sections
            {

                #pragma omp section
                for(auto it: nodes){
                    if (temp_gen[it]) {
                        it->Data = temp_gen[it];
                    }
                }
                #pragma omp section
                status_iteration();
            }
        }
    }

    void node_evolution_thread(){
        int my_rank = omp_get_thread_num();
        node* cnode = T->G->nodes[my_rank];
        N current = cnode->Data;
        N copy_current = new ch(*current);
        double rand_crossing = reindom();
        double rand_mutation = reindom();
        bool is_save = false;
        if(rand_crossing < this->prob_crossing){
            N other = roulette_selection(cnode);
            N copy_other = new ch(*other);
            copy_current->crossing(*copy_other, optimization);
            is_save = true;
        }
        if(rand_mutation < this->prob_mutation){
            // Modified:
            // Was copy_current->mutation(t_mutation);
            copy_current->mutation();
            is_save = true;
        }
        double fitness_current = 0;
        #pragma omp critical
        fitness_current = mfitness[current];
        if(cmp(eval_function(copy_current->get_real()), fitness_current)){
            #pragma omp critical
            temp_gen[cnode] = copy_current;
        }
        if(is_save) temp_gen[cnode] = copy_current;
    }

    void node_evolution_de_thread(){
        int my_rank = omp_get_thread_num();
        node* cnode = T->G->nodes[my_rank];
        N current = cnode->Data;
        N copy_current = new ch(*current);
        vector<N> rnodes = random_nodes(cnode);
        N v = rnodes[0]->mutation_de(*rnodes[1], *rnodes[2]);
        copy_current->crossing_de(*v, 0.85);
        double fitness_current = 0;
        #pragma omp critical
        fitness_current = mfitness[current];
        if(cmp(eval_function(copy_current->get_real()), fitness_current)){
            #pragma omp critical
            temp_gen[cnode] = copy_current;
        }
    }

    void node_evolution(node* cnode){
        if(!worsts[cnode]) {
            return;
        }
        N current = cnode->Data;
        N copy_current = new ch(*current);
        bool will_update = false;
        if(reindom() < this->prob_crossing){
            N other = tournament_selection(cnode);
            N copy_other = new ch(*other);
            // Modified: old was:
            // copy_current->crossing(*copy_other, optimization);
            copy_current->crossing(*copy_other);
            will_update = true;
        }
        if(reindom() < this->prob_mutation){
            // Modified:
            // Was copy_current->mutation(t_mutation);
            copy_current->mutation();
            will_update = true;
        }
        // if(cmp(eval_function(copy_current->get_real()), mfitness[current])){
        //     temp_gen[cnode] = copy_current;
        // }
        if(will_update) temp_gen[cnode] = copy_current;
    }

    void node_evolution_de(node* cnode){
        N current = cnode->Data;
        vector<N> rnodes = random_nodes(cnode);
        N copy_current = new ch(*current);
        N v = rnodes[0]->mutation_de(*rnodes[1], *rnodes[2]);
        copy_current->crossing_de(*v, 0.85);
        if (cmp(eval_function(copy_current->get_real()), mfitness[current])) {
            temp_gen[cnode] = copy_current;
        }
    }

    void plot_bests(ostream & os){
        int size = bests.size();
        for(int i=0; i<size; ++i){
            os << i << " " << bests[i] << " " << averages[i] << endl;
        }
    }

    vector<vector<double> > population_points(){
        vector<vector<double> > resp;
        for(auto it: T->G->nodes){
            resp.push_back(it->Data->get_real());
        }
        return resp;
    }

    void precalculate_fitnesses(){
        temp_gen.clear();
        mfitness.clear();
        worsts.clear();
        multimap<double, node*, c_cmp> tmp_worst;
        for(auto it: T->G->nodes){
            N data = it->Data;
            double value = eval_function(data->get_real());
            mfitness[data] = value;
            tmp_worst.insert(pair<double, node*>(value, it));
            node_chromosome[data] = it;
        }
        int count = 0;
        // int len_worsts = (T->cols*T->rows*0.3);
        int len_worsts = int(T->cols * T->rows * len_worsts_inc);
        for(auto it = tmp_worst.rbegin(); it != tmp_worst.rend(); it++){
            worsts[it->second] = true;
            if(count > len_worsts) break;
            count++;
        }
        len_worsts_inc *= 0.99;
    }

    vector<N> random_nodes(node* cnode, int size_N=3){
        EdgeSeq es = cnode->edges;
        int size_e = es.size();
        vector<N> vresult;
        if(size_e < size_N){
            cout << "[ERROR] Neighborhood too small" << endl;
            return vresult;
        }
        // map<double, N, c_cmp> result;
        double n_val = float(size_N)/es.size();
        map<N, bool> n_selected;
        int count = 0;
        // for(int i=0; i<size_N;){
        while(n_selected.size() < size_N){
            EdgeIte it = es.begin();
            advance(it, reindint(0, size_e-1));
            N n = (*it)->nodes[1]->Data;
            double prob_chosen = 1.0/(*it)->Data;
            if(reindom()<prob_chosen){
                n_selected[n];
            }
            // cout <<  prob <<" " << reindom()<< endl;
            // if(reindom()<prob)
        }
        for(auto it: n_selected)
            vresult.push_back(it.first);
        // for(auto it=result.rbegin(); it != result.rend(); it++)
        //     vresult.push_back(it->second);
        return vresult;
    }

    N roulette_selection(node* cnode){
        N current = cnode->Data;
        EdgeSeq es = cnode->edges;
        double fs = 0;
        double minimun = mfitness[es.front()->nodes[1]->Data];
        double maximum = minimun;
        for(EdgeIte it=es.begin(); it!=es.end(); it++){
            double value = mfitness[(*it)->nodes[1]->Data];
            fs += value;
            if(maximum < value){
                maximum = value;
            }
            if(minimun > value){
                minimun = value;
            }
        }
        double p = reindom(0, fs);
        N other;
        double sum_min_max = minimun + maximum;

        if(cmp(1, 0)){
            for(auto it: cnode->edges){
                other = it->nodes[1]->Data;
                p -= mfitness[other];
                if(p <= 0) break;
            }
        } else{
            for(auto it: cnode->edges){
                other = it->nodes[1]->Data;
                p -= sum_min_max - mfitness[other];
                if(p <= 0) break;
            }
        }
        return other;
    }

    N tournament_selection(node* cnode){
        N best = cnode->edges.front()->nodes[1]->Data;
        for(auto it: cnode->edges){
            N other = it->nodes[1]->Data;
            if(cmp(mfitness[other], mfitness[best])) best = other;
        }
        return best;
    }

    void status_iteration(){
        int t = bests.size();
        double avg = 0;
        N best = mfitness.begin()->first;
        for(auto it: mfitness){
            avg += it.second;
            if(cmp(mfitness[it.first], mfitness[best])){
                best = it.first;
            }
        }
        if(mfitness[best] == optimal && !has_converge){
            has_converge = true;
            iteration_convergence = bests.size();
        }
        node* n_best = node_chromosome[best];
        avg = avg/mfitness.size();
        bests.push_back(mfitness[best]);
        averages.push_back(avg);
        double sum_avgs = 0;
        double sum_best = 0;
        for (double n : averages) sum_avgs += n;
        for (double n : bests) sum_best += n;
        double delta = float(1)/(t+1);
        online.push_back(delta*sum_avgs);
        offline.push_back(delta*sum_best);
        // cout << "Best: " << mfitness[best]  << " (";
        // cout << n_best->x << ", " << n_best->y << ")" << endl;
        // best->print_value();
        // cout << avg << endl;
        // cout << endl;
        // cout << "-> Avg: " << avg;
        // cout << endl << endl;

    }
};

#endif // FINEGRAINED_H
