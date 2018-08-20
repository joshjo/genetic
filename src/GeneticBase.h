#ifndef GENETIC_BASE_H
#define GENETIC_BASE_H
#include "Chromosome.h"
#include "Genetic.h"

class GeneticBase: public Genetic<Chromosome>{
public:
    GeneticBase(
        int population_size,
        int dimension,
        opt_function fun,
        NumberParameters & number_params,
        bool type_opt
    ) : Genetic<Chromosome > (
        population_size,
        dimension,
        fun,
        number_params,
        type_opt
    ) {}

    void evolution(int num_iterations=10, double prob_crossing=.8, double prob_mutation=.25){
        for(int i = 0; i < num_iterations;++i){
            precalculate_fitnesses();
            chromosome_vector selected_population = tournament_selection();
            int size, j;
            size = selected_population.size();
            chromosome_vector new_population(size, 0);

            # pragma omp parallel for num_threads(4) private(j) \
                shared(new_population)
            for (j = 0; j < size; ++j) {
                double rand1 = reindom();
                double rand2 = reindom();
                Chromosome* current_copy = new Chromosome(
                    *selected_population[j]);
                if(rand1 < prob_crossing){
                    Chromosome * other = new Chromosome(
                        *population[reindint(len_worsts - 1)]);
                    current_copy->crossing(*other);
                }
                if(rand2 < prob_mutation){
                    current_copy->mutation();
                }
                new_population[j] = current_copy;
            }
            replace_worsts(new_population);
            status_iteration();
        }
    }
};

#endif // GENETIC_BASE_H
