#include <stdio.h>
#include <stdlib.h>
#include "topology.h"
#include "testfunction.h"
// #include "FineGrained.h"
#include "GeneticBase.h"

typedef double dtype;

int thread_count = 4;


class Nice {
public:
    vector<int> * array;
    int repeat;
    Nice(int repeat=10) {
        this->repeat = repeat;
        array = new vector<int>(repeat, 0);
    }

    void parallel () {
        // for (int j = 0; j < 3; j += 1) {
            # pragma omp parallel for num_threads(4)
            for (int i = 0; i < repeat; ++i) {
                printf("Hellow %d\n", i);
            }
            printf("< === === === === === === >\n");
        // }
    }
};



int main(){
    // Nice N;
    // N.parallel();

    ofstream file;
    vector<double> lefts = {-5.12, -5.12};
    vector<double> rights = {5.12, 5.12};
    NumberParameters np(&lefts, &rights, 5);
    GeneticBase G(200, 2, drop_wave, np, MINIMIZATION);
    int iterations = 200;

    double timer = omp_get_wtime() - timer;
    G.evolution(iterations, 1, 0.35);
    timer = omp_get_wtime() - timer;

    printf("With %zu threads, the computation time is: %g\n", thread_count, timer);

    // cout << G.population_points().size() << endl;
    print_vector(G.population_points().back());

    // cout << G.offline.back() << endl;
    // cout << G.bests.back() << endl;
    // cout << G.online[iterations-1] << endl;
    // cout << G.offline[iterations-1] << endl;
    // cout << G.bests[iterations-1] << endl;
    // print_vector(G.online, true);
    // print_vector(G.offline, true);
    // print_vector(G.bests, true);
    // file.open("canonical_dispersion.txt");
    // G.plot_population_points(file);
    return 0;
}
