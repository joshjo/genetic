#include <cstdlib>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include <algorithm>
#include "functions.h"

typedef unsigned long size_n;

struct Attributes {
    int rows;
    int cols;
    size_n N;
    size_n l_x;
    size_n l_y;
    double prob_crossing;
    double prob_mutation;

    size_n array_size;

    size_n get_array_size() {
        array_size = N * (l_x + l_y);
        return array_size;
    }
};


double limit_x [] = {-5.12, 5.12};
double limit_y [] = {-5.12, 5.12};


__host__ __device__ double drop_wave(double x, double y){
    return -((1+cos(12*sqrt(pow(x, 2)+pow(y, 2))))/(0.5*pow(x, 2)+pow(y, 2)+2));
}

__constant__ Attributes attrs;

__host__ __device__ double get_real(
    bool * array,
    int length,
    double left,
    double right,
    size_n start,
    size_n end
){
    int sum = 0;
    int j = 0;
    for(size_n i = start; i < end; ++i, ++j) {
        // printf("%d:", array[i]);
        sum += array[i] * pow(2, length - j - 1);
    }
    // printf("\n");
    return left + sum * ((right - left) / (pow(2, length) - 1));
}

__device__ float generate( curandState* globalState, int ind )
{
    curandState localState = globalState[ind];
    float RANDOM = curand_uniform( &localState );
    globalState[ind] = localState;
    return RANDOM;
}

__global__ void setup_kernel ( curandState * state, unsigned long seed )
{
    int id = threadIdx.x;
    curand_init( seed, id, 0, &state[id] );
}


__device__ int chooseNeighbor (int idx, curandState* st, int cols, int rows, int c, int r) {
    int k = generate(st, idx) * 100000;
    k = (k / 100) % 4;
    size_n odx = 0;

    if (k == 0) {
        if ((c + 1) < cols) {
            odx = idx + 1;
        } else {
            odx = idx - cols + 1;
        }
    } else if (k == 1) {
        if ((c - 1) >= 0) {
            odx = idx - 1;
        } else {
            odx = idx + cols - 1;
        }
    } else if (k == 2) {
        if ((r + 1) < cols) {
            odx = idx + cols;
        } else {
            odx = idx % cols;
        }
    } else if (k == 3) {
        if ((r - 1) >= 0) {
            odx = idx - cols;
        } else {
            odx = (cols * (rows - 1)) + idx;
        }
    }
    return odx;
}


__global__ void tournament(bool* a, bool* o, curandState* st) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    int cols = attrs.cols;
    int rows = attrs.rows;
    int r = idx / rows;
    int c = idx % cols;
    if (idx >= attrs.N) {
        return;
    }

    int odx = chooseNeighbor(idx, st, cols, rows, c, r);

    int ai = idx * (attrs.l_x + attrs.l_y);
    int oi = odx * (attrs.l_x + attrs.l_y);

    double ax = get_real(a, attrs.l_x, -5.12, 5.12, ai, ai + attrs.l_x);
    double ay = get_real(
        a, attrs.l_y, -5.12, 5.12, ai + attrs.l_x, ai + attrs.l_x + attrs.l_y);

    double ox = get_real(a, attrs.l_x, -5.12, 5.12, oi, oi + attrs.l_x);
    double oy = get_real(
        a, attrs.l_y, -5.12, 5.12, oi + attrs.l_x, oi + attrs.l_x + attrs.l_y);

    // printf("idx: %d (%f, %f) -> odx: %d (%f, %f)\n", idx, ax, ay, odx, ox, oy);
    double fa = drop_wave(ax, ay);
    double fo = drop_wave(ox, oy);

    int start;
    if (fa > fo) {
        start = ai;
    } else {
        start = oi;
    }
    for (int i = 0; i < (attrs.l_x + attrs.l_y); i++) {
        o[ai + i] = a[start + i];
    }
}

__global__ void reproduction(bool * a, bool * o, curandState* st) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    int cols = attrs.cols;
    int rows = attrs.rows;
    int r = idx / rows;
    int c = idx % cols;
    if (idx >= attrs.N) {
        return;
    }

    double p_mutation = generate(st, idx);
    double p_crossing = generate(st, idx);

    if (p_mutation < attrs.prob_mutation) {
        int k = generate(st, idx) * (attrs.l_x + attrs.l_y) * 100;
        k = k % (attrs.l_x + attrs.l_y);
        int ai = idx * (attrs.l_x + attrs.l_y);
        o[ai + k] = !a[ai + k];
    }
    if (p_crossing < attrs.prob_crossing) {
        int k = generate(st, idx) * (attrs.l_x + attrs.l_y) * 100;
        k = k % (attrs.l_x + attrs.l_y);
        int odx = chooseNeighbor(idx, st, cols, rows, c, r);
        int ai = idx * (attrs.l_x + attrs.l_y);
        int ao = odx * (attrs.l_x + attrs.l_y);

        for (int i = 0; i < (attrs.l_x + attrs.l_y); i++) {
            if (i < k) {
                o[ai + i] = a[ao + i];
            } else {
                o[ai + i] = a[ai + i];
            }

        }
    }
}


void printMatrix(bool * & mat, Attributes & attrs, int limit = 10) {
    size_n inc = attrs.l_x + attrs.l_y;
    int i = 0;
    for (int x = 0; x < attrs.array_size; x += inc, i++) {
        if (i > limit) {
            break;
        }
        double realx = get_real(
            mat, attrs.l_x, limit_x[0], limit_x[1], x, x + attrs.l_x);
        double realy = get_real(
            mat, attrs.l_y, limit_y[0], limit_y[1], x + attrs.l_x, x + attrs.l_x + attrs.l_y);
        printf("%d: %f, %f\n", i, realx, realy);
    }
}


void random_bitvector(bool * & arr, size_n size) {
    for (size_n i = 0; i < size; i++) {
        int number = rand() % 2;
        arr[i] = number;
    }
}

void printBitVector(bool * & arr, size_n size, int limit = 10) {
    for (size_n i = 0; i < size; i++) {
        if (i > limit) break;
        printf("%d - ", arr[i]);
    }
}


int main () {
    srand((unsigned)time(0));
    int iterations = 100;

    size_n cols = 31;
    size_n rows = 31;
    size_n precission = 10;
    Attributes A;
    A.cols = cols;
    A.rows = rows;
    A.prob_crossing = 0.75;
    A.prob_mutation = 0.10;
    A.N = cols * rows;
    A.l_x = get_l(limit_x[0], limit_x[1], precission);
    A.l_y = get_l(limit_y[0], limit_y[1], precission);

    size_n arraySize = A.get_array_size();

    bool * a = new bool [arraySize];
    bool * x = new bool [arraySize];
    bool * da, * dx;

    size_n size = (arraySize) * sizeof(bool);

    cudaMemcpyToSymbol(attrs, &A, sizeof(Attributes));
    cudaMalloc((void **) &da, size);
    cudaMalloc((void **) &dx, size);

    curandState *d_state;

    cudaMalloc(&d_state, A.N * sizeof(curandState));

    random_bitvector(a, arraySize);


    int maxThreads = (A.N > 128) ? 128 : A.N;
    int blocks = (A.N + maxThreads - 1) / maxThreads;

    printMatrix(a, A);

    cudaMemcpy(da, a, size, cudaMemcpyHostToDevice);
    cudaMemcpy(dx, x, size, cudaMemcpyHostToDevice);

    setup_kernel <<< blocks, maxThreads >>> (d_state, unsigned(time(NULL)) );

    for (int i = 0; i < iterations; i++) {
        tournament <<< blocks, maxThreads >>> (da, dx, d_state);
        cudaMemcpy(da, dx, size, cudaMemcpyDeviceToDevice);
        reproduction <<< blocks, maxThreads >>> (da, dx, d_state);
        cudaMemcpy(da, dx, size, cudaMemcpyDeviceToDevice);
    }

    cudaMemcpy(a, da, size, cudaMemcpyDeviceToHost);

    printf("\n\n");
    printMatrix(a, A);

    printf("%f\n", drop_wave(-5.120000, -5.120000));

    return 0;
}
