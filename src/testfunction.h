#ifndef TESTFUNCTION_H
#define TESTFUNCTION_H

#include "includes.h"
#define PI 3.141592653589793238462643383279502884L


struct TestFunction{
    opt_function foo;
    vector<double> lefts;
    vector<double> rights;
    double optimal;
    string name;

    TestFunction(opt_function f, string n, vector<double> l, vector<double> r,
            double o): foo(f), name(n), lefts(l), rights(r), optimal(o){}
};

double ackley(vector<double> x){
    double a = 20;
    double b = 0.2;
    double c = 2*PI;
    double first_sum = 0;
    double second_sum = 0;
    int n = x.size();
    for(int i=0; i<n; ++i){
        first_sum += pow(x[i], 2);
        second_sum += cos(c*x[i]);
    }
    return -a*exp(-b*sqrt(first_sum/n))-exp(second_sum/n)+a+exp(1);
}

double bukin6(vector<double> x){
    return 100*sqrt(abs(x[1]-0.01*pow(x[0],2)))+0.01*abs(x[0]+10);
}

double cross_in_tray(vector<double> x){
    return -0.0001*pow(abs(sin(x[0])*sin(x[1])*exp(abs(100-sqrt(pow(x[0], 2)+pow(x[1], 2))/PI )) )+1, 0.1);
}

double cross_leg_table(vector<double> x){
    return -1/(pow(abs(exp(abs(100-(sqrt(pow(x[0], 2)*pow(x[1], 2))/PI)))*sin(x[0])*sin(x[1]))+1, 0.1));
}

double dejong5(vector<double> x){
    vector< vector<double> > A = {{-32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32}, {-32, -32, -32, -32, -32, -16, -16, -16, -16, -16, 0, 0, 0, 0, 0, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32}};
    double resp = 0;
    for(int i=0; i<25; ++i){
        double a0i = A[0][i];
        double a1i = A[1][i];
        int term1 = i+1;
        double term2 = pow(x[0]-a0i, 6);
        double term3 = pow(x[1]-a1i, 6);
        resp += float(1)/(term1+term2+term3);
    }
    return 0.002 + float(1)/(0.002 + resp);
}

double rastrigin(vector<double> x){
    int d = x.size();
    int resp = 10*d;
    for(int i=0; i<d; ++i){
        resp += pow(x[i], 2)-10*cos(2*PI*x[i]);
    }
    return resp;
}

double damavandi(vector<double> x){
    return (1-pow(abs((sin(PI*(x[0]-2))*sin(PI*(x[1]-2)))/(pow(PI, 2)*(x[0]-2)*(x[1]-2))),5))*(2+pow(x[0]-7, 2)+pow(x[1]-7, 2));
}

double drop_wave(vector<double> x){
    return -((1+cos(12*sqrt(pow(x[0], 2)+pow(x[1], 2))))/(0.5*pow(x[0], 2)+pow(x[1], 2)+2));
}

double easom(vector<double> x){
    return -cos(x[0])*cos(x[1])*exp(-pow(x[0]-PI, 2)-pow(x[1]-PI, 2));
}

double fn6(vector<double> x){
    return 0.5-(pow(sin(sqrt(pow(x[0],2)+pow(x[1], 2))), 2)-0.5)/pow(1+.001*(pow(x[0], 2)+pow(x[1], 2)), 2);
}

double goldsteinprice(vector<double> x){
    return (1+pow(x[0]+x[1]+1, 2)*(
        19-14*x[0]+3*pow(x[0],2)-14*x[1]+6*x[0]*x[1]+3*pow(x[1],2)))*(
        30+pow(2*x[0]-3*x[1], 2)*(
            18-32*x[0]+12*pow(x[0], 2)+48*x[1]-36*x[0]*x[1]+27*pow(x[1],2))
    );
}

double levy13(vector<double> x){
    return pow(sin(3*PI*x[0]), 2)+pow(x[0]-1, 2)*(1+pow(sin(3*PI*x[1]), 2))+pow(x[1]-1, 2)*(1+pow(sin(2*PI*x[1]), 2));
}

double michalewicz(vector<double> x){
    int d = x.size();
    double m = 10;
    double r = 0;
    for(int i=0; i<d; ++i)
        r += sin(x[i])*pow(sin((i+1)*pow(x[i], 2)/PI),2*m);
    return -r;
}

double rosenbrock(vector<double> x){
    return 100*pow((x[1]-pow(x[0], 2)), 2)+pow(x[0]-1, 2);
}

double schaffer2(vector<double> x){
    return 0.5+(pow(sin(pow(x[0], 2)-pow(x[1], 2)), 2)-0.5)/(
        1+0.001*(pow(x[0], 2)+pow(x[1], 2)));
}

double schwefel(vector<double> x){
    double tmp = 0;
    for(int i=0; i<x.size(); ++i){
        tmp += x[i]*sin(sqrt(abs(x[i])));
    }
    return (418.9829*2)-tmp;
}

double sphere(vector<double> x){
    return pow(x[0], 2) + pow(x[1], 2);
}

double sqsums(vector<double> x){
    return 1*pow(x[0], 2)+2*pow(x[1], 2);
}

double zakharov(vector<double> x){
    return (pow(x[0], 2)+pow(x[1], 2))+pow(0.5*1*pow(x[0], 2)+0.5*2*pow(x[1], 2), 2)+pow(0.5*1*pow(x[0], 2)+0.5*2*pow(x[1], 2), 4);
}

#endif // TESTFUNCTION_H
