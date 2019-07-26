#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

double sqr(double x){
    return x*x;
}

double math_expected_value(double *arr, int n)
{
    int i;
    double E=0;
    for(i=0;i<n;i++){
        E+=*(arr+i);	
    }
    return E/n;
}

double std_deviation(double E, double *arr, int n)
{

    int i;
    double D;
    D=0;
    for(i=0;i<n;i++){
        D+=sqr(E - *(arr+i));
    }
    return sqrt(D/n);
}


double uniform(double a, double b)
{  
  // return rand() / (RAND_MAX + 1.0) * (b - a) + a;
    return (genrand_int32()+0.5)/4294967296.0 * (b - a) + a;
}

double gauss(double mu,double sigma)
{
    double x1, x2, w, y1, y2;


    do {
        x1 = 2.0 * uniform(0,1) - 1.0;
        x2 = 2.0 * uniform(0,1) - 1.0;
        w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    return mu+sigma*y1;	
}

double rgamma(double a,double b) {
    double d,c,x,v,u;
    if(a>=1){
        d = a-1./3.; 
        c = 1./sqrt(9.*d);
        while(1){
            do {
                x=gauss(0,1.0); 
                v=1.+c*x;
            } while(v<=0.);
            v = v * v * v; 
            u = uniform(0,1);
            if( u < 1.0-0.0331*(x*x)*(x*x) ){
                return d*v*b;
            }
            if( log(u) < 0.5*x*x+d*(1.0-v+log(v)) ){
                return d*v*b;
            }
        }
    } else {
        x = rgamma(a+1,b);
        x = x * pow(uniform(0,1), 1.0/a); 
        return x;
    }
}

double gamma_math_expected_value(double a,double b)
{
    return a*b;
}
double gamma_std_deviation(double a,double b)
{
    return sqrt(a*b*b);
}

#ifdef MAIN
int main()
{
    srand(time(NULL));
    int i;
    int n=1000;
    double arr[n];
    for(i=0;i<1000;i++){
        arr[i] = rgamma(4, 1);
        printf("%f\n",arr[i]);
    }
    double E,D;
    E = math_expected_value(arr,n);
    D = std_deviation(E, arr, n);
    printf("----------------------\n");
    printf("E=%f&&D=%f\n",E,D);
    printf("----------------------\n");
    printf("E=%f&&D=%f\n",gamma_math_expected_value(4, 1),gamma_std_deviation(4, 1));
    return 0;
}
#endif // MAIN
