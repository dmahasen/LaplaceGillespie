#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iostream>

unsigned int GetUint(unsigned int& u, unsigned int& v)
{
    v = 36969*(v & 65535) + (v >> 16);
    u = 18000*(u & 65535) + (u >> 16);
    return (v << 16) + u;
}



double GetUniform(unsigned int& u, unsigned int& v)
{
    // 0 <= u <= 2^32
    unsigned int z = GetUint(u, v);
    // The magic number is 1/(2^32 + 1) and so result is positive and less than 1.
    return z*2.328306435996595e-10;
}


// Get normal (Gaussian) random sample with specified mean and standard deviation
double GetNormal(double mean = 0.0, double standardDeviation = 1.0)
{
    if (standardDeviation <= 0.0)
    {
		std::stringstream os;
        os << "Standard deviation must be positive." << "\n" 
           << "Received standard deviation " << standardDeviation;
        throw std::invalid_argument( os.str() );
    }
	    // Use Box-Muller algorithm
    double u1 = GetUniform();
    double u2 = GetUniform();
    double r = sqrt( -2.0*log(u1) );
    double theta = 2.0*PI*u2;
    return mean + standardDeviation*r*sin(theta);
}

double  rgamma(double shape, double scale)
{
    // Implementation based on "A Simple Method for Generating Gamma Variables"
    // by George Marsaglia and Wai Wan Tsang.  ACM Transactions on Mathematical Software
    // Vol 26, No 3, September 2000, pages 363-372.

    double d, c, x, xsquared, v, u;

    if (shape >= 1.0)
    {
        d = shape - 1.0/3.0;
        c = 1.0/sqrt(9.0*d);
        for (;;)
        {
            do
            {
                x = GetNormal();
                v = 1.0 + c*x;
            }
            while (v <= 0.0);
            v = v*v*v;
            u = GetUniform();
            xsquared = x*x;
            if (u < 1.0 -.0331*xsquared*xsquared || log(u) < 0.5*xsquared + d*(1.0 - v + log(v)))
                return scale*d*v;
        }
    }
    else if (shape <= 0.0)
    {
		std::stringstream os;
        os << "Shape parameter must be positive." << "\n"
           << "Received shape parameter " << shape;
        throw std::invalid_argument( os.str() );
    }
    else
    {
        double g = rgamma(shape+1.0, 1.0);
        double w = GetUniform();
        return scale*g*pow(w, 1.0/shape);
    }
}
