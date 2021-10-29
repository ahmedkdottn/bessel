#include <iostream>
#define _USE_MATH_DEFINES
#include <tr1/cmath>
#include <complex>
#include "bessel.tcc"

int main()
{
    std::cout << bessel_j(0.,std::complex<double>(0,0)) << std::endl;
    std::cout << bessel_j(-3.5,std::complex<double>(0,6.5)) << std::endl;
    std::cout << bessel_j(3.f,std::complex<float>(-3.5f,-0.9f)) << std::endl;
    std::cout << bessel_j(-24.5,std::complex<double>(30.5,30)) << std::endl;
    
    std::cout << neumann(2.,std::complex<double>(0,-4.05)) << std::endl;
    std::cout << neumann(-60.5,std::complex<double>(13.45,45.25)) << std::endl;
    std::cout << neumann(7.,std::complex<double>(0,0)) << std::endl;
}