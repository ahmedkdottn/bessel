Motivation:
The C++ standard library header <cmath> uses only reals as arguments of the bessel functions (https://en.cppreference.com/w/cpp/numeric/special_functions). 
Although useful and quite accurate (https://www.boost.org/doc/libs/1_67_0/libs/math/doc/html/math_toolkit/bessel/bessel_first.html), the boost library will not accept complex arguments for the bessel functions neither. 
We present hereby an implementation of the Bessel functions of the first and second kinds for complex arguments. 
We provide also some test cases in the "test.cpp" file.
  
-------------------------------------------------------------------------------
Please note: 
While the implementation provides a good approximation of the exact result for most cases, it is yet to be improved (for large arguments and orders in particular). 
  
-------------------------------------------------------------------------------
How to use:

-Please use #include "bessel.tcc" in your main file. 

-The file bessel.tcc contains templatized functions for 

the Bessel function of the first kind,

and the Neumann function (Bessel function of the second kind).

-For integer orders, please use the "Tp" type. e.g: bessel_j(3.L,std::complex<long double>(-3.3L,9.3L))   
