// HEADER FILE 
// BESSEL FUNCTIONS OF FIRST AND SECOND KINDS FOR COMPLEX ARGUMENTS 

#ifndef HEADER
#define HEADER


#define egamma_v  0.5772156649L

template<typename Tp>
    struct numeric_constants
    {
        static Tp pi() throw()
        { return static_cast<Tp>(3.1415926535897932384626433832795029L); }
        static Tp pi_2() throw()
        { return static_cast<Tp>(1.5707963267948966192313216916397514L); }
    };

/*
Returns the cylindrical Bessel functions of order nu : J_{nu} or  I_{nu} by series expansion. 
    nu  The order of the Bessel function.
    x   The argument of the Bessel function.
    sgn  The sign of the alternate terms
        -1 for the Bessel function of the first kind.
        +1 for the modified Bessel function of the first kind.
    max_iter maximum number of series terms
    return  The output Bessel function.
*/
template <typename Tp>
std::complex<Tp> cyl_bessel_ij_series(Tp nu, std::complex<Tp> x, Tp sgn, unsigned int max_iter)
{
    Tp flag;
    flag = nu < Tp(0) && (std::fmod(trunc(nu),Tp(2)) != Tp(0)) ? Tp(-1) : Tp(1); // Checks for 
    // J_{-nu}(x) = (-1)^{nu} J_{nu}(x)             if nu is an integer,
    // Complex logarithm calculations consistency   if nu is a non-integer.
    if (floor(nu)==nu)
        nu = std::abs(nu);
    if (x == std::complex<Tp> (0,0))
	    return nu == Tp(0) ? Tp(1) : Tp(0);
    std::complex<Tp> x2 = x / Tp(2);
    std::complex<Tp> fact = nu * std::log(x2);
    fact -= lgamma(nu + Tp(1));
    fact = flag * std::exp(fact);
    std::complex<Tp> xx4 = sgn * x2 * x2;
    std::complex<Tp> Jn = Tp(1);
    std::complex<Tp> term = Tp(1);
    for (unsigned int i = 1; i < max_iter; ++i)
        {
            term *= xx4 / (Tp(i) * (nu + Tp(i)));
            Jn += term;
            if (std::abs(term /Jn) < std::numeric_limits<Tp>::epsilon())
                break;
        }
    return fact * Jn ;
}

/*
Computes the asymptotic cylindrical Bessel and Neumann functions of order nu: J_{nu} ,N_{nu} .
    nu  The order of the Bessel functions.
    x   The argument of the Bessel functions.
    Jnu  The output Bessel function of the first kind.
    Nnu  The output Neumann function (Bessel function of the second kind).
*/
template<typename Tp>
void cyl_bessel_jn_asymp(Tp nu, std::complex<Tp> x, std::complex<Tp>& Jnu, std::complex<Tp>& Nnu)
    {
        const Tp mu   = Tp(4) * nu * nu;
        const Tp m1 = mu - Tp(1);
        const Tp m9 = mu - Tp(9);
        const Tp m25 = mu - Tp(25);
        const Tp m49 = mu - Tp(49);
        std::complex<Tp> xx = Tp(64) * x * x;
        std::complex<Tp> P = Tp(1) - m1 * m9 / (Tp(2) * xx) * (Tp(1) - m25 * m49 / (Tp(12) * xx));
        std::complex<Tp> Q = m1 / (Tp(8) * x) * (Tp(1) - m9 * m25 / (Tp(6) * xx));
        std::complex<Tp> chi = x - (nu + Tp(0.5L)) * numeric_constants<Tp>::pi_2();
        std::complex<Tp> c = std::cos(chi);
        std::complex<Tp> s = std::sin(chi);
        std::complex<Tp> coef = std::sqrt(Tp(2) / (numeric_constants<Tp>::pi() * x));
        Jnu = coef * (c * P - s * Q);
        Nnu = coef * (s * P + c * Q);
        return;
    }

/*
Returns the Bessel function of order mu: J_{mu}(z) .
    mu  The order of the Bessel function.
    z   The argument of the Bessel function.
    return  The output Bessel function.
*/
template <typename Tp>
std::complex<Tp> bessel_j(Tp mu, std::complex<Tp> z)
{ 
    if (std::isnan(mu) || std::isnan(z.real()) || std::isnan(z.imag()) )
        return std::numeric_limits<Tp>::quiet_NaN();  
    else if (std::abs(z)> Tp(50))
    {
        std::complex<Tp> J_nu, N_nu;
        cyl_bessel_jn_asymp(mu, z, J_nu, N_nu);
        return J_nu;
    }
    else 
        return cyl_bessel_ij_series(mu, z, Tp(-1), 200);        
}

/*
Returns the neumann function of order nu : N_{nu} by series expansion. 
    nu  The order of the Bessel function.
    x   The argument of the Bessel function.
    max_iter maximum number of series terms
    return  The output Bessel function.
*/
template<typename Tp>
std::complex<Tp> neumann_series(Tp nu, std::complex<Tp> x, unsigned int max_iter)
{
    if (x==std::complex<Tp>(0,0)) 
        {
            if (nu <= Tp(0) && std::fmod(trunc(nu),Tp(2)) != Tp(0))
                return std::numeric_limits<Tp>::infinity();
            else
                return -std::numeric_limits<Tp>::infinity();   
        }       
    else
    {
        if (floor(nu)==nu)
            {   
                Tp flag;
                flag = nu < Tp(0) && (std::fmod(nu,Tp(2)) != Tp(0)) ? Tp(-1) : Tp(1);
                nu = std::abs(nu);
                std::complex<Tp> x2 = x/Tp(2);
                std::complex<Tp> fact_ = nu * std::log(x2);
                std::complex<Tp> fact = fact_ - lgamma(nu + Tp(1));
                fact = std::exp(fact);
                std::complex<Tp> xx4 = Tp(-1)*x2*x2;
                std::complex<Tp> term = Tp(1);
                std::complex<Tp> Jn = Tp(1);
                Tp t1 = Tp(1);
                Tp fact1 = Tp(1);
                Tp t = -Tp(2)*egamma_v;
                for (unsigned int i =1; i<nu; i++)
                    {
                        fact1 *= Tp(i);  
                        t += Tp(1)/Tp(i); 
                    }
                std::complex<Tp> _Nn = Tp(-1);
                std::complex<Tp> term1 = Tp(1);
                t += Tp(1)/nu;
                std::complex<Tp> fact2 = std::exp(fact_)/(fact1*nu);
                std::complex<Tp> Nn_ = -t;
                Tp tt = Tp(0);
                for (unsigned int i = 1; i<max_iter; i++)
                    {
                        term *= xx4 / (Tp(i) * (nu + Tp(i)));
                        Jn += term;
                        if (i < nu)
                            {
                                term1 *= x2*x2/((nu - Tp(i)) * Tp(i));
                                _Nn -= term1; 
                            }
                        tt +=  (Tp(1) / (nu + Tp(i)) + Tp(1) / Tp(i));
                        Nn_ -= (term * (t +  tt));  
                    }
                Jn *= Tp(2)*fact*(fact_/nu);
                _Nn *= (fact1/std::exp(fact_));  
                Nn_ *= fact2; 
                return flag * (Tp(1)/ numeric_constants<Tp>::pi())* (Jn + Nn_ + _Nn);
            }
        else
            {
                return (bessel_j(nu, x)*std::cos(numeric_constants<Tp>::pi()*nu)-bessel_j(-nu, x))/std::sin(numeric_constants<Tp>::pi()*nu);
            }
    }   
}

/*
Returns the neumann function (Bessel function of the second kind) of order mu: N_{mu}(z) .
    mu  The order of the neumann function .
    z   The argument of the neumann function.
    return  The output neumann function.
*/
template <typename Tp>
std::complex<Tp> neumann(Tp mu, std::complex<Tp> z)
{
    if (std::isnan(mu) || std::isnan(z.real()) || std::isnan(z.imag()) )
        return std::numeric_limits<Tp>::quiet_NaN();  

    else if (std::abs(z)> Tp(50))
    {
        std::complex<Tp> J_nu, N_nu;
        cyl_bessel_jn_asymp(mu, z, J_nu, N_nu);
        return N_nu;
    }
    else 
        return neumann_series(mu, z, 200);
}

#endif 