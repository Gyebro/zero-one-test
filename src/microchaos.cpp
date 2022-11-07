#include "microchaos.h"
#include <cmath>

double Coth(double v) { // Hyperbolic cotangent
    return 1.0/tanh(v);
}

double Csch(double v) { // Hyperbolic cosecant
    return 1.0/sinh(v);
}

double sgn(double x) {
    if (x == 0.0) return 0;
    else return ((x > 0) ? 1.0 : -1.0);
}

MicroChaosParamDomain::MicroChaosParamDomain(double alpha, double delta) : a(alpha), d(delta) {
    pmin    = pow(a,2);
    pmax    = pow(a,2)*(1 - (8*(pow(a,2) + pow(a,2)*pow(d,2))*exp(2*a*d))/(2*pow(a,2) + 2*pow(a,2)*pow(d,2) + 2*(pow(a,2) + pow(a,2)*pow(d,2))*exp(2*a*d) - 2*(pow(a,2) + pow(a,2)*pow(d,2))*exp(a*d)*(exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2))))));
    dmin    = -(a*d) - (2*sqrt(pow(a,2) + pow(a,2)*pow(d,2))*exp(a*d))/(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))) + sqrt(pow(a,2) + pow(a,2)*pow(d,2))/(exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))*(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2))))) + (sqrt(pow(a,2) + pow(a,2)*pow(d,2))*exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2))))/(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2))));
    dswitch = ((-4*pow(pow(a,2) + pow(a,2)*pow(d,2),1.5)*exp(3*a*d))/(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))) + (-2*pow(a,2) - 2*pow(a,2)*pow(d,2))*(a*d - (sqrt(pow(a,2) + pow(a,2)*pow(d,2))*(exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))))/(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2))))) - 2*(pow(a,2) + pow(a,2)*pow(d,2))*exp(2*a*d)*(-3*a*d + (sqrt(pow(a,2) + pow(a,2)*pow(d,2))*(exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))))/(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2))))) + (2*sqrt(pow(a,2) + pow(a,2)*pow(d,2))*exp(a*d)*(pow(a,2) + pow(a,2)*pow(d,2) + 3*(pow(a,2) + pow(a,2)*pow(d,2)) + a*d*sqrt(pow(a,2) + pow(a,2)*pow(d,2))*(-exp(-2*sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(2*sqrt(pow(a,2) + pow(a,2)*pow(d,2)))) - (pow(a,2) + pow(a,2)*pow(d,2))*(exp(-2*sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(2*sqrt(pow(a,2) + pow(a,2)*pow(d,2))))))/(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))))/(2*pow(a,2) + 2*pow(a,2)*pow(d,2) + 2*(pow(a,2) + pow(a,2)*pow(d,2))*exp(2*a*d) - 2*(pow(a,2) + pow(a,2)*pow(d,2))*exp(a*d)*(exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))));
    dmax    = -(a*d) + (2*sqrt(pow(a,2) + pow(a,2)*pow(d,2))*(exp(a*d) + (exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2))))/2.))/(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2))));
}
double MicroChaosParamDomain::a1condP(double D) const {
    return (pow(a,2)*exp(a*d)*(sqrt(pow(a,2) + pow(a,2)*pow(d,2))*(-exp(-(a*d)) + exp(a*d)) + (D*(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))))/2.))/(-sqrt(pow(a,2) + pow(a,2)*pow(d,2)) + exp(a*d)*(-(a*d*(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))))/2. + (sqrt(pow(a,2) + pow(a,2)*pow(d,2))*(exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))))/2.));
}
double MicroChaosParamDomain::a0condP(double D) const {
    return (pow(a,2)*((D*(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))))/2. - sqrt(pow(a,2) + pow(a,2)*pow(d,2))*((exp(-(a*d)) + exp(a*d))/2. + (exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2))))/2.)))/((sqrt(pow(a,2) + pow(a,2)*pow(d,2))*(-exp(-(a*d)) + exp(a*d)))/2. - (a*d*(-exp(-sqrt(pow(a,2) + pow(a,2)*pow(d,2))) + exp(sqrt(pow(a,2) + pow(a,2)*pow(d,2)))))/2.);
}
double MicroChaosParamDomain::a1condD(double P) const {
    return (-(P*d) + sqrt(1 + pow(d,2))*Csch(a*sqrt(1 + pow(d,2)))*(-(P*cosh(a*d)) + P*cosh(a*sqrt(1 + pow(d,2))) + (P - 2*pow(a,2))*sinh(a*d)))/a;
}
double MicroChaosParamDomain::a0condD(double P) const {
    return (-(P*d) + sqrt(1 + pow(d,2))*Csch(a*sqrt(1 + pow(d,2)))*(pow(a,2)*(cosh(a*d) + cosh(a*sqrt(1 + pow(d,2)))) + P*sinh(a*d)))/a;
}
double MicroChaosParamDomain::getBestP() const {
    return (pow(a,2)*(1 - 2*exp(a*d)*cosh(a*sqrt(1 + pow(d,2)))))/(1 + exp(2*a*d) - 2*exp(a*d)*cosh(a*sqrt(1 + pow(d,2))));
}
double MicroChaosParamDomain::getBestD() const {
    return (a*(-d + 2*d*exp(a*d)*cosh(a*sqrt(1 + pow(d,2))) +  sqrt(1 + pow(d,2))*(Coth(a*sqrt(1 + pow(d,2))) - exp(a*d)*cosh(2*a*sqrt(1 + pow(d,2)))*Csch(a*sqrt(1 + pow(d,2))))))/(1 + exp(2*a*d) - 2*exp(a*d)*cosh(a*sqrt(1 + pow(d,2))));
}
std::vector<vec2> MicroChaosParamDomain::getStableParams(const double r, const double padding) {
    std::vector<vec2> pars;
    double d;
    for (d = dmin; d <= dswitch; d+=r) {
        for (double p = pmin+padding; p <= a1condP(d); p+=r) {
            pars.emplace_back(vec2({p, d}));
        }
    }
    for (; d <= dmax; d+=r) {
        for (double p = pmin+padding; p <= a0condP(d); p+=r) {
            pars.emplace_back(vec2({p, d}));
        }
    }
    return pars;
}
std::vector<vec2> MicroChaosParamDomain::getStableParamsConstD(const double d, const double r, const double padding) {
    std::vector<vec2> pars;
    if (d < dswitch) {
        for (double p = pmin+padding; p <= a1condP(d); p+=r) {
            pars.emplace_back(vec2({p, d}));
        }
    } else {
        for (double p = pmin+padding; p <= a0condP(d); p+=r) {
            pars.emplace_back(vec2({p, d}));
        }
    }
    return pars;
}
std::vector<vec2> MicroChaosParamDomain::getStableParamsConstP(const double p, const double r, const double padding) {
    std::vector<vec2> pars;
    double dmin = a1condD(p);
    double dmax = a0condD(p);
    for (double d = dmin+padding; d <= dmax-padding; d+=r) {
        pars.emplace_back(vec2({p, d}));
    }
    return pars;
}


MicroChaosMapStatic::MicroChaosMapStatic(double P, double D, double alpha, double delta) {
    MicroChaosMapStatic::alpha = alpha;
    MicroChaosMapStatic::delta = delta;
    MicroChaosMapStatic::P = P;
    MicroChaosMapStatic::D = D;
    Gamma = sqrt(1.0+delta*delta);
    U = Us(1.0);
    b = bs(1.0);
    k[0] = P;
    k[1] = D;
}

mat2 MicroChaosMapStatic::Us(double s) const {
    mat2 Uss;
    double chags = cosh(alpha*Gamma*s);
    double shags = sinh(alpha*Gamma*s);
    double eads = exp(alpha*delta*s);
    double emads = 1.0/eads;
    Uss(0,0) = (delta*shags+Gamma*chags)*emads/Gamma;
    Uss(0,1) = (shags/alpha)*emads/Gamma;
    Uss(1,0) = (alpha*shags)*emads/Gamma;
    Uss(1,1) = (Gamma*chags-delta*shags)*emads/Gamma;
    return Uss;
}

vec2 MicroChaosMapStatic::bs(double s) const {
    vec2 bss;
    double chags = cosh(alpha*Gamma*s);
    double shags = sinh(alpha*Gamma*s);
    double eads = exp(alpha*delta*s);
    double emads = 1.0/eads;
    bss[0]   = (Gamma-emads*(Gamma*chags+delta*shags))/(Gamma*alpha*alpha);
    bss[1]   = -shags*emads/(alpha*Gamma);
    return bss;
}

vec2 MicroChaosMapStatic::step(const vec2& y0) const {
    double M;
    M = SymmetricFloor(k*y0); // Rounding at output
    vec2 y1 = U*y0 + b*M;
    return y1;
}

MicroChaosMapStaticRoundIn::MicroChaosMapStaticRoundIn(double P, double D, double alpha, double delta) {
    MicroChaosMapStaticRoundIn::alpha = alpha;
    MicroChaosMapStaticRoundIn::delta = delta;
    MicroChaosMapStaticRoundIn::P = P;
    MicroChaosMapStaticRoundIn::D = D;
    Gamma = sqrt(1.0+delta*delta);
    U = Us(1.0);
    b = bs(1.0);
    k[0] = P;
    k[1] = D;
}

mat2 MicroChaosMapStaticRoundIn::Us(double s) const {
    mat2 Uss;
    double chags = cosh(alpha*Gamma*s);
    double shags = sinh(alpha*Gamma*s);
    double eads = exp(alpha*delta*s);
    double emads = 1.0/eads;
    Uss(0,0) = (delta*shags+Gamma*chags)*emads/Gamma;
    Uss(0,1) = (shags/alpha)*emads/Gamma;
    Uss(1,0) = (alpha*shags)*emads/Gamma;
    Uss(1,1) = (Gamma*chags-delta*shags)*emads/Gamma;
    return Uss;
}

vec2 MicroChaosMapStaticRoundIn::bs(double s) const {
    vec2 bss;
    double chags = cosh(alpha*Gamma*s);
    double shags = sinh(alpha*Gamma*s);
    double eads = exp(alpha*delta*s);
    double emads = 1.0/eads;
    bss[0]   = (Gamma-emads*(Gamma*chags+delta*shags))/(Gamma*alpha*alpha);
    bss[1]   = -shags*emads/(alpha*Gamma);
    return bss;
}

vec2 MicroChaosMapStaticRoundIn::step(const vec2& y0) const {
    double M;
    vec2 y0r = {SymmetricFloor(y0[0]), SymmetricFloor(y0[1])};
    M = k*y0r; // Rounding at input
    vec2 y1 = U*y0 + b*M;
    return y1;
}

double MicroChaosMapStaticRoundIn::getP() const {
    return P;
}

void MicroChaosMapStaticRoundIn::setP(double p) {
    P = p;
    k[0] = P;
}

double MicroChaosMapStaticRoundIn::getD() const {
    return D;
}

void MicroChaosMapStaticRoundIn::setD(double d) {
    D = d;
    k[1] = D;
}
