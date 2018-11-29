#include <pybind11/pybind11.h>
#include <cmath>

static inline double polevl(double x, double coef[], int N)
{
    double ans;
    int i;
    double *p;

    p = coef;
    ans = *p++;
    i = N;

    do
	ans = ans * x + *p++;
    while (--i);

    return (ans);
}

static inline double p1evl(double x, double coef[], int N)
{
    double ans;
    double *p;
    int i;

    p = coef;
    ans = x + *p++;
    i = N - 1;

    do
	ans = ans * x + *p++;
    while (--i);

    return (ans);
}

/* A[]: Stirling's formula expansion of log Gamma
 * B[], C[]: log Gamma function between 2 and 3
 */
static double A[] = {
    8.11614167470508450300E-4,
    -5.95061904284301438324E-4,
    7.93650340457716943945E-4,
    -2.77777777730099687205E-3,
    8.33333333333331927722E-2
};

static double B[] = {
    -1.37825152569120859100E3,
    -3.88016315134637840924E4,
    -3.31612992738871184744E5,
    -1.16237097492762307383E6,
    -1.72173700820839662146E6,
    -8.53555664245765465627E5
};

static double C[] = {
    /* 1.00000000000000000000E0, */
    -3.51815701436523470549E2,
    -1.70642106651881159223E4,
    -2.20528590553854454839E5,
    -1.13933444367982507207E6,
    -2.53252307177582951285E6,
    -2.01889141433532773231E6
};

/* log( sqrt( 2*pi ) ) */
const double LS2PI = 0.91893853320467274178;
const double LOGPI = log(M_PI);

#define MAXLGM 2.556348e305


double lgam_sgn(double x, int *sign)
{
    double p, q, u, w, z;
    int i;

    *sign = 1;

    if (!isfinite(x))
	return x;

    if (x < -34.0) {
	q = -x;
	w = lgam_sgn(q, sign);
	p = floor(q);
	if (p == q) {
	  lgsing:
	    return (INFINITY);
	}
	i = p;
	if ((i & 1) == 0)
	    *sign = -1;
	else
	    *sign = 1;
	z = q - p;
	if (z > 0.5) {
	    p += 1.0;
	    z = p - q;
	}
	z = q * sin(M_PI * z);
	if (z == 0.0)
	    goto lgsing;
	/*     z = log(NPY_PI) - log( z ) - w; */
	z = LOGPI - log(z) - w;
	return (z);
    }

    if (x < 13.0) {
	z = 1.0;
	p = 0.0;
	u = x;
	while (u >= 3.0) {
	    p -= 1.0;
	    u = x + p;
	    z *= u;
	}
	while (u < 2.0) {
	    if (u == 0.0)
		goto lgsing;
	    z /= u;
	    p += 1.0;
	    u = x + p;
	}
	if (z < 0.0) {
	    *sign = -1;
	    z = -z;
	}
	else
	    *sign = 1;
	if (u == 2.0)
	    return (log(z));
	p -= 2.0;
	x = x + p;
	p = x * polevl(x, B, 5) / p1evl(x, C, 6);
	return (log(z) + p);
    }

    if (x > MAXLGM) {
	return (*sign * INFINITY);
    }

    q = (x - 0.5) * log(x) - x + LS2PI;
    if (x > 1.0e8)
	return (q);

    p = 1.0 / (x * x);
    if (x >= 1000.0)
	q += ((7.9365079365079365079365e-4 * p
	       - 2.7777777777777777777778e-3) * p
	      + 0.0833333333333333333333) / x;
    else
	q += polevl(p, A, 4) / x;
    return (q);
}

double lgam(double x)
{
    int sign;
    return lgam_sgn(x, &sign);
}

double logbin(int k, int n, double p) {
  double combiln = (lgam(n+1) - (lgam(k+1) + lgam(n-k+1)));
  return combiln + k * log(p) + (n - k) * log1p(-p);
}

double binomial_(int k, int n, double p) {
  return exp(logbin(k, n, p));
}

namespace py = pybind11;
PYBIND11_MODULE(binomial, m) {
        m.doc() = "pybind11 example plugin"; // optional module docstring
        m.def("binomial", &binomial_,
              py::arg("k"), py::arg("n"), py::arg("p"));
}
