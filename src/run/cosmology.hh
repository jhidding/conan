#pragma once
#include <cmath>

namespace Conan
{
	class Cosmology
	{
		double H0, Omega_m, Omega_L, Omega_k, grav_cst;

		public:
			Cosmology(double H0_, double m_, double L_):
				H0(H0_), Omega_m(m_), Omega_L(L_),
				Omega_k(1 - m_ - L_), grav_cst(3./2 * m_ * H0_ * H0_) 
			{}

			double operator()(double a) const
			{
				return H0 * a * sqrt(
					Omega_L + 
					Omega_m * pow(a, -3) +
					Omega_k * pow(a, -2));
			}

			double G() const
			{
				return grav_cst;
			}
	};
}

