/* fourier.h
 *
 * contains methods for fourier space filtering
 */

#pragma once
#include <functional>
#include <cmath>
#include "fft.hh"
#include "mvector.hh"
#include "mdrange.hh"
#include "map.hh"

namespace Fourier
{
	typedef std::complex<double> complex64;
	extern complex64 math_i;

	template <unsigned R>
	using dVector = System::mVector<double, R>;

	inline std::function<double (size_t)> wave_number(size_t N, double L)
	{
		return [N, L] (size_t i)
		{
			long a(i);
			return (i > N/2 ? a - long(N) : a) * (2 * M_PI / L);
		};
	}

	inline std::function<double (std::complex<double> const &)> real_part(double size)
	{
		return [size] (std::complex<double> const &z)
		{
			return z.real() / size;
		};
	}

	template <unsigned R>
	using KSpace = System::Map<System::MdRange<R>, System::mVector<double, R>>;

	template <unsigned R>
	KSpace<R> kspace(unsigned N, double L)
	{
		typename KSpace<R>::arg_type shape(N);
		System::MdRange<R> X(shape);

		return KSpace<R>(X, [N, L] (typename KSpace<R>::arg_type const &x)
		{ 
			typename KSpace<R>::value_type k; 
			for (unsigned i = 0; i < R; ++i) 
				k[i] = (unsigned(x[i]) > N/2 ? x[i] - long(N) : x[i]) * (2 * M_PI / L);
			return k; 
		});
	}

	typedef std::function<double (double)> RealFunction;

	template <unsigned R>
	using FilterBase = std::function<complex64 (dVector<R> const &)>;

	template <unsigned R>
	class Filter: public FilterBase<R>
	{
		public:
			Filter(FilterBase<R> const &f):
				FilterBase<R>(f) {}

			Filter operator*(Filter const &o) const
			{
				Filter p(*this);
				return Filter([p, o] (dVector<R> const &K)
				{
					return p(K) * o(K);
				});
			}	
	};

	template <unsigned R>
	inline Filter<R> scale(double t)
	{
		return Filter<R>([t] (dVector<R> const &K)
		{
			double v = 1.0;
			for (auto k : K)
			{
				v *= exp(t * (cos(k) - 1));
			}
			return complex64(v);
		});
	}

	template <unsigned R>
	inline Filter<R> power_spectrum(RealFunction const &P)
	{
		return Filter<R>([P] (dVector<R> const &K)
		{
			return complex64(sqrt(P(K.norm())));
		});
	}

	template <unsigned R>
	inline Filter<R> potential()
	{
		return Filter<R>([] (dVector<R> const &K)
		{
			return complex64(-1. / K.sqr());
		});
	}

	template <unsigned R>
	inline Filter<R> derivative(unsigned i)
	{
		return Filter<R>([i] (dVector<R> const &K)
		{
			return math_i * sin(K[i]);
		});
	}

	template <unsigned R>
	inline std::function<complex64 (complex64, dVector<R> const &)>
	filter(Filter<R> const &f)
	{
		return [f] (complex64 z, dVector<R> const &K)
		{
			return z * f(K);
		};
	}
}

// vim:ts=4:sw=4:tw=80
