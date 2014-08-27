#pragma once
#include "system.hh"

namespace Conan
{
	template <unsigned R>
	class Spline;

	template <>
	class Spline<1>
	{
		enum { R = 1 };

		static double const A[16];
		BoxPtr<R> 	box;
		Array<double>	data;

		public:
			Spline(BoxPtr<R> box_, Array<double> f_):
				box(box_), data(f_)
			{}
			
			Array<double> calc_coeff(iVector<R> const &S) const;
			double operator()(dVector<R> const &x) const;
			double f(dVector<R> const &x) const;
			dVector<R> df(dVector<R> const &x) const;
			std::pair<double, dVector<R>> fdf(dVector<R> const &x) const;
	};

	template <>
	class Spline<2>
	{
		enum { R = 2 };

		static double const A[256];
		BoxPtr<R> 	box;
		Array<double>	data;

		public:
			Spline(BoxPtr<R> box_, Array<double> f_):
				box(box_), data(f_)
			{}
			
			Array<double> calc_coeff(iVector<R> const &S) const;
			double operator()(dVector<R> const &x) const;
			double f(dVector<R> const &x) const;
			dVector<R> df(dVector<R> const &x) const;
			std::pair<double, dVector<R>> fdf(dVector<R> const &x) const;
	};

	template <>
	class Spline<3>
	{
		enum { R = 3 };

		static double const A[4096];	
		BoxPtr<R> 	box;
		Array<double>	data;

		public:
			Spline(BoxPtr<R> box_, Array<double> f_):
				box(box_), data(f_)
			{}

			Array<double> calc_coeff(iVector<R> const &S) const;
			double operator()(dVector<R> const &x) const;
			double f(dVector<R> const &x) const;
			dVector<R> df(dVector<R> const &x) const;
			std::pair<double, dVector<R>> fdf(dVector<R> const &x) const;
	};
}

