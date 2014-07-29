#pragma once
#include "system.hh"

namespace Conan
{
	template <unsigned R>
	class Limits: public std::pair<dVector<R>, dVector<R>>
	{
		typedef std::pair<dVector<R>, dVector<R>> Base;

		public:
			using Base::Base;

			std::pair<Limits, Limits> split(unsigned dim, double p)
			{
				dVector<R> a = Base::first,
				           b = Base::second;

				a[dim] = p;
				b[dim] = p;

				return make_pair(Limits(Base::first, b), 
				                 Limits(a, Base::second));
			}
	};

	template <unsigned R>
	inline Function<R>
	argument(Function<R> phi, double D, dVector<R> const &x)
	{
		return [&] (dVector<R> const &q)
		{
			return -q.dot(q)/2 + D * phi(q) + q.dot(x);
		}
	}

	template <unsigned R>
	std::pair<dVector<R>, double>
	find_maximum(BoxPtr<R> box, Array<double> phi, dVector<R> const &x, Limits<R> lim)
	{
		double res = box->res();
		iVector<R> a = ceil_cast(lim->first / res),
		           b = floor_cast(lim->second / res),
			   shape = b - a;

		MdRange<R> search_range(shape);
		auto func = map(search_range, [&] (iVector<R> const &qi)
		{
			dVector<R> q = dVector<R>(qi + a) * res;
			return -q.dot(q)/2 + D * phi[box->idx(qi + a)] + q.dot(x);
		});

		auto max = std::max_element(func.begin(), func.end());
		dVector<R> q = dVector<R>(max.arg() + a) * res;

	}

	template <unsigned R>
	std::pair<Array<dVector<R>>, Array<double>>
	flt(BoxPtr<R> L, BoxPtr<R> E, Array<double> phi, double D)
	{
		// progressively split
		unsigned N = E->N();
	}
	
	template <unsigned R>
	void run_nv(std::string const &id, BoxMaker L_, BoxMaker E_,
		Array<double> phi, FILE_FORMAT fmt)
	{
	}
}

