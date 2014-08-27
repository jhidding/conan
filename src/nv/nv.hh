#pragma once
#include <queue>
#include <fstream>
#include "system.hh"
#include "minimizor.hh"
#include "spline.hh"

namespace Conan
{
	template <unsigned R>
	using Limits = std::pair<dVector<R>, dVector<R>>;

	/*template <unsigned R>
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
	};*/

	template <unsigned R>
	class FLT_iter_pars
	{
		typedef std::pair<Maybe<FLT_iter_pars>, Maybe<FLT_iter_pars>> pair;

		public:
			iVector<R>	x, min, max;

			FLT_iter_pars(iVector<R> const &x_, iVector<R> const &min_, iVector<R> const &max_):
				x(x_), min(min_), max(max_) {}
			
			pair split(unsigned sdim) const
			{
				if (max[sdim] - min[sdim] <= 2)
				{
					return pair(Nothing, Nothing);
				}

				if (x[sdim] - min[sdim] <= 1)
				{
					iVector<R> nx = x, nmin = min;
					nx[sdim] = (x[sdim] + max[sdim])/2;
					nmin[sdim] = x[sdim];

					return pair(Just<FLT_iter_pars>(nx, nmin, max), Nothing);
				}

				if (max[sdim] - x[sdim] <= 1) // shouldn't happen, just in case
				{
					iVector<R> nx = x, nmax = max;
					nx[sdim] = (x[sdim] + min[sdim])/2;
					nmax[sdim] = x[sdim];

					return pair(Just<FLT_iter_pars>(nx, min, nmax), Nothing);
				}

				if (x[sdim] == 0) // the point is the first one in this
						 // direction, special action is needed
				{
					iVector<R> nx = x, nmax = max, nmin = min;
					nx[sdim] = max[sdim];
					nmax[sdim] = max[sdim]*2;
					nmin[sdim] = 0;

					return pair(Just<FLT_iter_pars>(nx, nmin, nmax), Nothing);
				}

				iVector<R> b_min = min, a_max = max,
					   a_x   = x,   b_x   = x;

				a_max[sdim] = x[sdim];
				b_min[sdim] = x[sdim];

				a_x[sdim] = (x[sdim] + min[sdim]) / 2;
				b_x[sdim] = (x[sdim] + max[sdim]) / 2;

				return pair(Just<FLT_iter_pars>(a_x, min, a_max),
					    Just<FLT_iter_pars>(b_x, b_min, max));
			}

			Limits<R> search_limits(BoxPtr<R> E, Array<dVector<R>> psi) const
			{
				dVector<R> q_min, q_max;

				for (unsigned k = 0; k < R; ++k)
				{
					iVector<R> a = x, b = x;
					a[k] = min[k]; b[k] = max[k];

					dVector<R> qa = dVector<R>(a)*E->res() + psi[E->idx(a)],
					           qb = dVector<R>(b)*E->res() + psi[E->idx(b)];

					q_min[k] = qa[k]; q_max[k] = qb[k];
				}

				return Limits<R>(q_min, q_max);
			}
	};

	template <unsigned R>
	class Maximant
	{
		BoxPtr<R>	box;
		Spline<R> 	phi;
		double		D;

		dVector<R>	x;

		public:
			Maximant(BoxPtr<R> box_, Array<double> phi_, double D_, dVector<R> const &x_):
				box(box_), phi(box_, phi_), D(D_), x(x_) {}

			double operator()(dVector<R> const &q) const
			{
				double v = f(q);
				return v;
			}

			double f(dVector<R> const &q) const
			{
				return q.dot(q)/2 - D * phi.f(q) - q.dot(x);
			}

			dVector<R> df(dVector<R> const &q) const
			{
				return q - phi.df(q) * D - x;
			}

			std::pair<double, dVector<R>> fdf(dVector<R> const &q) const
			{
				auto a = phi.fdf(q);
				return std::make_pair(
					q.dot(q)/2 - a.first * D - q.dot(x),
					q - a.second * D - x);
			}
	};

	template <unsigned R>
	std::pair<double, dVector<R>>
	find_maximum(BoxPtr<R> box, Array<double> phi, double D, dVector<R> const &x, Limits<R> lim, bool enhance = true)
	{
		double res = box->res();
		iVector<R> a = floor_cast(lim.first / res),
		           b = ceil_cast(lim.second / res),
			   shape = b - a;

		MdRange<R> search_range(shape);
		auto func = map(search_range, [&] (iVector<R> const &qi)
		{
			dVector<R> q = dVector<R>(qi + a) * res;
			return q.dot(q)/2 - D * phi[box->idx(qi + a)] - q.dot(x);
		});

		auto intres = std::min_element(func.begin(), func.end());
		dVector<R> q = dVector<R>(intres.arg() + a) * res;
		
		if (not enhance)
			return std::make_pair(*intres, q);

		Maximant<R> M(box, phi, D, x);
		Minimizor<Maximant<R>, R> minimize(M);
		minimize.set(q);

		unsigned i = 0;
		while (minimize.should_continue() and i < 100)
		{
			minimize.iterate();
			++i;
		}

		return std::make_pair(minimize.value(), minimize.pos());
	}

	template <unsigned R>
	std::pair<Array<double>, Array<dVector<R>>>
	flt(BoxPtr<R> L, BoxPtr<R> E, Array<double> phi_0, double D, bool enhance = true)
	{
		typedef FLT_iter_pars<R> Fip;

		// progressively split: take a position x in the middle of
		// the domain; search limits are determined by taking the found
		// previous result. Every new point splits the search erea for
		// the next ones in two.
	
		Array<dVector<R>> psi(E->size(), dVector<R>(0.));	
							// store the displacement, this
							// cures all ails with periodics
		Array<double>	  phi(E->size());	// potential
		
		int N = E->N();

		// first make some borders
		auto search_fip = [&] (ptr<Fip> f)
		{
			size_t i = E->idx(f->x);
			dVector<R> x = dVector<R>(f->x) * E->res();

			auto limits = f->search_limits(E, psi);
			auto z = find_maximum(L, phi_0, D, x, limits, enhance);

			phi[i] = (z.first + (x.dot(x))/2) / D;
			psi[i] = z.second - x;
		};

		// 0-th element
		std::vector<std::queue<ptr<Fip>>> search_queue(R);
		search_queue[0].push(make_ptr<Fip>(iVector<R>(0), iVector<R>(-N/2), iVector<R>(N/2)));

		int dim = 0;
		while (dim >= 0)
		{
			auto front = search_queue[dim].front();
			search_queue[dim].pop();
			search_fip(front);

			for (unsigned k = dim; k < R; ++k)
			{
				auto S = front->split(k);
				if (S.first)
					search_queue[k].push(S.first);

				if (S.second)
					search_queue[k].push(S.second);

				if (not search_queue[k].empty())
					dim = k;
			}

			while (dim >= 0 and search_queue[dim].empty())
				--dim;
		}

		return std::make_pair(phi, psi);
	}


	template <typename T, unsigned R>
	void write_matrix_array_txt(std::ostream &fo, BoxPtr<R> box, Array<T> data)
	{
		for (size_t i = 0; i < box->size(); ++i)
		{
			fo << data[i] << " ";
			
			for (unsigned k = 0; k < R-1; ++k)
			{
				if (box->I[i][k] == (box->shape()[k] - 1))
				{	
					fo << std::endl;
				}
			}
		}
	}

	template <typename T, unsigned R>
	void write_list_array_txt(std::ostream &fo, BoxPtr<R> box, Array<T> data)
	{
		for (size_t i = 0; i < box->size(); ++i)
		{
			fo << box->G[i] << " " << data[i] << "\n";
		}
	}

	template <unsigned R>
	void nv_run(std::string const &id, BoxMaker L_, BoxMaker E_,
		Array<double> phi, double D, FILE_FORMAT fmt, bool enhance)
	{
		auto E = E_.box<R>(), L = L_.box<R>();
		auto result = flt<R>(L, E, phi, D, enhance);

		std::ofstream fo;
		switch (fmt)
		{
			case FMT_ASCII:
				fo.open(System::timed_filename(id, "flt", D, ".txt"));

				fo << "# potential\n";
				write_matrix_array_txt(fo, E, result.first);

				fo << "\n\n# eulerian map\n";	
				write_list_array_txt(fo, E, result.second);

				fo.close();
				break;

			case FMT_CONAN:
				fo.open(System::timed_filename(id, "flt", D));

				System::save_to_file(fo, result.first, "potential");
				System::save_to_file(fo, result.second, "eulerian_map");

				fo.close();
				break;
		}
	}
}

