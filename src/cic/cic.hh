#pragma once
#include "../base/system.hh"
#include "../base/progress.hh"
#include "../run/numeric.hh"

namespace Conan
{
	template <unsigned R, typename Container>
	void md_cic_T(BoxPtr<R> box, Array<double> result, Container X, double m = 1.0)
	{
		System::fill(result, 0.0);
		Misc::ProgressBar pb(X.size());
		#pragma omp parallel for
		for (size_t i = 0; i < X.size(); ++i)
		{
			dVector<R> p0 = (X[i] % box->L()) / box->res();
			iVector<R> origin;
			dVector<R> A[2];

			for (unsigned k = 0; k < R; ++k)
			{
				origin[k] = (int)(p0[k]);
				A[1][k] = p0[k] - origin[k];
				A[0][k] = 1 - A[1][k];
			}

			size_t i0 = box->idx(origin);
			for (unsigned j = 0; j < (1 << R); ++j)
			{
				double z = m;
				for (unsigned k = 0; k < R; ++k) 
					z *= A[box->block[j][k]][k];

				size_t idx = box->idx(origin + box->block[j]);
				#pragma omp critical (A)
				{
					result[idx] += z;
				}
			}

			#pragma omp critical (B)
			{
				pb.tic();
			}
		}
		pb.finish();
	}

	template <unsigned R>
	void cic_run(std::istream &fi, std::ostream &fo, System::BoxMaker b_m_, 
		     System::BoxMaker b_f_, unsigned s)
	{
		BoxPtr<R> b_m = b_m_.box<R>(),
		          b_f = b_f_.box<R>();
		
		Array<dVector<R>> X = System::load_from_file<dVector<R>>(fi, "pos");
		Array<double> rho(b_f->size(), 0.0);

		if (s == 1)
		{
			double a = pow(double(b_f->N())/double(b_m->N()), -float(R));
			md_cic_T(b_f, rho, X);
		}
		else
		{
			for (size_t i = 0; i < X.size(); ++i)
				X[i] -= b_m->G[i];

			LinearInterpolation<Array<dVector<R>>, R> X_li(b_m, X);
			auto grid = System::grid_space<R>(b_m->N() * s, b_m->L());
			auto X_mp = System::map(grid, [&] (dVector<R> const &y)
			{
				return y + X_li(y);
			});

			double a = pow(double(s * b_m->N())/double(b_f->N()), -float(R));
			md_cic_T(b_f, rho, X_mp, a);
		}

		System::save_to_file(fo, rho, "density");
	}
}

