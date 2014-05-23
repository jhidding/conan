#pragma once
#include "system.hh"

namespace Conan
{
	template <typename Q, unsigned R>
	class LinearInterpolation
	{
		BoxPtr<R>	box;
		Q		f;

		public:
			LinearInterpolation(BoxPtr<R> box_, Q f_):
				box(box_), f(f_)
			{}

			typename Q::value_type operator()(dVector<R> const &x_) const
			{
				dVector<R> x = x_ / box->res();
				iVector<R> origin;
				dVector<R> A[2];
				//double r = box->scale();

				for (unsigned k = 0; k < R; ++k)
				{
					origin[k] = (int)(x[k]);
					A[1][k] = x[k] - origin[k];
					A[0][k] = 1 - A[1][k];
				}

				size_t s = box->idx(origin);
				typename Q::value_type v(0);

				for (unsigned i = 0; i < (1 << R); ++i)
				{
					double z = 1;
					for (unsigned k = 0; k < R; ++k) 
						z *= A[box->block[i][k]][k];

					v += f[box->idx(origin + box->block[i])] * z;
				}

				return v;
			}
	};


	template <unsigned R>
	void gradient_2nd_order(BoxPtr<R> box, Array<dVector<R>> result, Array<double> f)
	{
		System::MdRange<R> N(iVector<R>(box->N()));

		for (size_t x = 0; x < box->size(); ++x)
		{
			iVector<R> p0 = N[x];
			dVector<R> v(0);

			for (unsigned k = 0; k < R; ++k)
			{
				size_t  iy = box->idx(p0 - box->dx[k]*2),
					iz = box->idx(p0 - box->dx[k]),
					ia = box->idx(p0 + box->dx[k]),
					ib = box->idx(p0 + box->dx[k]*2);

				v[k] = f[iy]/12. - 2*f[iz]/3. + 2*f[ia]/3. - f[ib]/12.;
			}

			result[x] = v / box->res();
		}
	}

	template <unsigned R>
	class Gradient
	{
		BoxPtr<R> box;
		Array<double> f;
		System::MdRange<R> N;

		public:
			typedef dVector<R> value_type;

			Gradient(BoxPtr<R> box_, Array<double> f_):
				box(box_), f(f_), N(iVector<R>(box->N()))
			{}

			dVector<R> operator[](size_t i) const
			{
				iVector<R> p0 = N[i];
				dVector<R> v(0);

				for (unsigned k = 0; k < R; ++k)
				{
					size_t  iy = box->idx(p0 - box->dx[k]*2),
						iz = box->idx(p0 - box->dx[k]),
						ia = box->idx(p0 + box->dx[k]),
						ib = box->idx(p0 + box->dx[k]*2);

					v[k] = f[iy]/12. - 2*f[iz]/3. + 2*f[ia]/3. - f[ib]/12.;
				}

				return v / box->res();
			}	
	};


	template <unsigned R>
	void md_cic(BoxPtr<R> box, Array<double> result, Array<dVector<R>> X, double m = 1.0)
	{
		System::fill(result, 0.0);

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
				#pragma omp critical
				{
					result[idx] += z;
				}
			}
		}
	}
}

