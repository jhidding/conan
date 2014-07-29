#pragma once
#include "system.hh"

namespace Conan
{
	template <typename Q, unsigned R>
	class SplineBase
	{
		BoxPtr<R> 	box;
		Q		f;

		public:
			Spline(BoxPtr<R> box_, Q f_):
				box(box_), f(f_)
			{}

			typename Q::value_type operator()(dVector<R> const &x_) const
			{
				dVector<R> x = x_ / box->res();
				iVector<R> origin = floor_cast(x);
				dVector<R> A[2];

				transform(x, origin, A[1], [] (double q, int o) 
					{ return q - o; });
				transform(A[1], A[0], [] (double g)
					{ return 1 - g; });

				typename Q::value_type v(0);

				for (auto dx : box->block)
				{
					double z = 1;
					for (unsigned k = 0; k < R; ++k) 
						z *= A[dx[k]][k];

					v += f[box->idx(origin+dx)] * z;
				}

				return v;
			}
	};

	template <typename Q>
	class Spline<Q, 2>: public SplineBase<2>
	{
		enum { R = 2 };

		BoxPtr<R> 	box;
		Q		f;

		public:
			Spline(BoxPtr<R> box_, Q f_):
				box(box_), f(f_)
			{}
			
			Array<double> calc_coeff(iVector<R> const &S)
			{
				Array<double> b(16);
				auto const &dx = box->dx;

				auto g = [&] (iVector<R> const &X)
				{
					return f[box->idx(X)];
				};

				for (unsigned i = 0; i < 4; ++i)
				{
					iVector<R>	X = S + box->block[i];

					b[i]		= g(X);


					b[i + 4]	= (g(X + dx[0]) - g(X - dx[0])) / 2;

					b[i + 8]	= (g(X + dx[1]) - g(X - dx[1])) / 2;

					b[i + 12]	= ((g(X + dx[1] + dx[0]) - g(X - dx[1] + dx[0])) -
					                   (g(X + dx[1] - dx[0]) - g(X - dx[1] - dx[0]))) / 4;
				}

				Array<double> alpha(16, 0.0);
				MdRange<R> Z(iVector<R>({16, 16}));

				for (unsigned i = 0; i < 256; ++i)
				{
					iVector<R> z = Z[i];
					alpha[z[1]] += b[z[0]] * A[i];
				}

				return alpha;
			}

			double operator()(dVector<R> const &x) const
			{
				return f(x);
			}

			double f(dVector<R> const &x) const
			{
			}

			dVector<R> df(dVector<R> const &x) const
			{
			}

			std::pair<double, dVector<R>> fdf(dVector<R> const &x) const
			{
			}
	};

	template <typename Q>
	class Spline<Q, 3>
	{
		enum { R = 3 };

		BoxPtr<R> 	box;
		Q		f;

		public:
			Spline(BoxPtr<R> box_, Q f_):
				box(box_), f(f_)
			{}

			Array<double> calc_coeff(iVector<R> const &S)
			{
				Array<double> b(64);
				auto const &dx = box->dx;

				auto g = [&] (iVector<R> const &X)
				{
					return f[box->idx(X)];
				};

				for (unsigned i = 0; i < 8; ++i)
				{
					iVector<R>	X = S + box->block[i];

					b[i]		= g(X);

					b[i + 8]	= (g(X + dx[0]) - g(X - dx[0])) / 2;

					b[i + 16]	= (g(X + dx[1]) - g(X - dx[1])) / 2;

					b[i + 24]	= (g(X + dx[2]) - g(X - dx[2])) / 2;

					b[i + 32]	= ((g(X + dx[1] + dx[0]) - g(X - dx[1] + dx[0])) -
					                   (g(X + dx[1] - dx[0]) - g(X - dx[1] - dx[0]))) / 4;

					b[i + 40]	= ((g(X + dx[2] + dx[1]) - g(X - dx[2] + dx[1])) -
					                   (g(X + dx[2] - dx[1]) - g(X - dx[2] - dx[1]))) / 4;

					b[i + 48]	= ((g(X + dx[0] + dx[2]) - g(X - dx[0] + dx[2])) -
					                   (g(X + dx[0] - dx[2]) - g(X - dx[0] - dx[2]))) / 4;
					
					b[i + 56]	= (g(X + dx[0] + dx[1] + dx[2]) - g(X + dx[0] + dx[1] - dx[2])
							-  g(X + dx[0] - dx[1] + dx[2]) + g(X + dx[0] - dx[1] - dx[2])
							-  g(X - dx[0] + dx[1] + dx[2]) + g(X - dx[0] + dx[1] - dx[2])
							+  g(X - dx[0] - dx[1] + dx[2]) - g(X - dx[0] - dx[1] - dx[2])) / 8;
				}

				Array<double> alpha(64, 0.0);
				MdRange<R> Z(iVector<R>({64, 64}));

				for (unsigned i = 0; i < 4096; ++i)
				{
					iVector<R> z = Z[i];
					alpha[z[1]] += b[z[0]] * A[i];
				}

				return alpha;
			}

			double operator()(dVector<R> const &x) const
			{
				return f(x);
			}

			double f(dVector<R> const &x) const
			{
			}

			dVector<R> df(dVector<R> const &x) const
			{
			}

			std::pair<double, dVector<R>> fdf(dVector<R> const &x) const
			{
			}
	};
}

