#include "spline.hh"
#include "../base/mdrange.hh"

using namespace Conan;

//=============================================================================
Array<double> Spline<1>::calc_coeff(iVector<1> const &S) const
{
	Array<double> b(4);
	auto const &dx = box->dx;

	auto g = [&] (iVector<R> const &X)
	{
		return data[box->idx(X)];
	};

	for (unsigned i = 0; i < 2; ++i)
	{
		iVector<R>	X = S + box->block[i];
		b[i]   		= g(X);
		b[i+2] 		= (g(X + dx[0]) - g(X - dx[0]))/2;
	}

	Array<double> alpha(4, 0.0);
	MdRange<2> Z(iVector<2>(4));

	for (unsigned i = 0; i < 16; ++i)
	{
		iVector<2> z = Z[i];
		alpha[z[1]] += b[z[0]] * A[i];
	}

	return alpha;
}

double Spline<1>::operator()(dVector<1> const &x) const
{
	return f(x);
}

double Spline<1>::f(dVector<1> const &x_) const
{
	dVector<R> x = x_ / box->res();
	iVector<R> origin = floor_cast(x);
	dVector<R> P = x - dVector<R>(origin);

	auto c = calc_coeff(origin);
	MdRange<R> loop(iVector<R>(4));

	// P(x, y, z) = Sum_i Sum_j Sum_k [c_ijk * x^i y^j z^k]
	double s = 0;
	for (unsigned i = 0; i < 4; ++i)
	{
		iVector<R> e = loop[i];
		s += c[i] * pow(P[0], e[0]);
	}

	return s;
}

dVector<1> Spline<1>::df(dVector<1> const &x_) const
{
	dVector<R> x = x_ / box->res();
	iVector<R> origin = floor_cast(x);
	dVector<R> P = x - dVector<R>(origin);

	auto c = calc_coeff(origin);
	MdRange<R> loop(iVector<R>(4));

	dVector<R> v(0);
	for (unsigned n = 0; n < 4; ++n)
	{
		iVector<R> e = loop[n];

		for (unsigned k = 0; k < R; ++k)
		{
			unsigned i = (k + 1) % R;

			if (e[k] != 0)
			{
				v[k] += c[n] * e[k] * pow(P[k], e[k] - 1);
			}
		}
	}

	return v;
}

std::pair<double, dVector<1>> Spline<1>::fdf(dVector<1> const &x_) const
{
	dVector<R> x = x_ / box->res();
	iVector<R> origin = floor_cast(x);
	dVector<R> P = x - dVector<R>(origin);

	auto c = calc_coeff(origin);
	MdRange<R> loop(iVector<R>(4));

	double s = 0;
	dVector<R> v(0);
	for (unsigned n = 0; n < 4; ++n)
	{
		iVector<R> e = loop[n];
		s += c[n] * pow(P[0], e[0]) * pow(P[1], e[1]);

		for (unsigned k = 0; k < R; ++k)
		{
			unsigned i = (k + 1) % R;

			if (e[k] != 0)
			{
				v[k] += c[n] * e[k] * pow(P[k], e[k] - 1);
			}
		}
	}

	return std::make_pair(s, v);
}
//=============================================================================

//=============================================================================
Array<double> Spline<2>::calc_coeff(iVector<2> const &S) const
{
	Array<double> b(16);
	auto const &dx = box->dx;

	auto g = [&] (iVector<R> const &X)
	{
		return data[box->idx(X)];
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
	MdRange<2> Z(iVector<2>(16));

	for (unsigned i = 0; i < 256; ++i)
	{
		iVector<2> z = Z[i];
		alpha[z[1]] += b[z[0]] * A[i];
	}

	return alpha;
}

double Spline<2>::operator()(dVector<2> const &x) const
{
	return f(x);
}

double Spline<2>::f(dVector<2> const &x_) const
{
	dVector<R> x = x_ / box->res();
	iVector<R> origin = floor_cast(x);
	dVector<R> P = x - dVector<R>(origin);

	auto c = calc_coeff(origin);
	MdRange<R> loop(iVector<R>(4));

	// P(x, y, z) = Sum_i Sum_j Sum_k [c_ijk * x^i y^j z^k]
	double s = 0;
	for (unsigned i = 0; i < 16; ++i)
	{
		iVector<R> e = loop[i];
		s += c[i] * pow(P[0], e[0]) * pow(P[1], e[1]) * pow(P[2], e[2]);
	}

	return s;
}

dVector<2> Spline<2>::df(dVector<2> const &x_) const
{
	dVector<R> x = x_ / box->res();
	iVector<R> origin = floor_cast(x);
	dVector<R> P = x - dVector<R>(origin);

	auto c = calc_coeff(origin);
	MdRange<R> loop(iVector<R>(4));

	dVector<R> v(0);
	for (unsigned n = 0; n < 16; ++n)
	{
		iVector<R> e = loop[n];

		for (unsigned k = 0; k < R; ++k)
		{
			unsigned i = (k + 1) % R;

			if (e[k] != 0)
			{
				v[k] += c[n] * e[k] * pow(P[k], e[k] - 1) 
				      * pow(P[i], e[i]);
			}
		}
	}

	return v;
}

std::pair<double, dVector<2>> Spline<2>::fdf(dVector<2> const &x_) const
{
	dVector<R> x = x_ / box->res();
	iVector<R> origin = floor_cast(x);
	dVector<R> P = x - dVector<R>(origin);

	auto c = calc_coeff(origin);
	MdRange<R> loop(iVector<R>(4));

	double s = 0;
	dVector<R> v(0);
	for (unsigned n = 0; n < 16; ++n)
	{
		iVector<R> e = loop[n];
		s += c[n] * pow(P[0], e[0]) * pow(P[1], e[1]);

		for (unsigned k = 0; k < R; ++k)
		{
			unsigned i = (k + 1) % R;

			if (e[k] != 0)
			{
				v[k] += c[n] * e[k] * pow(P[k], e[k] - 1) 
				      * pow(P[i], e[i]);
			}
		}
	}

	return std::make_pair(s, v);
}
//=============================================================================

//=============================================================================
Array<double> Spline<3>::calc_coeff(iVector<3> const &S) const
{
	Array<double> b(64);
	auto const &dx = box->dx;

	auto g = [&] (iVector<R> const &X)
	{
		return data[box->idx(X)];
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
	MdRange<2> Z(iVector<2>(64));

	for (unsigned i = 0; i < 4096; ++i)
	{
		iVector<2> z = Z[i];
		alpha[z[1]] += b[z[0]] * A[i];
	}

	return alpha;
}

double Spline<3>::operator()(dVector<3> const &x) const
{
	return f(x);
}

double Spline<3>::f(dVector<3> const &x_) const
{
	dVector<R> x = x_ / box->res();
	iVector<R> origin = floor_cast(x);
	dVector<R> P = x - dVector<R>(origin);

	auto c = calc_coeff(origin);
	MdRange<R> loop(iVector<R>(4));

	// P(x, y, z) = Sum_i Sum_j Sum_k [c_ijk * x^i y^j z^k]
	double s = 0;
	for (unsigned i = 0; i < 64; ++i)
	{
		iVector<R> e = loop[i];
		s += c[i] * pow(P[0], e[0]) * pow(P[1], e[1]) * pow(P[2], e[2]);
	}

	return s;
}

dVector<3> Spline<3>::df(dVector<3> const &x_) const
{
	dVector<R> x = x_ / box->res();
	iVector<R> origin = floor_cast(x);
	dVector<R> P = x - dVector<R>(origin);

	auto c = calc_coeff(origin);
	MdRange<R> loop(iVector<R>(4));

	dVector<R> v(0);

	for (unsigned n = 0; n < 64; ++n)
	{
		iVector<R> e = loop[n];

		for (unsigned k = 0; k < R; ++k)
		{
			unsigned i = (k + 1) % R;
			unsigned j = (k + 2) % R;

			if (e[k] != 0)
			{
				v[k] += c[n] * e[k] * pow(P[k], e[k] - 1) 
				      * pow(P[i], e[i]) * pow(P[j], e[j]);
			}
		}
	}

	return v;
}

std::pair<double, dVector<3>> Spline<3>::fdf(dVector<3> const &x_) const
{
	dVector<R> x = x_ / box->res();
	iVector<R> origin = floor_cast(x);
	dVector<R> P = x - dVector<R>(origin);

	auto c = calc_coeff(origin);
	MdRange<R> loop(iVector<R>(4));

	double s = 0;
	dVector<R> v(0);
	for (unsigned n = 0; n < 64; ++n)
	{
		iVector<R> e = loop[n];
		s += c[n] * pow(P[0], e[0]) * pow(P[1], e[1]) * pow(P[2], e[2]);

		for (unsigned k = 0; k < R; ++k)
		{
			unsigned i = (k + 1) % R;
			unsigned j = (k + 2) % R;

			if (e[k] != 0)
			{
				v[k] += c[n] * e[k] * pow(P[k], e[k] - 1) 
				      * pow(P[i], e[i]) * pow(P[j], e[j]);
			}
		}
	}

	return std::make_pair(s, v);
}
//=============================================================================

