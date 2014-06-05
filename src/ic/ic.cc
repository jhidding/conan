#include "ic.hh"
#include "../base/fourier.hh"
#include "../base/access.hh"

using namespace System;
using namespace Conan;

template <typename A>
typename A::value_type mean(A const &a)
{
	return std::accumulate(a.begin(), a.end(), typename A::value_type(0)) / a.size();
}

template <unsigned R>
typename Fourier<R>::Filter CDM(Header const &H)
{
	double const 
		e 		= exp(1),
		Theta_CMB 	= 2.7255/2.7,
		Omega0 		= H.get<double>("Omega0"),
		h		= H.get<double>("H0") / 100.0,
		ns		= H.get<double>("ns"),
		A		= 1122670;

	return Fourier<R>::power_spectrum(
		[=] (double k)
	{
		double  q  = k * pow(Theta_CMB, 2)/(Omega0 * h),
			L0 = log(2*e + 1.8*q),
			C0 = 14.2 + 731.0/(1 + 62.5*q),
			T0 = L0 / (L0 + C0 * pow(q, 2));

		return A * pow(k, ns) * pow(T0, 2);
	});
}

template <unsigned R>
Array<double> _generate_random_field(Header const &C)
{
	unsigned    mbits = C.get<unsigned>("mbits");
	double	        L = C.get<double>("size");
	double 	    slope = C.get<double>("ns");
	bool       smooth = C.get<bool>("smooth");
	double      sigma = C.get<double>("scale");
	double       norm = C.get<double>("sigma");
	unsigned     seed = C.get<unsigned>("seed");

	size_t N = size_t(1) << mbits;

	mVector<int, R> shape(N);
	size_t size = product(shape);
	Transform fft(std::vector<int>(R, N));

	Array<double> dens(size);
	generate(dens, Gaussian_white_noise(seed));
	copy(dens, fft.in);

	if (C.get<bool>("scale-free"))
	{
		auto P = Fourier<R>::power_spectrum(
			[slope] (double k) { return pow(k, slope); });
		auto S = Fourier<R>::scale(sigma * N/L);
		auto K = kspace<R>(N, N);

		fft.forward();
		transform(fft.out, K, fft.in, Fourier<R>::filter(P * S));
		fft.in[0] = 0;
		fft.backward();
		copy(dens, fft.in);
		transform(fft.out, dens, real_part(size));

		double var = mean(map(dens, [] (double a) { return a*a; }));

		fft.forward();
		transform(fft.out, K, fft.in, Fourier<R>::filter((smooth ? P * S : P)));

		fft.in[0] = 0;
		fft.backward();
		transform(fft.out, dens, real_part(size * sqrt(var) / norm));
		
		return dens;
	}
	else
	{
		auto K = kspace<R>(N, L);
		auto P = CDM<R>(C);
		fft.forward();
		
		transform(fft.out, K, fft.in, Fourier<R>::filter(P));
		fft.in[0] = 0;
		fft.backward();
		copy(dens, fft.in);
		transform(fft.out, dens, real_part(size));

		return dens;
	}
}

template <unsigned R>
void _compute_potential(Header const &C, Array<double> density)
{
	double L = C.get<double>("size");
	unsigned mbits = C.get<unsigned>("mbits");
	unsigned N = 1 << mbits;
	size_t size = 1U << (mbits * R);

	Conan::Transform fft(std::vector<int>(R, N));
	auto K = kspace<R>(N, L);
	auto F = Fourier<R>::potential();
	copy(density, fft.in);
	fft.forward();
	transform(fft.out, K, fft.in, Fourier<R>::filter(F));
	fft.in[0] = 0;
	fft.backward();
	transform(fft.out, density, real_part(size));
}

template <unsigned R>
Array<mVector<double, R>> _compute_displacement(Header const &C, Array<double> potential)
{
	double L = C.get<double>("size");
	unsigned mbits = C.get<unsigned>("mbits");
	unsigned N = 1 << mbits;
	size_t size = 1U << (mbits * R);

	Array<complex64> phi_f(size);
	Transform fft(std::vector<int>(R, N));
	auto K = kspace<R>(N, N);
	copy(potential, fft.in);
	fft.forward();
	copy(fft.out, phi_f);

	Array<mVector<double, R>> psi(size);
	for (unsigned k = 0; k < R; ++k)
	{
		auto F = Fourier<R>::derivative(k);
		auto psi_k = access(psi, [k] (mVector<double, R> &x) -> double&
			{ return x[k]; });

		transform(phi_f, K, fft.in, Fourier<R>::filter(F));
		fft.backward();
		transform(fft.out, psi_k, real_part(size / L * N));
	}
	return psi;
}

Array<double> Conan::generate_random_field(Header const &C)
{
	unsigned      dim = C.get<unsigned>("dim");

	switch (dim)
	{
		case 1: return _generate_random_field<1>(C);
		case 2: return _generate_random_field<2>(C);
		case 3: return _generate_random_field<3>(C);
		case 4: return _generate_random_field<4>(C);
	}

	throw "only 2 and 3 dimensions supported.";
}

void Conan::compute_potential(Header const &C, Array<double> data)
{
	unsigned      dim = C.get<unsigned>("dim");

	switch (dim)
	{
		case 1: _compute_potential<1>(C, data); return;
		case 2: _compute_potential<2>(C, data); return;
		case 3: _compute_potential<3>(C, data); return;
		case 4: _compute_potential<4>(C, data); return;
	}

	throw "only 2 and 3 dimensions supported.";
}

void Conan::compute_displacement(Header const &C, Array<double> data, std::ostream &fo)
{
	unsigned dim = C.get<unsigned>("dim");

	switch (dim)
	{
		case 2: save_to_file(fo, _compute_displacement<2>(C, data), "displacement");
			return;
		case 3: save_to_file(fo, _compute_displacement<3>(C, data), "displacement");
			return;
	}
	throw "only 2 and 3 dimensions supported.";
}

