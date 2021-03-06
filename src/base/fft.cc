#include "fft.hh"
#include <vector>
#include <omp.h>

using namespace Fourier;

template <typename T>
inline typename T::value_type _product(T a)
{
	typename T::value_type v = 1;
	for (auto x : a) v *= x;
	return v;
}

Transform::Transform(std::vector<int> const &shape):
	size(_product(shape)), in(size), out(size)
{
	std::vector<int> ishape(shape.begin(), shape.end());
	fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());

	d_plan_fwd = fftw_plan_dft(ishape.size(), ishape.data(),
		reinterpret_cast<fftw_complex *>(in.data()), 
		reinterpret_cast<fftw_complex *>(out.data()), FFTW_FORWARD, FFTW_ESTIMATE);

	d_plan_bwd = fftw_plan_dft(ishape.size(), ishape.data(),
		reinterpret_cast<fftw_complex *>(in.data()),
		reinterpret_cast<fftw_complex *>(out.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
}

Transform::~Transform()
{
	fftw_destroy_plan(d_plan_fwd);
	fftw_destroy_plan(d_plan_bwd);
}

void Transform::forward() { fftw_execute(d_plan_fwd); }
void Transform::backward() { fftw_execute(d_plan_bwd); }

