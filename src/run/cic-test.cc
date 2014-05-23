#ifdef UNITTEST
#include "system.hh"
#include "numeric.hh"
#include "../base/unittest.hh"
#include <iostream>

using namespace System;
using namespace Conan;

Test::Unit cic_test("0101",
	"test cic on a random sample of points.",
	[] ()
{
	auto b = make_ptr<Box<2>>(64, 5.0); unsigned N = 256;
	
	auto noise = Gaussian_white_noise(0);
	Array<dVector<2>> X(N);
	for (unsigned i = 0; i < N; ++i)
	{
		X[i][0] = noise() + 2.5; X[i][1] = noise() + 2.5;

		std::cout << X[i] << std::endl;
	}
	std::cout << "\n\n";

	Array<double> A(b->size());
	md_cic(b, A, X);

	size_t o = 0;
	for (unsigned i = 0; i < b->N(); ++i)
	{
		for (unsigned j = 0; j < b->N(); ++j, ++o)
		{
			std::cout << A[o] << " ";
		}
		std::cout << std::endl;
	}

	return true;
});

#endif

