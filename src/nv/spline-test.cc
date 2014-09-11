#ifdef UNITTEST

#include "../base/unittest.hh"
#include "spline.hh"
#include "system.hh"
#include "minimizor.hh"
#include "../base/misc.hh"
#include <iostream>
#include <fstream>

using namespace System;
using namespace Conan;

inline std::function<dVector<2> ()> poisson_noise(unsigned long seed)
{
	std::shared_ptr<std::mt19937> 
		random(new std::mt19937(seed));
	std::shared_ptr<std::uniform_real_distribution<double>> 
		dist(new std::uniform_real_distribution<double>(0.0, 1.0));

	return [random, dist] () 
	{
		dVector<2> a;
		for (double &x : a)
			x = (*dist)(*random);
		return a;
	};
}

Test::Unit Spline_test("SPL",
	"gives interpolated results in 1-d and 2-d cases, plot them "
	"to see if results make sense.",
	[] ()
{
	std::ofstream fo("test.dat");

	auto noise = Gaussian_white_noise(0);

	//===============================================
	auto src_box_1 = make_ptr<Box<1>>(10, 1.0),
	     tgt_box_1 = make_ptr<Box<1>>(500, 1.0);

	Array<double> test_data_1(src_box_1->size(), 0.0);
	for (double &x : test_data_1) x = noise();

	Spline<1> S_1(src_box_1, test_data_1);

	for (unsigned i : Range<unsigned>(src_box_1->size()))
	{
		auto x = src_box_1->G[i];
		fo << x << " " << test_data_1[i] << "\n";
	}
	fo << "\n\n";

	fo << "# 1-d example\n";
	for (unsigned i : Range<unsigned>(tgt_box_1->size()))
	{
		auto x = tgt_box_1->G[i];
		fo << x << " " << S_1(x) << "\n";
	}
	fo << "\n\n";
	//===============================================

	//===============================================
	auto src_box_2 = make_ptr<Box<2>>(10, 1.0),
	     tgt_box_2 = make_ptr<Box<2>>(100, 1.0);

	Array<double> test_data_2(src_box_2->size(), 0.0);
	for (double &x : test_data_2) x = noise();

	Spline<2> S_2(src_box_2, test_data_2);

	for (unsigned i : Range<unsigned>(src_box_2->size()))
	{
		auto x = src_box_2->G[i];
		fo << x << " " << test_data_2[i] << "\n";
	}
	fo << "\n\n";

	fo << "# 2-d example\n";
	unsigned N = tgt_box_2->N();
	for (unsigned i : Range<unsigned>(tgt_box_2->size()))
	{
		auto x = tgt_box_2->G[i];
		auto y = S_2.fdf(x);

		fo << x << " " << y.first << " " << y.second << "\n";
		if (i % N == (N-1)) fo << "\n";
	}
	fo << "\n\n";
	//===============================================
	

	// test minimisation ============================
	fdfMinimizor<Spline<2>, 2> M(S_2);
	auto point_noise = poisson_noise(0);

	for (unsigned i : Range<unsigned>(50))
	{
		dVector<2> p = point_noise();
		M.set(p, src_box_2->res()/5, 0.1);
		fo << M.pos() << " " << "0" << std::endl;

		unsigned j = 0;
		while(j < 20 and M.should_continue())
		{
			M.iterate();
			++j;

			fo << M.pos() << " " << j << std::endl;
		}

		fo << std::endl;
	}
	fo.close();

	return true;
});

#endif

