#include "../decomp/decomposition.hh"
#include "system.hh"

namespace Conan {
	template <unsigned R>
	void sg_run(std::istream &fi, std::ostream &fo, BoxMaker const &box_, FILE_FORMAT fmt)
	{
		auto box = box_.box<R>();
		std::cerr << box->N() << " " << box->L() << " " << box->size() << std::endl;

		Array<dVector<R>> psi = System::load_from_file<dVector<R>>(fi, "eulerian-map");

		Array<double> rho(box->size());
		for (size_t i = 0; i < box->size(); ++i)
		{
			iVector<R> x = box->G[i];
			bool parity = (x.sum() % 2 == 0);

			Array<dVector<R>> d(2 << R);
			for (unsigned j = 0; j < (2 << R); ++j)
			{
				iVector<R> y = x + box->block[j];
				d[j] = psi[box->idx(y)];
			}

			rho[i] = 1./Decomposition<R>(parity, d).total_volume();
		}

		switch (fmt)
		{
			case FMT_ASCII:
				write_matrix_array_txt(fo, box, rho); break;

			case FMT_CONAN:
				System::save_to_file(fo, rho, "density"); break;
		}
	}
}

