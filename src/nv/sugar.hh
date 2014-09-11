#include "../decomp/decomposition.hh"
#include "system.hh"

namespace Conan {
	template <unsigned R>
	inline double calc_density(BoxPtr<R> box, Array<dVector<R>> data, size_t idx);

	template <>
	inline double calc_density<2>(BoxPtr<2> box, Array<dVector<2>> data, size_t idx)
	{
		dVector<2> x[4];

		for (unsigned i = 0; i < 4; ++i)
		{
			iVector<2> q = box->I[idx] + box->block[i];
			x[i] = dVector<2>(q) * box->res() + data[box->idx(q)];
			//std::cerr << box->idx(q) << " | " << x[i] << "> ";
		}

		double rho1 = det2(x[1] - x[0], x[3] - x[1]),
		       rho2 = det2(x[2] - x[3], x[0] - x[2]);

		//std::cerr << rho1 << " " << rho2 << std::endl;

		return (rho1 + rho2) / 2.0;
	}

	template <>
	inline double calc_density<3>(BoxPtr<3> box, Array<dVector<3>> data, size_t idx)
	{
		Array<dVector<3>> x(8);

		bool parity = ((box->I[idx].sum() % 2) == 0);
		for (unsigned i = 0; i < 8; ++i)
		{
			iVector<3> q = box->I[idx] + box->block[i];
			x[i] = dVector<3>(q[i]) * box->res() + data[box->idx(q)];
		}

		return Decomposition<3>(parity, x).total_volume();
	}

	template <unsigned R>
	void sg_run(std::istream &fi, std::ostream &fo, BoxMaker const &box_, FILE_FORMAT fmt)
	{
		auto box = box_.box<R>();
		std::cerr << box->N() << " " << box->L() << " " << box->size() << std::endl;

		Array<dVector<R>> psi = System::load_from_file<dVector<R>>(fi, "eulerian-map");

		Array<double> rho(box->size());
		for (size_t i = 0; i < box->size(); ++i)
		{
			rho[i] = calc_density<R>(box, psi, i);
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

