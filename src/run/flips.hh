#pragma once
#include "system.hh"

namespace Conan
{
	template <unsigned R>
	class Dimension_traits
	{};

	template <>
	class Dimension_traits<2>
	{
		public:
		enum { R = 2, block_size = 4, n_cells = 2 };
	};

	template <>
	class Dimension_traits<3>
	{
		public:
		enum { R = 3, block_size = 8, n_cells = 5 };
	};

	template <unsigned R>
	void find_flip_flops(BoxPtr<R> box, Array<dVector<R>>, Array<unsigned> ff)
	{
		typedef Dimension_traits<R> traits;

		#pragma omp parallel for
		for (size_t i = 0; i < box->size(); ++i)
		{
			iVector<R> x = box->I[i];
			bool idx_parity = (sum(x) % 2 == 0),
			     ff_even    = (ff[i] % 2 == 0);

			Array<dVector<R>> pts(traits::block_size);

			// copy a 4-neighbourhood of points, correct for periodic boundaries
			
			double v = Decomposition<R>(idx_parity, pts).total_volume();

			if ((v < 0.0 and ff_even) or (v > 0.0 and not ff_even))
				++(ff[i]);
		}
	}
}
