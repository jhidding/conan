#pragma once

#include "../base/system.hh"
#include <functional>

namespace Conan
{
	template <unsigned R>
	using BoxPtr = System::ptr<System::Box<R>>;

	template <unsigned R>
	using dVector = System::mVector<double, R>;

	template <unsigned R>
	using iVector = System::mVector<int, R>;

	using System::Array;
	using System::ptr;
	using System::BoxMaker;

	enum FILE_FORMAT { FMT_ASCII, FMT_IFRIT, FMT_CONAN };

	template <unsigned R>
	using Function = std::function<double (dVector<R> const &)>;

	template <unsigned R>
	inline iVector<R> floor_cast(dVector<R> const &x)
	{
		iVector<R> y;
		for (unsigned i = 0; i < R; ++i) 
			y[i] = static_cast<int>(x[i] < 0 ? x[i] - 1.0 : x[i]);
		return y;
	}

	template <unsigned R>
	inline iVector<R> ceil_cast(dVector<R> const &x)
	{
		iVector<R> y;
		for (unsigned i = 0; i < R; ++i) 
			y[i] = static_cast<int>(x[i] < 0 ? x[i] : x[i] + 1.0);
		return y;
	}
}

