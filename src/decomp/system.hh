#pragma once

#include "../base/system.hh"

namespace Conan
{
	using System::mVector;
	using System::round_up;
	using System::round_down;

	template <unsigned R>
	using BoxPtr = System::ptr<System::Box<R>>;

	template <unsigned R>
	using dVector = mVector<double, R>;

	template <unsigned R>
	using iVector = mVector<int, R>;

	using System::Array;
	using System::ptr;
}

