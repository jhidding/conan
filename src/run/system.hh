#pragma once

#include "../base/system.hh"

namespace Conan
{
	using System::FILE_FORMAT;
	using System::FMT_ASCII;
	using System::FMT_CONAN;
	using System::FMT_IFRIT;

	template <unsigned R>
	using BoxPtr = System::ptr<System::Box<R>>;

	template <unsigned R>
	using dVector = System::mVector<double, R>;

	template <unsigned R>
	using iVector = System::mVector<int, R>;

	using System::Array;
	using System::ptr;
}

