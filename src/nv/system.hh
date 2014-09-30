#pragma once

#include "../base/system.hh"
#include <iostream>
#include <iomanip>
#include <functional>

namespace Conan
{
	using System::FILE_FORMAT;
	using System::FMT_ASCII;
	using System::FMT_CONAN;
	using System::FMT_IFRIT;

	using System::Maybe;
	using System::Just;
	using System::Nothing;

	template <unsigned R>
	using BoxPtr = System::ptr<System::Box<R>>;

	template <unsigned R>
	using dVector = System::mVector<double, R>;

	template <unsigned R>
	using iVector = System::mVector<int, R>;

	using System::Array;
	using System::ptr;
	using System::make_ptr;
	using System::BoxMaker;
	using System::MdRange;

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

	template <typename T, unsigned R>
	void write_matrix_array_txt(std::ostream &fo, BoxPtr<R> box, Array<T> data)
	{
		for (size_t i = 0; i < box->size(); ++i)
		{
			fo << std::setw(5) << data[i] << " ";
			
			for (unsigned k = 0; k < R-1; ++k)
			{
				if (box->I[i][k] == (box->shape()[k] - 1))
				{	
					fo << std::endl;
				}
			}
		}
	}

	template <typename T, unsigned R>
	void write_list_array_txt(std::ostream &fo, BoxPtr<R> box, Array<T> data)
	{
		for (size_t i = 0; i < box->size(); ++i)
		{
			fo << std::setw(5) << box->G[i] << " " << data[i] << "\n";
		}
	}
}

