#pragma once

#include <memory>

namespace System
{
	template <typename T>
	using ptr = std::shared_ptr<T>;

	template <typename T, typename ...Args>
	ptr<T> make_ptr(Args &&...args)
	{ return std::make_shared<T>(std::forward<Args>(args)...); }
}

