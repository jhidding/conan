#pragma once
#include "ptr.hh"
#include <iostream>

namespace System
{
	template <typename T>
	class Maybe
	{
		ptr<T> X;

		public:
			Maybe() {}
			Maybe(T const &X_): X(new T(X_)) {}
			Maybe(ptr<T> X_): X(X_) {}

			operator bool() const { return static_cast<bool>(X); }
			T const &operator*() const { return *X; }
			T const *operator->() const { return X.get(); }
			operator ptr<T>() const { return X; }
	};

	class MaybeNothing
	{
		public:
			template <typename T>
			operator Maybe<T>() const { return Maybe<T>(); }
	};
	
	extern MaybeNothing Nothing;

	template <typename T>
	Maybe<T> Just(T const &X) { return Maybe<T>(X); }

	template <typename T, typename ...Args>
	Maybe<T> Just(Args &&...args)
	{
		return Maybe<T>(make_ptr<T>(std::forward<Args>(args)...)); 
	}

	template <typename T>
	Maybe<T> Just(ptr<T> X) { return Maybe<T>(X); }

	template <typename T>
	std::ostream &operator<<(std::ostream &out, Maybe<T> const &X)
	{
		if (X)
			return out << *X;
		else
			return out << "Nothing";
	}
}

