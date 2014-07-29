#pragma once

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "system.hh"

namespace Conan
{
	namespace _Minimizor
	{
		template <typename T, int rank>
		double functorTerm(const gsl_vector *x, void *pars)
		{
			T *F = reinterpret_cast<T *>(pars);

			mVector<double, rank> A;
			for (unsigned i = 0; i < rank; ++i)
				A[i] = gsl_vector_get(x, i);

			return (*F)(A);
		}

		template <typename T, int rank>
		double f_term(const gsl_vector *x, void *pars)
		{
			T *F = reinterpret_cast<T *>(pars);

			mVector<double, rank> A;
			for (unsigned i = 0; i < rank; ++i)
				A[i] = gsl_vector_get(x, i);

			return F->f(A);
		}

		template <typename T, int rank>
		void df_term(const gsl_vector *x, void *pars, gsl_vector *g)
		{
			T *F = reinterpret_cast<T *>(pars);

			mVector<double, rank> A;
			for (unsigned i = 0; i < rank; ++i)
				A[i] = gsl_vector_get(x, i);

			for (unsigned i = 0; i < rank; ++i)
				gsl_vector_set(g, i, F->df(i, A));
		}

		template <typename T, int rank>
		void fdf_term(const gsl_vector *x, void *pars, double *v, gsl_vector *g)
		{
			T *F = reinterpret_cast<T *>(pars);

			mVector<double, rank> A;
			for (unsigned i = 0; i < rank; ++i)
				A[i] = gsl_vector_get(x, i);

			for (unsigned i = 0; i < rank; ++i)
				gsl_vector_set(g, i, F->df(i, A));

			*v = F->f(A);
		}
	}

	template <typename T, int rank>
	class fdfMinimizor
	{
		gsl_multimin_function_fdf	fn;
		gsl_multimin_fdfminimizer	*min;

		gsl_vector	*x;
		T		fctor;

		public:
			fdfMinimizor(T functor):
				fctor(functor)
			{
				fn.f = &_Minimizor::f_term<T, rank>;
				fn.df = &_Minimizor::df_term<T, rank>;
				fn.fdf = &_Minimizor::fdf_term<T, rank>;
				fn.n = rank;
				fn.params = reinterpret_cast<void *>(&fctor);

				const gsl_multimin_fdfminimizer_type *mint = gsl_multimin_fdfminimizer_vector_bfgs2;

				min = gsl_multimin_fdfminimizer_alloc(mint, rank);
				x = gsl_vector_alloc(rank);
			}

			~fdfMinimizor()
			{
				gsl_multimin_fdfminimizer_free(min);
				gsl_vector_free(x);
			}

			fdfMinimizor(fdfMinimizor const &o):
				fctor(o.fctor)
			{
				fn.f = &_Minimizor::f_term<T, rank>;
				fn.df = &_Minimizor::df_term<T, rank>;
				fn.fdf = &_Minimizor::fdf_term<T, rank>;
				fn.n = rank;
				fn.params = reinterpret_cast<void *>(&fctor);

				const gsl_multimin_fdfminimizer_type *mint = gsl_multimin_fdfminimizer_vector_bfgs2;
				min = gsl_multimin_fdfminimizer_alloc(mint, rank);
				x = gsl_vector_alloc(rank);
			}

			void set(mVector<double, rank> const &X, mVector<double, rank> const &Q)
			{
				fctor.set(X);

				for (unsigned i = 0; i < rank; ++i)
				{
					gsl_vector_set(x, i, Q[i]);
				}

				gsl_multimin_fdfminimizer_set(min, &fn, x, 0.5, 0.1);
			}

			void iterate()
			{
				gsl_multimin_fdfminimizer_iterate(min);
			}

			bool should_continue()
			{
				return gsl_multimin_test_gradient(
					min->gradient, 1e-4) != GSL_SUCCESS;
			}

			mVector<double, rank> pos()
			{
				gsl_vector *y;
				y = gsl_multimin_fdfminimizer_x(min);

				mVector<double, rank> A;
				for (unsigned i = 0; i < rank; ++i)
					A[i] = gsl_vector_get(y, i);

				return A;
			}

			double value()
			{
				return gsl_multimin_fdfminimizer_minimum(min);
			}
	};

	template <typename T, int rank>
	class Minimizor
	{
		gsl_multimin_function	fn;
		gsl_multimin_fminimizer	*min;

		gsl_vector	*step_size;
		gsl_vector	*x;

		public:
			Minimizor(T *functor)
			{
				fn.f = &_Minimizor::functorTerm<T, rank>;
				fn.n = rank;
				fn.params = reinterpret_cast<void *>(functor);

				const gsl_multimin_fminimizer_type *mint = gsl_multimin_fminimizer_nmsimplex;

				min = gsl_multimin_fminimizer_alloc(mint, rank);
				x = gsl_vector_alloc(rank);
				step_size = gsl_vector_alloc(rank);
			}

			void set(mVector<double, rank> const &A)
			{
				for (unsigned i = 0; i < rank; ++i)
				{
					gsl_vector_set(x, i, A[i]);
					gsl_vector_set(step_size, i, 1.0);
				}

				gsl_multimin_fminimizer_set(min, &fn, x, step_size);
			}

			void iterate()
			{
				gsl_multimin_fminimizer_iterate(min);
			}

			mVector<double, rank> pos()
			{
				gsl_vector *y;
				y = gsl_multimin_fminimizer_x(min);

				mVector<double, rank> A;
				for (unsigned i = 0; i < rank; ++i)
					A[i] = gsl_vector_get(y, i);

				return A;
			}

			double value()
			{
				return gsl_multimin_fminimizer_minimum(min);
			}

			double size()
			{
				return gsl_multimin_fminimizer_size(min);
			}

			bool should_continue()
			{
				return gsl_multimin_test_size(size(), 1e-4) == GSL_CONTINUE;
			}

			~Minimizor()
			{
				gsl_multimin_fminimizer_free(min);
				gsl_vector_free(x);
				gsl_vector_free(step_size);
			}		
	};
}

