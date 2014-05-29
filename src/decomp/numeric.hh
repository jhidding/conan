#pragma once

#include "system.hh"
#include <vector>

namespace Conan
{
	template <typename T, typename Func>
	inline void for_each_decrement(
		std::vector<T> const &A,
		Func f)
	{
		for (uint8_t i = 0; i < A.size(); ++i)
		{
			std::vector<T> B(A.size() - 1);

			uint8_t k = 0;
			for (uint8_t j = 0; j < A.size(); ++j)
			{
				if (j == i) continue;
				B[k] = A[j];
				++k;
			}

			f(B, A[i]);
		}
	}

	template <typename T, unsigned R>
	mVector<T,R> cross(mVector<T,R> const &a, mVector<T,R> const &b)
	{
		return mVector<T,R>({ a[1]*b[2] - a[2]*b[1],
				      a[2]*b[0] - a[0]*b[2],
				      a[0]*b[1] - a[1]*b[0] });
	}

	template <typename T>
	inline double det2(
		mVector<T, 2> const &a,
		mVector<T, 2> const &b)
	{
		return (a[0] * b[1] - b[0] * a[1]);
	}

	template <typename T>
	inline double det3(
		mVector<T, 3> const &a,
		mVector<T, 3> const &b,
		mVector<T, 3> const &c)
	{
		return (  a[0] * (b[1] * c[2] - c[1] * b[2])
			- b[0] * (a[1] * c[2] - c[1] * a[2])
			+ c[0] * (a[1] * b[2] - b[1] * a[2]));
	}


	template <typename T>
	inline mVector<T, 3> calc_q(mVector<T, 3> const &n, int i, int j, int k)
	{
		mVector<T, 3> q;
		// projection along axis k, vector lying in jk-plane
		// q[i] = 0; n . q = 0; q^2 = 1.
		q[i] = 0;
		q[j] = sqrt(1. / (1 + pow(n[j] / n[k], 2)));
		q[k] = - n[j] * q[j] / n[k];
		return q;
	}

	template <typename T>
	inline mVector<T, 3> calc_p(mVector<T, 3> const &n, int i, int j, int k)
	{
		mVector<T, 3> p;
		// projection along axis k, vector perpendicular to q
		// n . p = 0; p . q = 0; p^2 = 1.
		p[i] = sqrt(pow(n[j], 2) + pow(n[k], 2));
		p[j] = - n[j] * n[i] / p[i];
		p[k] = - n[k] * n[i] / p[i];
		return p;
	}

	// give doubles a, b and function to be called on each integer in (a, b)
	// and "time" t where a lies at t = 0 and b at t = 1.
	template <typename Func>
	void for_all_integers_in_interval(double a, double b, Func f)
	{
		double v = b - a;
		a = (v < 0.0 ? round_down(a) : round_up(a));

		for (double t = 1.0 - (b - a) / v; t < 1.0; t += fabs(1./v))
			f(t);
	}

	template <typename T, unsigned R, typename Func>
	void for_all_line_intersections(mVector<T, R> const &a, mVector<T, R> const &b,
		Func f)
	{
		mVector<T, R> v = b - a;
		for (unsigned k = 0; k < R; ++k)
		{
			for_all_integers_in_interval(a[k], b[k],
				[&] (double t)
			{
				f(a + v * t, k);
			});
		}
	}

	/*           p1
	 *           *,,
	 *          /|  ''',,
	 *         / |     ,,;* p2
	 *        /  | ,,''
	 *       / ,,'' 
	 *      *''     k-axis --> 
	 *      p0
	 * t = |---------------|
	 *     0     t_*       1
	 */

	template <typename T, int R, typename Func1, typename Func2>
	void for_all_ef_intersections(std::vector<mVector<T, R>> const &P,
		Func1 on_edge, Func2 on_face)
	{
		for (unsigned m = 0; m < R; ++m)
		{
			// project down to dimension m, 
			// call other directions k and l
			unsigned
				k = (m + 1) % R,
				l = (m + 2) % R;

			// sort points on axis k
			int pmt[3] = {0, 1, 2};
			std::sort(pmt, pmt + 3, 
				[P, k] (int i, int j)
			{
				return P[i][k] < P[j][k];
			});

			// calculate relative vectors, might be most expensive part
			// possible room for optimisation...
			mVector<T, R> 
				v_02 = P[pmt[2]] - P[pmt[0]], vu_02 = v_02 / v_02.norm(),
				v_01 = P[pmt[1]] - P[pmt[0]], vu_01 = v_01 / v_01.norm(),
				v_12 = P[pmt[2]] - P[pmt[1]], vu_12 = v_12 / v_12.norm();

			// calculate t_*
			double t_star = v_01[k] / v_02[k];

			for_all_integers_in_interval(P[pmt[0]][k], P[pmt[2]][k],
				[&] (double t)
			{
				mVector<T, R>
					c1 = P[pmt[0]] + v_02 * t,
					c2 = (t < t_star ?
						P[pmt[0]] + v_01 * (t / t_star) :
						P[pmt[1]] + v_12 * (t - t_star) / (1.0 - t_star)),
					dc = c2 - c1;

				on_edge(c1, vu_02, P[pmt[1]], k);
				on_edge(c2, (t < t_star ? vu_01 : vu_12), 
					    (t < t_star ? P[pmt[2]] : P[pmt[0]]), k);

				for_all_integers_in_interval(c1[l], c2[l],
					[&] (double u)
				{
					on_face(c1 + dc * u, m);
				});
			});
		}
	}
}

