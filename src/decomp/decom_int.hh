#pragma once
#include "system.hh"
#include "numeric.hh"

namespace Conan
{
	template <typename T, unsigned R>
	int arg_min(mVector<T, R> const &v);

	template <typename T>
	int arg_min(mVector<T, 2> const &v)
	{ return (v[0] < v[1] ? 0 : 1); }

	template <typename T>
	int arg_min(mVector<T, 3> const &v)
	{ return (v[0] < v[1] ? (v[0] < v[2] ? 0 : 2) : (v[1] < v[2] ? 1 : 2)); }

	inline int fsign(double a)
	{
		if (a > 0.0) return 1;
		if (a < 0.0) return -1;
		return 0;
	}

	template <typename Func>
	void Decomposition<2>::for_each_intersection(Func blub) const
	{
		iVector<2> dx[2];
		for (unsigned k = 0; k < 2; ++k)
			dx[k][k] = 1;

		// for each cell
		for (uint8_t c_ = 0;  c_ < 2; ++c_)
		{
			double density = D.C[c_].unit_volume / d_cell_volume[c_];
			Decom_data_2::Cell const &c = D.C[c_];
			std::cerr << (d[c.vertices[0]] + d[c.vertices[1]] + d[c.vertices[2]])/3. << std::endl;

			for_each_decrement(c.vertices, [&]
				(Vec const &A, uint8_t a)
			{
				for_each_decrement(A, [&]
					(Vec const &B, uint8_t b)
				{
					mVector<double, 2> 
						c = d[B[0]],
						l = d[b] - d[a],
						p;

					l /= l.norm();
					p[0] = l[1]; p[1] = -l[0];
					p *= (p.dot(c - d[a]) > 0 ? 1 : -1);
					blub(d[a].floor(), density * d[a].dot(l) * d[a].dot(p) / 2.);
				});

				mVector<double, 2>
					l = d[A[1]] - d[A[0]],
					p;

				l /= l.norm();
				p[0] = l[1]; p[1] = -l[0];
				// if sign is positive, triangle is on the right of the line A[0]->A[1]
				int sign = (p.dot(d[a] - d[A[0]]) > 0 ? 1 : -1);

				// p points into the triangle
				p *= sign;

				for_all_line_intersections(d[A[0]], d[A[1]],
					[&] (mVector<double, 2> const &x, unsigned k)
				{
					unsigned m = (k + 1) % 2;
					mVector<int, 2> X; X[k] = round(x[k]); X[m] = floor(x[m]);
					//std::cerr << x << std::endl << X << "\n\n\n";
					//std::cerr << x << " " << l << " " << 1 << std::endl;
					//std::cerr << x << " " << p << " " << 2 << std::endl;

					blub(X, density * x[k] * x[m] * fsign(p[m]) / 2.0);
					blub(X - dx[k], - density * x[k] * x[m] * fsign(p[m]) / 2.0);
					blub(X, density * x.dot(l) * x.dot(p) * fsign(l[k]) / 2.0);
					blub(X - dx[k], - density * x.dot(l) * x.dot(p) * fsign(l[k]) / 2.0);
				});
			});
		}
	}

	template <typename Func>
	void Decomposition<2>::loop_vefas_for_vertex(uint8_t v_, Func blub) const
	{
		Decom_data_2::Vertex const &v = D.V[v_];
		// for all incident edges
		for (auto e_ = v.begin(); e_ != v.end(); ++e_)
		{
			// e_ points to pair(target vertex idx, resulting edge)
			Decom_data_2::Edge const &e = D.E[e_->second];
			mVector<double, 2> 
				l = d[e_->first] - d[v_],
				p;

			l /= l.norm();
			// p perpendicular to l -> quarter turn
			p[0] = l[1];
			p[1] = -l[0];

			double z = d[v_].dot(l) * d[v_].dot(p);
			// both adjacent faces/cells (nomenclature compromise?)
			for (auto c_ = e.begin(); c_ != e.end(); ++c_)
			{
				Decom_data_2::Cell const &c = D.C[c_->second];
				double density = c.unit_volume / d_cell_volume[c_->second];
				blub(density * z * (p.dot(d[c_->first] - d[v_]) > 0 ? 1 : -1) / 2.0);
			}
		}
	}

	template <typename Func>
	void Decomposition<3>::loop_vefas_for_vertex(uint8_t v_, Func blub) const
	{
		Decom_data_3::Vertex const &v = D.V[v_];
		for (auto e_ = v.begin(); e_ != v.end(); ++e_)
		{
			Decom_data_3::Edge const &e = D.E[e_->second];
			mVector<double, 3> 
				l = d[e_->first] - d[v_];

			for (auto f_ = e.begin(); f_ != e.end(); ++f_)
			{
				Decom_data_3::Face const &f = D.F[f_->second];
				mVector<double, 3>
					n = d_normals[f_->second],
					p = cross(n, l);

				p *= (p.dot(d[f_->first] - d[v_]) > 0 ? 1 : -1);

				for (auto c_ = f.begin(); c_ != f.end(); ++c_)
				{
					Decom_data_3::Cell const &c = D.C[c_->second];
					double density = c.unit_volume / d_cell_volume[c_->second];
					n *= (n.dot(d[c_->first] - d[v_]) > 0 ? 1 : -1);

					blub(density * d[v_].dot(l) * d[v_].dot(p) * d[v_].dot(n));
				}
			}
		}
	}

	template <typename Func>
	void Decomposition<3>::for_each_intersection(Func blub) const
	{
		mVector<int, 3> dx[3];
		for (unsigned k = 0; k < 3; ++k)
			dx[k][k] = 1;

		// several problems force me to loop first over the cells
		// for all segments
		for (uint8_t c_ = 0; c_ < 5; ++c_)
		{
			double density = D.C[c_].unit_volume / d_cell_volume[c_];

			// for all faces
			// for each cell
			for_each_decrement(D.C[c_].vertices,
				[&] (Vec const &A, uint8_t a)
			{
				// add vertex [a], edge [ab], face [abc], adjacency
				// one time for each permutation of [abcd]
				for_each_decrement(A,
					[&] (Vec const &B, uint8_t b)
				{
					Point l = d[b] - d[a];
					l /= l.norm();

					for_each_decrement(B,
						[&] (Vec const &C, uint8_t c)
					{
						uint8_t dd = C[0]; 
						
						Point n = cross(d[b] - d[a], d[c] - d[a]);
						n /= n.norm();

						Point q = cross(l, n);
						q *= (q.dot(d[c] - d[a]) > 0 ? 1 : -1); // pointing into face
						n *= (n.dot(d[dd] - d[a]) > 0 ? 1 : -1); // pointing into cell

						blub(d[a].floor(), - density * d[a].dot(l) * d[a].dot(q) * d[a].dot(n) / 6.);
					});
				});

				// determine normal of face;
				Point n = cross(d[A[1]] - d[A[0]], d[A[2]] - d[A[0]]);
				n /= n.norm();
				n *= (n.dot(d[a] - d[A[0]]) > 0 ? 1 : -1); // pointing into cell

				std::vector<Point> V(3);
				V[0] = d[A[0]]; V[1] = d[A[1]]; V[2] = d[A[2]];
				for_all_ef_intersections(V,
					[&] (Point const &x, Point const &l, Point const &b, int k) // line intersection
				{
					int i = (k + 1) % 3, j = (k + 2) % 3;
					mVector<int, 3> X; X[k] = round_down(x[k] + 0.5); 
					X[i] = round_down(x[i]); X[j] = round_down(x[j]);

					Point u = cross(n, l);
					u *= (Conan::cross(l, b - x).dot(n) > 0 ? 1 : -1); // pointing into face

					double z = -density * x.dot(u) * x.dot(l) * x.dot(n) * fsign(l[k]) / 6.;
					blub(X, z);
					blub(X - dx[k], -z);

					Point q = calc_q(n, k, i, j), p = cross(n, q);
					q *= (q.dot(u) > 0 ? 1 : -1);
					z = -density * x.dot(q) * x.dot(p) * x.dot(n) * fsign(p[k]) / 6.;
					blub(X, z);
					blub(X - dx[k], -z);

					p = cross(q, Point(dx[k]));
					p *= (p.dot(n) > 0 ? 1 : -1);
					z = -density * x.dot(q) * x.dot(p) * x[k] / 6.; 
					blub(X, z);
					blub(X - dx[k], -z);
				},
					[&] (Point const &x, int k) // face intersection
				{
					int i = (k + 1) % 3, j = (k + 2) % 3;
					mVector<int, 3> X; X[k] = round_down(x[k]); 
					X[i] = round_down(x[i] + 0.5); X[j] = round_down(x[j] + 0.5);
					Point p, q;

					q = calc_q(n, i, j, k);
					q *= fsign(q[j]);
					p = cross(n, q);
					p *= fsign(p[i]);

					double z = - density * x.dot(q) * x.dot(p) * x.dot(n) / 6.;
					blub(X, z);
					blub(X - dx[i], -z);
					blub(X - dx[j], -z);
					blub(X - dx[j] - dx[i], z);

					p = cross(q, Point(dx[i]));
					p *= (p.dot(n) > 0 ? 1 : -1);
					z = - density * x.dot(q) * x.dot(p) * x[i] / 6.;
					blub(X, z);
					blub(X - dx[i], -z);
					blub(X - dx[j], -z);
					blub(X - dx[j] - dx[i], z);

					q = calc_q(n, j, i, k);
					q *= fsign(q[i]);
					p = cross(n, q);
					p *= fsign(p[j]);

					double w = - density * x.dot(q) * x.dot(p) * x.dot(n) / 6.;
					blub(X, w);
					blub(X - dx[j], -w);
					blub(X - dx[i], -w);
					blub(X - dx[i] - dx[j], w);

					p = cross(q, Point(dx[j]));
					p *= (p.dot(n) > 0 ? 1 : -1);
					w = - density * x.dot(q) * x.dot(p) * x[j] / 6.;
					blub(X, w);
					blub(X - dx[j], -w);
					blub(X - dx[i], -w);
					blub(X - dx[i] - dx[j], w);

					w = - density * x[k] * fsign(n[k]) * x[i] * x[j] / 3.;
					blub(X, w);
					blub(X - dx[j], -w);
					blub(X - dx[i], -w);
					blub(X - dx[i] - dx[j], w);
				});
			});
		}
	}
}
