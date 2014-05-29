#include "decomposition.hh"
#include "decom_data.hh"

using namespace Conan;

Decom_data_2::Cell_structure _p2_odd[2] = {
	{ 1./2, { 0, 2, 1 } }, 
	{ 1./2, { 3, 1, 2 } } };

Decom_data_2::Cell_structure _p2_even[2] = {
	{ 1./2, { 1, 0, 3 } }, 
	{ 1./2, { 2, 3, 0 } } };


Decom_data_3::Cell_structure _p3_odd[5] = {
	{ 1./6, {7, 5, 6, 3} }, 
	{ 1./6, {4, 6, 5, 0} }, 
	{ 1./6, {1, 5, 3, 0} }, 
	{ 1./6, {2, 6, 0, 3} },
	{ 1./3, {0, 5, 3, 6} } };

Decom_data_3::Cell_structure _p3_even[5] = {
	{ 1./6, {0, 1, 2, 4} }, 
	{ 1./6, {3, 2, 1, 7} }, 
	{ 1./6, {5, 1, 4, 7} }, 
	{ 1./6, {6, 2, 7, 4} }, 
	{ 1./3, {7, 1, 4, 2} } };


Decom_data_2 const Decomposition<2>::dd_odd(_p2_odd);
Decom_data_2 const Decomposition<2>::dd_even(_p2_even);
Decom_data_3 const Decomposition<3>::dd_odd(_p3_odd);
Decom_data_3 const Decomposition<3>::dd_even(_p3_even);

dVector<3> Decomposition<3>::face_normal(Vec const &v) const
{
	dVector<3> c = cross(d[v[1]] - d[v[0]], d[v[2]] - d[v[0]]);
	return c / c.norm();
}

double Decomposition<2>::calc_volume(Vec const &v) const
{
	return 1./2 * det2(d[v[1]] - d[v[0]], d[v[2]] - d[v[0]]);
}

//double Decomposition<2>::point2face_volume(mVector<double, 2> const &p, Vec const &v) const
//{
//	return 1./2 * det2(d[v[0]] - p, d[v[1]] - p);
//}

double Decomposition<3>::calc_volume(Vec const &v) const
{
	return 1./6 * det3(d[v[1]] - d[v[0]], d[v[2]] - d[v[0]], d[v[3]] - d[v[0]]);
}

//double Decomposition<3>::point2face_volume(mVector<double, 3> const &p, Vec const &v) const
//{
//	return 1./6 * det3(d[v[0]] - p, d[v[1]] - p, d[v[2]] - p);
//}

bool Decomposition<2>::cell_contains(uint8_t c, dVector<2> const &p) const
{
	bool q = true;
	for_each_decrement(D.C[c].vertices,
		[&] (Vec const &A, uint8_t b)
	{
		double u = det2(d[A[0]] - d[b], d[A[1]] - d[b]) * det2(d[A[0]] - p, d[A[1]] - p);
		if (u < 0) q = false;
	});
	return q;
}

bool Decomposition<3>::cell_contains(uint8_t c, dVector<3> const &p) const
{
	bool q = true;
	for_each_decrement(D.C[c].vertices,
		[&] (Vec const &A, uint8_t b)
	{
		double u = det3(d[A[0]] - d[b], d[A[1]] - d[b], d[A[2]] - d[b]) 
			* det3(d[A[0]] - p, d[A[1]] - p, d[A[2]] - p);
		if (u < 0) q = false;
	});
	return q;
}


