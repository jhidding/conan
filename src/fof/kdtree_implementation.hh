//#include "kdtree.hh"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iterator>
#include <iostream>

namespace kdTree {

template <typename Point>
inline std::function<bool (Point const &, Point const &)> compare_element(
	unsigned i,
	std::function<double (Point const &, int)> &index)
{
	return [i, &index] (Point const &A, Point const &B) -> bool 
		{ return index(A, i) < index(B, i); };
}

/*
 * BoundingBox ---------------------------------------------------
 */

template <typename Point, unsigned R>
BoundingBox<Point, R>::BoundingBox(
	mVector<double, R> const &X1_, 
	mVector<double, R> const &X2_,
	std::function<double (Point const &, int)> index_)
	: X1(X1_), X2(X2_), index(index_) {}

template <typename Point, unsigned R>
template <typename Iter>
BoundingBox<Point, R>::BoundingBox(Iter begin, Iter end,
	std::function<double (Point const &, int)> index_):
	index(index_)
{
	for (unsigned i = 0; i < R; ++i)
	{
		X1[i] = index(*std::min_element(begin, end, compare_element<Point>(i, index)), i);
		X2[i] = index(*std::max_element(begin, end, compare_element<Point>(i, index)), i);
	}

	std::cout << X1 << "\n" << X1[0] << " " << X2[1] << "\n"
	  << X2 << "\n" << X2[0] << " " << X1[1] << "\n" << X1 << "\n\n";
}

template <typename Point, unsigned R>
double BoundingBox<Point, R>::half(unsigned i)
{
	return (X2[i] + X1[i]) / 2;
}

/*
 * Tree implementation -------------------------------------------------
 */

template <typename Point, unsigned R>
template <typename Iter>
Tree<Point, R>::Tree(Iter begin, Iter end, 
	std::function<double (Point const &, int)> idx, 
	int dim, int depth)
	: Branch<Point, R>(begin, end, idx)
{
	double boundary = BoundingBox<Point, R>::half(dim);
	Iter mid = std::partition(begin, end,
		[&] (Point const &p)
	{
		return idx(p, dim) < boundary;
	});
	
	if ((end - mid) < 32 or depth > 128)
		B.second = Branch_ptr(dynamic_cast<Branch<Point, R> *>(
			new Leaf<Point, R, Iter>(mid, end, idx)));
	else
		B.second = Branch_ptr(dynamic_cast<Branch<Point, R> *>(
			new Tree<Point, R>(mid, end, idx, (dim + 1) % R, depth+1)));

	if ((mid - begin) < 32 or depth > 128)
		B.first = Branch_ptr(dynamic_cast<Branch<Point, R> *>(
			new Leaf<Point, R, Iter>(begin, mid, idx)));
	else
		B.first = Branch_ptr(dynamic_cast<Branch<Point, R> *>(
			new Tree<Point, R>(begin, mid, idx, (dim + 1) % R, depth+1)));
}


template <typename Point, unsigned R>
size_t Tree<Point, R>::count_if(Predicate<Point, R> const &pred) const
{
	if (pred(*this))
		return  B.first->count_if(pred) +
			B.second->count_if(pred);
	else
		return 0;
}

template <typename Point, unsigned R>
void Tree<Point, R>::traverse(Visitor<Point> &visit, Predicate<Point, R> const &pred)
{
	if (not pred(*this)) return;

	B.first->traverse(visit, pred);
	B.second->traverse(visit, pred);
}

/*
 * Leaf implementation ---------------------------------------------
 */

template <typename Point, typename F>
std::function<bool (Point const &)> Pack(F const &pred)
{
	return [&pred] (Point const &A) -> bool { return pred(A); };
}

template <typename Point, unsigned R, typename Iter>
Leaf<Point, R, Iter>::Leaf(Iter begin, Iter end,
	std::function<double (Point const &, int)> idx)
	: Branch<Point, R>(begin, end, idx), I(begin, end)
{}

template <typename Point, unsigned R, typename Iter>
size_t Leaf<Point, R, Iter>::count_if(Predicate<Point, R> const &pred) const
{
	if (pred(*this))
	{
		return std::count_if(I.first, I.second, Pack<Point>(pred));
	}
	else
	{
		return 0;
	}
}

template <typename Point, unsigned R, typename Iter>
void Leaf<Point, R, Iter>::traverse(Visitor<Point> &visit, Predicate<Point, R> const &pred)
{
	if (not pred(*this)) return;
	Iter pos = I.first;
	while ((pos = std::find_if(pos, I.second, Pack<Point>(pred))) != I.second)
	{
		visit(*pos);
		++pos;
	}
}

} // namespace KdTree

