#pragma once
#include <utility>
#include <memory>
#include <functional>
#include <algorithm>

#include "../base/mvector.hh"

/*
 * KdTree algorithm
 * To count the number of points that satisfy a certain
 * predicate, you have to define a function object
 * that derives from predicate.
 * This object needs not be constant. So that a comparison
 * with a bounding box can influence the way points within 
 * this box are tested.
 */

namespace kdTree {

using System::mVector;

// should change to some different template scheme

template <typename Point, unsigned R>
class BoundingBox
{
	std::function<double (Point const &, int)> index;

	public:
		mVector<double, R> X1, X2;
		BoundingBox(std::function<double (Point const &, int)> index_):
       			index(index_) {}

		BoundingBox(
			mVector<double, R> const &X1_, 
			mVector<double, R> const &X2_,
			std::function<double (Point const &, int)> index_);

		template <typename Iter>
		BoundingBox(Iter begin, Iter end,
			std::function<double (Point const &, int)> index_);

		double max_coord(unsigned i) const { return X2[i]; }
		double min_coord(unsigned i) const { return X1[i]; }
		double half(unsigned i);
};

template <typename Point, unsigned R>
class Predicate
{
	public:
		virtual bool operator()(Point const &) const = 0;
		virtual bool operator()(BoundingBox<Point, R> const &) const = 0;
		virtual ~Predicate() {}
};

template <typename Point>
class Visitor
{
	public:
		virtual void operator()(Point &) = 0;
		virtual ~Visitor() {}
};

template <typename Point, unsigned R>
class Branch: public BoundingBox<Point, R>
{
	public:
		Branch() {}

		template <typename Iter>
		Branch(Iter begin, Iter end,
			std::function<double (Point const &, int)> idx): 
			BoundingBox<Point, R>(begin, end, idx) {}

		virtual size_t count_if(Predicate<Point, R> const &pred) const = 0;
		virtual ~Branch() {}
		virtual void traverse(Visitor<Point> &visit, Predicate<Point, R> const &pred) = 0;
};

template <typename Point, unsigned R>
class Tree: public Branch<Point, R>
{
	typedef std::unique_ptr<Branch<Point, R>> Branch_ptr;
	std::pair<Branch_ptr, Branch_ptr> B;

	public:
		template <typename Iter>
		Tree(Iter begin, Iter end, 
			std::function<double (Point const &, int)> idx, 
			int dim = 0, int depth = 0);

		virtual size_t count_if(Predicate<Point, R> const &pred) const;
		virtual void traverse(Visitor<Point> &visit, Predicate<Point, R> const &pred);

//	private:
//		Tree(Tree const &); //NI
};

template <typename Point, unsigned R, typename Iter>
class Leaf: public Branch<Point, R>
{
	std::pair<Iter, Iter>	I;

	public:
		Leaf(Iter begin, Iter end,
			std::function<double (Point const &, int)> idx);

		virtual size_t count_if(Predicate<Point, R> const &pred) const;
		virtual void traverse(Visitor<Point> &visit, Predicate<Point, R> const &pred);
};

} // namespace kdTree

#include "kdtree_implementation.hh"

