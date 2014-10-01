#include "../base/system.hh"

using namespace System;

#include <stack>
#include "kdtree.hh"
#include "../base/progress.hh"

template <unsigned R>
using BoxPtr = System::ptr<System::Box<R>>;

template <unsigned R>
using dVector = System::mVector<double, R>;

using System::FILE_FORMAT;

template <typename T>
void flood_fill(
	std::function<Array<T> (T)> neighbours,
	std::function<bool (T)>     predicate,
	std::function<void (T)>     pass,
	T start)
{
	std::stack<T> S;
	pass(start);
	S.push(start);

	while (not S.empty())
	{
		T a = S.top();
		S.pop();

		for (auto b : neighbours(a))
			if (predicate(b))
			{
				pass(b);
				S.push(b);
			}
	}
}

class Collector: public kdTree::Visitor<size_t>
{
	Array<size_t>	A;

	public:
		Collector(): A(0) {}

		virtual void operator()(size_t &i)
		{
			A->push_back(i);
		}

		Array<size_t> yield() const
		{
			return A;
		}
};

template <unsigned R>
class PrintVector: public kdTree::Visitor<size_t>
{
	Array<dVector<R>> p;
	size_t i;

	public:
		PrintVector(Array<dVector<R>> pos_, size_t i_):
			p(pos_), i(i_) {}

		virtual void operator()(size_t &j)
		{
			std::cout << p[i] << " " << p[j] - p[i] << " " << (p[i] - p[j]).norm() << std::endl;
		}
};

inline double fmodulus(double a, double b)
{
	if (a < 0) return a - static_cast<int>(a / b - 1) * b;
	else return a - static_cast<int>(a / b) * b;
}

template <unsigned R>
class Proximity: public kdTree::Predicate<size_t, R>
{
	BoxPtr<R>		box;
	Array<dVector<R>> 	pos;
	dVector<R>		a;
	double			r, rsqr;

	public:
		Proximity(BoxPtr<R> box_, Array<dVector<R>> pos_, dVector<R> a_, double r_):
			box(box_), pos(pos_), a(a_), r(r_), rsqr(r_*r_) {}

		virtual bool operator()(size_t const &i) const
		{
			dVector<R>	b = pos[i] - a, d;
			double L = box->L();

			for (unsigned k = 0; k < R; ++k)
			{
				if (b[k] > L/2)
				{
					d[k] = b[k] - L;
					continue;
				}

				if (b[k] < -L/2)
				{
					d[k] = b[k] + L;
					continue;
				}

				d[k] = b[k];
			}

			return d.sqr() <= rsqr;
		}

		virtual bool operator()(kdTree::BoundingBox<size_t, R> const &bb) const
		{
			double L = box->L();
			double p, q, z;

			/*for (unsigned k = 0; k < R; ++k)
			{
				if (a[k] > (bb.max_coord(k) + r) or a[k] < (bb.min_coord(k) - r))
					return false;
			}
			return true;*/

			for (unsigned k = 0; k < R; ++k)
			{
				if (bb.max_coord(k) - bb.min_coord(k) + 2*r >= L)
					continue;

				z = fmodulus(a[k], L);
				p = fmodulus(bb.min_coord(k) - r, L);
				q = fmodulus(bb.max_coord(k) + r, L);
				//z = fmodulus(a[k], L); p = bb.min_coord(k) - r; q = bb.max_coord(k) + r;

				if (p < q and (z <= p  or q <= z)) return false;
				if (q < p and (q <= z and z <= p)) return false;
			}

			return true;
		}
};

#include <cstdlib>
inline int random_colour(size_t seed)
{
	srand(seed);
	return rand() % 1022;
}

template <unsigned R>
void fof_run(BoxMaker const &box_, std::istream &fi, std::ostream &fo, double ll, FILE_FORMAT file_format = System::FMT_CONAN)
{
	auto box = box_.box<R>(); double L = box->L();
	auto pos = load_from_file<dVector<R>>(fi, "pos");

	// store the kdTree in the form of a permutated list of
	// integers pointing to the correct location in [pos].
	std::cerr << "[creating index array] ";
	Array<size_t> indices(pos.size());
	copy(Range<size_t>(pos.size()), indices);

	std::cerr << "[building kdTree] ";
	kdTree::Tree<size_t, R> tree(indices.begin(), indices.end(),
		[pos] (size_t idx, unsigned k) { return pos[idx][k]; });

	Array<size_t> tags(pos.size(), 0);
	size_t current = 1;

	std::cerr << "[scanning particle neighbours] " << std::endl;
	Misc::ProgressBar pb(pos.size());
	for (size_t idx = 0; idx < pos.size(); ++idx)
	{
		pb.tic();

		/*
		if (tags[idx] == 0)
		{
			tags[idx] = current;
			++current;
		}

		Proximity<R> P(box, pos, pos[idx], ll);
		Collector C;
		tree.traverse(C, P);

		for (size_t j : C.yield())
		{
			if (tags[j] != tags[idx])
				std::cerr << "!";
			tags[j] = tags[idx];
		}*/

		if (tags[idx] != 0) continue;

		flood_fill<size_t>(
			// find neighbours
			[&] (size_t j) {
			/*	Proximity<R> P2(box, pos, pos[j], ll*2);
				PrintVector<R> C2(pos, j);
				tree.traverse(C2, P2); */

				Proximity<R> P(box, pos, pos[j], ll);
				Collector C;		 
				tree.traverse(C, P);
				return C.yield();
			},

			// predicate, true if element needs attention
			[&] (size_t j) -> bool {
				if (tags[j] != current and tags[j] != 0)
					std::cerr << "!";
				return tags[j] != current;
			},

			// paint action
			[&] (size_t j) {
				tags[j] = current;
			},

			// start
			idx);
		++current;
	}
	pb.finish();

	std::cerr << "[found " << current-1 << " groups] [post-processing] ";
	Array<Array<dVector<R>>> fof_pos(0, Array<dVector<R>>(0));
	for (size_t i = 0; i < (current-1); ++i) 
		fof_pos->push_back(Array<dVector<R>>(0));

	Array<size_t>		 fof_mass(current-1, 0);

	for (size_t idx = 0; idx < pos.size(); ++idx)
	{
		size_t t = tags[idx];
		if (t == 0 or t >= current) throw "error";
		fof_pos[t-1]->push_back(pos[idx]);
		++fof_mass[t-1];
	}

	for (size_t t = 0; t < fof_pos.size(); ++t)
	{
		if (fof_mass[t] < 10) continue;
		fo << fof_pos[t][0] % L << " " << fof_mass[t] << " " << t << std::endl;
	}
	fo << "\n\n";

	for (size_t t = 0; t < fof_pos.size(); ++t)
	{
		//if (fof_mass[t] < 10) continue;
		for (auto p : fof_pos[t])
			fo << p << " " << random_colour(t) << std::endl;
	}

	std::cerr << "[done]\n";
}

void cmd_fof(int argc, char **argv)
{
	Argv C = read_arguments(argc, argv,
		Option({0, "h", "help", "false",
			"print help on the use of this program."}),

		Option({Option::VALUED | Option::CHECK, "i", "id", date_string(),
			"identifier for filenames."}),

		Option({Option::VALUED | Option::CHECK, "t", "time", "1.0",
			"time of snapshot."}),

		Option({Option::VALUED | Option::CHECK, "l", "linking-length", "0.1",
			"linking length in units of Mpc h^-1."}),

		Option({Option::VALUED | Option::CHECK, "f", "filename", "",
			"if given, supersedes internal naming system."}),

		Option({0, "", "ascii", "true",
			"give output in ASCII."})
	);

	if (C.get<bool>("help"))
	{
		std::cout << "Cosmic workset Conan, by Johan Hidding.\n\n";
		C.print(std::cout);
		exit(0);
	}

	/************************************************************
	| Input data
	+----------------------------------------------------------*/
	double t = C.get<double>("time");

	//std::ifstream fi0(timed_filename(C["id"], "density", -1.0, ".conan"));
	std::string fi_fn = (C["filename"] == "" ? 
		timed_filename(C["id"], "pos", t, ".conan") :
		C["filename"]);

	std::ifstream fi(fi_fn);
	Header H(fi); H << C;
	History I(fi); I << C;

	/************************************************************
	| Box
	+----------------------------------------------------------*/
	unsigned N  = H.get<unsigned>("N"),
		 R  = (H.count("dim") == 1 ? H.get<unsigned>("dim") : 3);
	double   L  = H.get<double>("L"),
		 ll = H.get<double>("linking-length");

	BoxMaker box(N, L);

	/************************************************************
	| Output data
	+----------------------------------------------------------*/
	std::ofstream fo(timed_filename(C["id"], "fof", t, ".txt"));
	H.to_txt_file(fo);
	I.to_txt_file(fo);

	switch (R)
	{
		case 2: fof_run<2>(box, fi, fo, ll);
			break;
		case 3: fof_run<3>(box, fi, fo, ll);
			break;
	}
}

Global<Command> _FOF("fof", cmd_fof);

