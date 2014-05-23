#pragma once
#include "system.hh"

namespace Conan
{
	class Solver
	{
		public:
			virtual void init(double a_pos, double a_mom) = 0;
			virtual void kick(double a, double da) = 0;
			virtual void drift(double a, double da) = 0;

			virtual void save_snapshot(unsigned i) const
			{}
	};

	class Integrator
	{
		public:
			virtual void run(ptr<Solver>) const = 0;
	};

	class LeapFrog: public Integrator
	{
		double   a0, da;
		unsigned n;

		public:
			LeapFrog(double a0_, double da_, unsigned n_):
				a0(a0_), da(da_), n(n_)
			{}

			void run(ptr<Solver> S) const
			{
				S->init(a0, a0 + da/2);

				double a = a0;
				for (unsigned i = 0; i < n; ++i)
				{
					S->drift(a, da);
					S->kick(a + da/2, da);
					a += da;

					S->save_snapshot(i);
				}
			}
	};
}

