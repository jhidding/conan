#pragma once
#include "system.hh"
#include "numeric.hh"

#include <iostream>
#include <fstream>

namespace Conan
{
	using System::BoxMaker;

	template <unsigned R>
	void print_array(BoxPtr<R> box, std::ostream &fo, Array<double> f)
	{
		size_t o = 0;
		for (size_t i = 0; i < box->N(); ++i)
		{
			for (size_t j = 0; j < box->N(); ++j, ++o)
			{
				fo << f[o] << " ";
			}
			fo << std::endl;
		}
	}

	template <unsigned R>
	class Zeldovich
	{
		BoxPtr<R> box;
		Array<dVector<R>> u;

		public:
			Zeldovich(BoxPtr<R> box_, Array<double> phi):
				box(box_),
				u(box_->size())
			{
				gradient_2nd_order(box, u, phi);
			}

			void p(Array<dVector<R>> p_, double a)
			{
				for (size_t i = 0; i < box->size(); ++i)
					p_[i] = u[i] * -a;
			}

			void x(Array<dVector<R>> x_, double a)
			{
				for (size_t i = 0; i < box->size(); ++i)
					x_[i] = box->G[i] - u[i] * a;
			}
	};

	template <unsigned R>
	class Gravity: public Solver
	{
		BoxPtr<R>	mbox, fbox;
		double		mass;
		Array<double> 	phi, delta;
		Cosmology	cosmos;

		Fourier::Transform fft;

		Array<dVector<R>> X, P;
		Gradient<R>       A_l;

		LinearInterpolation<Gradient<R>, R> A;

		std::string 	id;

		public:
			Gravity(std::string const &id_, BoxPtr<R> mbox_, BoxPtr<R> fbox_, 
					Array<double> phi_, Cosmology cosmos_):
				mbox(mbox_), fbox(fbox_), phi(phi_), delta(fbox->size()), 
				cosmos(cosmos_), fft(std::vector<int>(R, fbox->N())), 
				X(mbox->size()), P(mbox->size()), A_l(fbox, delta),
				A(fbox, A_l), id(id_)
			{
				mass = pow(double(fbox->N()) / mbox->N(), R);
			}

			virtual void init(double a_pos, double a_mom)
			{
				Zeldovich<R> ZA(mbox, phi);
				ZA.x(X, a_pos);
				ZA.p(P, a_mom);
				phi->resize(fbox->size());
			}

			virtual void kick(double a, double da)
			{
				md_cic(fbox, delta, X, mass);

				std::cout << "."; std::cout.flush();

				transform(delta, fft.in, [] (double f)
					{ return f - 1.0; });
				fft.forward();
				p_transform(fft.out, fbox->K, fft.in, Fourier::Fourier<R>::filter(
					Fourier::Fourier<R>::potential()));
				fft.in[0] = 0;
				fft.backward();
				p_transform(fft.out, delta, Fourier::real_part(fbox->size() / cosmos.G() * a));

				double adot = cosmos(a);
				#pragma omp parallel for
				for (size_t i = 0; i < mbox->size(); ++i)
				{
					P[i] -= A(X[i]) * (da / adot);
				}
			}

			virtual void drift(double a, double da)
			{
				double adot = cosmos(a);
				
				#pragma omp parallel for
				for (size_t i = 0; i < mbox->size(); ++i)
				{
					X[i] += P[i] * (da / (a * a * adot));
				}
			}

			virtual void save_snapshot(unsigned i) const
			{
				if ((i+1) % 10 == 0)
				{
					std::cout << "|"; std::cout.flush();
					std::ofstream fo(Misc::format(id, ".pos.", (i+1)/10, ".conan"));
					for (size_t i = 0; i < mbox->size(); ++i)
					{
						fo << X[i] << " " << P[i] << std::endl;
					}
					fo.close();
				}
			}
	};

	template <unsigned R>
	void nbody_run(std::string const &id, BoxMaker mass_box_, BoxMaker force_box_,
		Cosmology cosmos, Integrator const &integr, Array<double> phi)
	{
		auto mass_box  =  mass_box_.box<R>(),
		     force_box = force_box_.box<R>();
		     
		ptr<Solver> solver(new Gravity<R>(id, mass_box, force_box, phi, cosmos));
		integr.run(solver);
	}
}

