#pragma once
#include "system.hh"
#include "numeric.hh"

#include <iostream>
#include <fstream>

namespace Conan
{
	enum FILE_FORMAT { FMT_ASCII, FMT_IFRIT, FMT_CONAN };

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
		bool            ff;
		Array<unsigned> flips;
		Cosmology	cosmos;

		Fourier::Transform fft;

		Array<dVector<R>> X, P;
		Gradient<R>       A_l;

		LinearInterpolation<Gradient<R>, R> A;

		std::string 	id;
		FILE_FORMAT	format;

		public:
			Gravity(std::string const &id_, BoxPtr<R> mbox_, BoxPtr<R> fbox_, 
					Array<double> phi_, Cosmology cosmos_, bool ff_,
					FILE_FORMAT format_):
				mbox(mbox_), fbox(fbox_), phi(phi_), delta(fbox->size()),
			       	ff(ff_), flips(ff ? mbox->size() : 0),
				cosmos(cosmos_), fft(std::vector<int>(R, fbox->N())), 
				X(mbox->size()), P(mbox->size()), A_l(fbox, delta),
				A(fbox, A_l), id(id_), format(format_)
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

			void save_ifrit(std::ostream &fo) const
			{
				if (R != 3) throw "IFRIT files only support 3D";

				std::vector<uint32_t> N(1, mbox->size());
				System::write_block_32(fo, N);

				std::vector<float> dim(6, 0.0);
				for (unsigned k = 0; k < 3; ++k)
					dim[k+3] = mbox->L();
				System::write_block_32(fo, dim);

				std::vector<float> data(mbox->size());
				for (unsigned k = 0; k < 3; ++k)
				{
					#pragma omp parallel for
					for (size_t i = 0; i < mbox->size(); ++i)
					{
						data[i] = X[i][k];
					}
					System::write_block_32(fo, data);
				}
				for (unsigned k = 0; k < 3; ++k)
				{
					#pragma omp parallel for
					for (size_t i = 0; i < mbox->size(); ++i)
					{
						data[i] = P[i][k];
					}
					System::write_block_32(fo, data);
				}
			}

			virtual void save_snapshot(unsigned i) const
			{
				if ((i+1) % 10 == 0)
				{
					std::cout << "|"; std::cout.flush();
					std::string fn;
					switch (format)
					{
						case FMT_ASCII: fn = Misc::format(id, ".pos.", (i+1)/10, ".txt");
								break;
						case FMT_IFRIT: fn = Misc::format(id, ".pos.", (i+1)/10, ".bin");
								break;
						case FMT_CONAN: fn = Misc::format(id, ".pos.", (i+1)/10, ".conan");
								break;
					}

					std::ofstream fo(fn);
					switch (format)
					{
						case FMT_ASCII:
							for (size_t i = 0; i < mbox->size(); ++i)
							{
								fo << X[i] << " " << P[i] << std::endl;
							} break;
						case FMT_IFRIT: 
							save_ifrit(fo);
							break;
						case FMT_CONAN:
							break;
					}
					fo.close();
				}
			}
	};

	template <unsigned R>
	void nbody_run(std::string const &id, BoxMaker mass_box_, BoxMaker force_box_,
		Cosmology cosmos, Integrator const &integr, Array<double> phi, 
		FILE_FORMAT format = FMT_CONAN, bool ff = false)
	{
		auto mass_box  =  mass_box_.box<R>(),
		     force_box = force_box_.box<R>();
		     
		ptr<Solver> solver(new Gravity<R>(id, mass_box, force_box, phi, cosmos, ff, format));
		integr.run(solver);
	}
}

