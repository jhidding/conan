#include "../base/system.hh"
#include "../base/format.hh"

#include "cosmology.hh"
#include "leapfrog.hh"
#include "nbody.hh"

#include <fstream>

using namespace System;
using namespace Conan;

void cmd_run(int argc, char **argv)
{
	std::ostringstream ss;
	ss << time(NULL);
	std::string timed_seed = ss.str();

	Argv C = read_arguments(argc, argv,
		Option({0, "h", "help", "false",
			"print help on the use of this program."}),

		Option({Option::VALUED | Option::CHECK, "i", "id", date_string(),
			"identifier for filenames."}),

		Option({Option::VALUED | Option::CHECK, "H0", "H0", "68.5",
			"Hubble parameter, default is point of overlap between "
			"WMAP9 and Planck."}),

		Option({Option::VALUED | Option::CHECK, "m", "Omega_m", "1.0",
			"Omega-matter. The default universe is EdS."}),

		Option({Option::VALUED | Option::CHECK, "L", "Omega_L", "0.0",
			"Omega-Lambda. The default universe is EdS."}),

		Option({0, "c", "concordance", "false",
			"Override cosmological parameters to have concordance (Planck) model: "
			"H0 = 68.5, Omega_L = 0.69, Omega_m = 0.31"}),

		Option({Option::VALUED | Option::CHECK, "f", "f", "2",
			"fractional force resolution. The default makes the"
			" force grid twice the size of the mass grid."}),

		Option({Option::VALUED | Option::CHECK, "a0", "a0", "0.02",
			"start time."}),

		Option({Option::VALUED | Option::CHECK, "da", "da", "0.02",
			"size of time step."}),

		Option({Option::VALUED | Option::CHECK, "n", "n", "50",
			"number of iterations."}));

	if (C.get<bool>("help"))
	{
		std::cout << "Cosmic workset Conan, by Johan Hidding.\n\n";
		C.print(std::cout);
		exit(0);
	}

	/************************************************************
	| Input data
	+----------------------------------------------------------*/	
	std::ifstream fi(C["id"] + ".density.init.conan"); 
	Header H(fi); H << C;
	History I(fi); I << C;
	Array<double> phi = load_from_file<double>(fi, "potential");
	fi.close();

	/************************************************************
	| Box
	+----------------------------------------------------------*/
	unsigned f = H.get<unsigned>("f"),
		 N = H.get<unsigned>("N");
	double   L = H.get<double>("size");

	BoxMaker mass_box(N, L),
		 force_box(f*N, L);

	/************************************************************
	| Cosmology
	+----------------------------------------------------------*/
	double  H0      = H.get<double>("H0"),
		Omega_m = H.get<double>("Omega_m"),
		Omega_L = H.get<double>("Omega_L");

	if (C.get<bool>("concordance"))
	{
		H0      = 68.5;
		Omega_m = 0.31;
		Omega_L = 0.69;
	}

	Cosmology  cosmos(H0, Omega_m, Omega_L);

	/************************************************************
	| Integrator
	+----------------------------------------------------------*/
	double  a0 = H.get<double>("a0"),
		da = H.get<double>("da"),
		n  = H.get<double>("n");

	LeapFrog   integr(a0, da, n);

	switch (H.get<int>("dim"))
	{
		case 1: nbody_run<1>(H["id"], 
				mass_box, force_box, 
				cosmos, integr, phi); break;

		case 2: nbody_run<2>(H["id"], 
				mass_box, force_box, 
				cosmos, integr, phi); break;

		case 3: nbody_run<3>(H["id"], 
				mass_box, force_box, 
				cosmos, integr, phi); break;
	}
}

Global<Command> _RUN("run", cmd_run);

