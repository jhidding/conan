#include "../base/system.hh"
#include "../base/format.hh"

#include "cic.hh"

#include <fstream>

using namespace System;
using namespace Conan;

void cmd_cic(int argc, char **argv)
{
	std::ostringstream ss;
	ss << time(NULL);
	std::string timed_seed = ss.str();

	Argv C = read_arguments(argc, argv,
		Option({0, "h", "help", "false",
			"print help on the use of this program."}),

		Option({Option::VALUED | Option::CHECK, "i", "id", date_string(),
			"identifier for filenames."}),

		Option({Option::VALUED | Option::CHECK, "f", "f", "2",
			"fractional Eulerian resolution. The default makes the"
			" Eulerian grid twice the size of the mass grid."}),

		Option({Option::VALUED | Option::CHECK, "m", "m", "5",
			"snap-shot"}),

		Option({Option::VALUED | Option::CHECK, "s", "sub-sample", "1",
			"subsample the grid and interpolate the position."}));

	if (C.get<bool>("help"))
	{
		std::cout << "Cosmic workset Conan, by Johan Hidding.\n\n";
		C.print(std::cout);
		exit(0);
	}

	/************************************************************
	| Input data
	+----------------------------------------------------------*/
	unsigned m = C.get<unsigned>("m");
	std::ifstream fi(Misc::format(C["id"], ".pos.", m, ".conan")); 
	Header H(fi); H << C;
	History I(fi); I << C;

	/************************************************************
	| Box
	+----------------------------------------------------------*/
	unsigned f = H.get<unsigned>("f"),
		 N = H.get<unsigned>("N");
	double   L = H.get<double>("L");
	unsigned s = H.get<unsigned>("sub-sample");

	if (H.count("dim") == 0)
		H["dim"] = "3";

	BoxMaker mass_box(N, L),
		 force_box(f*N, L);

	std::ofstream fo;
	fo.open(Misc::format(C["id"], ".rho.", m, ".conan"));
	H["N"] = System::to_string(N*f);
	H.to_file(fo);
	I.to_file(fo);

	switch (H.get<int>("dim"))
	{
		case 1: cic_run<1>(fi, fo, mass_box, force_box, s); 
			break;

		case 2: cic_run<2>(fi, fo, mass_box, force_box, s);
			break;

		case 3: cic_run<3>(fi, fo, mass_box, force_box, s);
			break;
	}

	fo.close();
}

Global<Command> _CMD_CIC("cic", cmd_cic);

