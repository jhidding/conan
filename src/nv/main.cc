#include "../base/system.hh"
#include "../base/format.hh"

#include "nv.hh"
#include "sugar.hh"

#include <fstream>

using namespace System;
using namespace Conan;

void cmd_nv(int argc, char **argv)
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

		Option({Option::VALUED | Option::CHECK, "D", "growing-mode", "1.0",
			"growing-mode solution."}),

		Option({0, "", "ascii", "false",
			"give output in ASCII (not recommended for 3D)."}),
		
		Option({0, "", "no-enh", "false",
			"do not enhance on sub-grid scale."}));

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
	double   D = H.get<double>("growing-mode");
	bool     enhance = not H.get<bool>("no-enh");

	BoxMaker mass_box(N, L),
		 force_box(f*N, L);

	FILE_FORMAT format = FMT_CONAN;
	if (H.get<bool>("ascii")) format = FMT_ASCII;
	std::ofstream fo;

	switch (format)
	{
		case FMT_CONAN:
			fo.open(System::timed_filename(H["id"], "flt", D));
			H["N"] = System::to_string(N*f);
			H.to_file(fo);
			I.to_file(fo);
			break;

		case FMT_ASCII:
			fo.open(System::timed_filename(H["id"], "flt", D, ".txt"));
			break;
	}

	switch (H.get<int>("dim"))
	{
		case 1: nv_run<1>(fo, 
				mass_box, force_box, 
				phi, D, format, enhance); break;

		case 2: nv_run<2>(fo, 
				mass_box, force_box, 
				phi, D, format, enhance); break;

		case 3: nv_run<3>(fo,
				mass_box, force_box, 
				phi, D, format, enhance); break;
	}

	fo.close();
}

Global<Command> _CMD_NV("nv", cmd_nv);

void cmd_sugar(int argc, char **argv)
{
	std::ostringstream ss;
	ss << time(NULL);
	std::string timed_seed = ss.str();

	Argv C = read_arguments(argc, argv,
		Option({0, "h", "help", "false",
			"print help on the use of this program."}),

		Option({Option::VALUED | Option::CHECK, "i", "id", date_string(),
			"identifier for filenames."}),

		Option({Option::VALUED | Option::CHECK, "D", "growing-mode", "1.0",
			"growing-mode solution."}),

		Option({0, "", "ascii", "false",
			"give output in ASCII (not recommended for 3D)."}));
		
	if (C.get<bool>("help"))
	{
		std::cout << "Cosmic workset Conan, by Johan Hidding.\n\n";
		C.print(std::cout);
		exit(0);
	}

	/************************************************************
	| Input data
	+----------------------------------------------------------*/
	std::ifstream fi(timed_filename(C["id"], "flt", C.get<double>("growing-mode")));
	std::ofstream fo(timed_filename(C["id"], "dens",  C.get<double>("growing-mode")));
	Header H(fi); H << C;
	History I(fi); I << C;

	/************************************************************
	| Box
	+----------------------------------------------------------*/
	unsigned N = H.get<unsigned>("N");
	double   L = H.get<double>("size");
	double   D = H.get<double>("growing-mode");

	BoxMaker box(N, L);

	FILE_FORMAT format = FMT_CONAN;
	if (H.get<bool>("ascii")) format = FMT_ASCII;

	if (format == FMT_CONAN)
	{
		H.to_file(fo);
		I.to_file(fo);
	}

	switch (H.get<int>("dim"))
	{
	//	case 1: sg_run<1>(fi, fo, box, format); break;

		case 2: sg_run<2>(fi, fo, box, format); break;

		case 3: sg_run<3>(fi, fo, box, format); break;
	}

	fi.close();
	fo.close();
}

Global<Command> _CMD_SUGAR("sugar", cmd_sugar);

