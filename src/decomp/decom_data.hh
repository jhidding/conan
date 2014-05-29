#pragma once
#include <vector>
#include <map>
#include "numeric.hh"

/* defines data formats for the decomposition
 * of squares and cubes.
 */

namespace Conan
{
	struct Decom_data_2
	{
		typedef std::vector<uint8_t> Vec;
		typedef std::map<uint8_t, uint8_t> Map;

		struct Cell_structure
		{
			double unit_volume;
			uint8_t vertices[3];
		};

		class Cell
		{
			public:
			double unit_volume;
			Vec vertices;
			Vec edges;
		};

		class Edge: public Map
		{
			public:
			Vec vertices;
		};

		class Vertex: public Map
		{
			public:
			Vec cells;
		};

		std::vector<Cell> C;
		std::vector<Edge> E;
		std::vector<Vertex> V;

		Decom_data_2(Cell_structure data[2]);
	};

	struct Decom_data_3
	{
		typedef std::vector<uint8_t> Vec;
		typedef std::map<uint8_t, uint8_t> Map;

		struct Cell_structure
		{
			double unit_volume;
			uint8_t vertices[4];
		};

		class Cell
		{
			public:
			double unit_volume;
			Vec vertices;
			Vec edges;
			Vec faces;
		};

		class Face: public Map
		{
			public:
			Vec vertices;
			Vec edges;
		};

		class Edge: public Map
		{
			public:
			Vec vertices;
			Vec cells;
		};

		class Vertex: public Map
		{
			public:
			Vec faces;
			Vec cells;
		};

		std::vector<Cell> C;
		std::vector<Face> F;
		std::vector<Edge> E;
		std::vector<Vertex> V;

		Decom_data_3(Cell_structure data[5]);
	};
}

