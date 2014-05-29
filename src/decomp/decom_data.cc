#include "decom_data.hh"
#include <algorithm>

using namespace Conan;

Decom_data_2::Decom_data_2(Cell_structure data[2])
{
	V.resize(4);

	// for all cells
	for (uint8_t i = 0; i < 2; ++i)
	{
		Vec Y(3);
		std::copy(data[i].vertices, data[i].vertices+3, Y.begin());

		// for all edges
		for_each_decrement(Y, [&] (Vec const &W, uint8_t j)
		{
			if (V[W[0]].count(W[1]) == 0)
			{
				Edge e;
				e.vertices = W;				// add vertices to edge
				uint8_t n = E.size();

				V[W[0]][W[1]] = n;			// add edge to vertices
				V[W[1]][W[0]] = n;

				E.push_back(e);
			}
		});

		Cell c;
		c.vertices = Y;
		uint8_t n = C.size();
		c.unit_volume = data[i].unit_volume;

		// for all edges, and opposed vertices
		for_each_decrement(Y, [&] (Vec const &W, uint8_t j)
		{
			V[j].cells.push_back(n);
			c.edges.push_back(V[W[0]][W[1]]);

			E[V[W[0]][W[1]]][j] = n;
		});

		C.push_back(c);
	}
}

Decom_data_3::Decom_data_3(Cell_structure data[5])
{
	// build vertices
	V.resize(8);

	// add cell
	for (uint8_t i = 0; i < 5; ++i)
	{
		Vec Y(4);			
		std::copy(data[i].vertices, data[i].vertices + 4, Y.begin());

		// ensure existence of faces
		for_each_decrement(Y, [&] (Vec const &X, uint8_t j)
		{
			// ensure existence of edges
			for_each_decrement(X, [&] (Vec const &W, uint8_t k)
			{
				if (V[W[0]].count(W[1]) == 0)
				{
					Edge e;
					e.vertices = W;				// add vertices to edge
					uint8_t n = E.size();
					E.push_back(e);

					V[W[0]][W[1]] = n;			// add edge to vertices
					V[W[1]][W[0]] = n;
				}
			});

			if (E[V[X[0]][X[1]]].count(X[2]) == 0)
			{
				Face f;
				uint8_t n = F.size();
				f.vertices = X;					// add vertices to face

				for_each_decrement(X, [&] (Vec const &W, uint8_t k)
				{
					f.edges.push_back(V[W[0]][W[1]]);	// add edges to face
					V[k].faces.push_back(n);		// add face to vertices
					E[V[W[0]][W[1]]][k] = n;		// add face to edges
				});

				F.push_back(f);
			}
		});
		
		Cell c;
		uint8_t n = C.size();
		c.vertices = Y;							// add vertices to cell
		c.unit_volume = data[i].unit_volume;

		for_each_decrement(Y, [&] (Vec const &X, uint8_t j)
		{
			c.faces.push_back(E[V[X[0]][X[1]]][X[2]]);		// add faces to cell

			for_each_decrement(X, [&] (Vec const &W, uint8_t k)
			{
				c.edges.push_back(V[W[0]][W[1]]);		// add edges to cell
				E[V[W[0]][W[1]]].cells.push_back(n);		// add cell to edges
			});

			V[j].cells.push_back(n);				// add cell to vertices
			F[E[V[X[0]][X[1]]][X[2]]][j] = n;			// add cell to faces
		});

		C.push_back(c);
	}
}

