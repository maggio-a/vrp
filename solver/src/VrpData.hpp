/************************************************************
  Autore: Andrea Maggiordomo - mggndr89@gmail.com
  Università di Pisa - Informatica Applicata - Logistica
  2012
*************************************************************/

#ifndef VrpData_H
#define VrpData_H

#include "VrpTypes.hpp"

#include <vector>

class VrpData {

public:
	int M;
	int N;
	int Cw;
	int Cv;
	std::vector<int> W;
	std::vector<int> V;

	std::vector<int> vBase;

	double *G;

	/** Constructor */
	VrpData( int M, int N, int Cw, int Cv, int *W, int *V, int *vBase );

	/** Copy constructor */
	VrpData( const VrpData & );

	/** Assignment operator */
	VrpData & operator=( const VrpData & );

	/** Destructor */
	~VrpData();

	/** Given a fractional solution from a solver builds the residual graph
	stored in net */
	void buildResidualGraph( const double *colSolution );

	double getCutValue( VertexSet &S );

	/** Depth first search on the stored graph */
	VertexSet *dfs( int source );

	/** Depth first search on the residual graph derived from a partial
	solution */
	static VertexSet *dfs( const double *Gcap, int source, int n );

	/** Get complement of vertex set S */
	VertexSet *getComplement( std::vector<int> *S );

	/** Row lb = ( sum of demands in S ) / ( capacity ) */
	int getSetLb( VertexSet *S );

	/** Retrieve routes from solution */
	RouteSet getRoutes( double eps=0 );

	/** Tell if a given set of routes is feasible */
	bool isFeasible( const RouteSet &S );
};

#endif
