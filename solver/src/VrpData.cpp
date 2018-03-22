/************************************************************
  Autore: Andrea Maggiordomo - mggndr89@gmail.com
  Università di Pisa - Informatica Applicata - Logistica
  2012
*************************************************************/

#include "VrpData.hpp"

#include <iostream>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;
using std::stack;
using std::max;
using std::vector;
using std::abs;

/** Constructor */
VrpData::VrpData( int M, int N, int Cw, int Cv, int *W, int *V, int *vBase ) :
	M( M ),
	N( N ),
	Cw( Cw ),
	Cv( Cv ),
	W( W, W+N ),
	V( V, V+N ),
	vBase( vBase, vBase+N-1 )
{
	VrpData::G = NULL;
}

/** Copy constructor */
VrpData::VrpData( const VrpData &rhs ) :
	M( rhs.M ),
	N( rhs.N ),
	Cw( rhs.Cw ),
	Cv( rhs.Cv ),
	W( rhs.W ),
	V( rhs.V ),
	vBase( rhs.vBase )
{
	if ( rhs.G != NULL ) {
		VrpData::G = new double [N*N];

		for ( int i = 0; i < N*N; i++ )
			VrpData::G[i] = rhs.G[i];
	}

	else
		VrpData::G = NULL;
}

/** Assignment operator */
VrpData & VrpData::operator=( const VrpData &rhs )
{
	if ( this != &rhs ) {
		M = rhs.M;
		N = rhs.N ;
		Cw = rhs.Cw;
		Cv = rhs.Cv,
		W = rhs.W;
		V = rhs.V;
		vBase = rhs.vBase;

		if ( G != NULL ) {
			delete [] G;
			G = NULL;
		}

		if ( rhs.G != NULL ) {
			G = new double [N*N];
			for ( int i = 0; i < N*N; i++ )
				G[i] =rhs.G[i];
		}
	}

	return *this;
}

/** Destructor */
VrpData::~VrpData()
{
	if ( G != NULL )
		delete [] G;
}

/** Given a fractional solution from a solver builds the residual graph
	stored in net */
void VrpData::buildResidualGraph( const double *colSolution )
{
	if ( G != NULL ) {
		cout << "Warning: deleting pre-existing net (VrpGraph)" << endl;
		delete [] G;
	}

	G = new double[N*N];

	// clear all values
	for ( int i = 0; i < N*N; i++ )
		G[i] = 0.0;

	// assign values
	for ( int i = 0; i < N-1; i++ ) {

		for( int j = i+1; j < N; j++ ) {
			// edges i-j
			int offset = j - i - 1;
			int colId = vBase[i] + offset;
			G[N*i+j] = colSolution[colId];
		}
	}
}

double VrpData::getCutValue( VertexSet &S )
{
	VertexSet *cmpl = getComplement( &S );

	VertexSet& out = *cmpl;

	double cutValue = 0.0;
	for ( VertexSet::size_type i = 0; i < S.size(); i++ ) {
		for ( VertexSet::size_type j = 0; j < out.size(); j++ ) {
			cutValue += G[ N*S[i] + out[j] ]; // cut out-going
			cutValue += G[ N*out[j] + S[i] ]; // cut in-going
		}
	}

	delete cmpl;

	return cutValue;
}

/** Depth first search on the stored graph */
VertexSet *VrpData::dfs( int source )
{
	if ( VrpData::G == NULL ) {
		cerr << "Error: No graph loaded (dfs)" << endl;
		exit( EXIT_FAILURE );
	}

	return dfs( VrpData::G, source, VrpData::N );
}

/** Depth first search on the residual graph derived from a partial
solution */
VertexSet *VrpData::dfs( const double *Gcap, int source, int n )
{
	if ( source < 0 || source >= n ) {
		cerr << "dfs: invalid source node" << endl;
		exit( EXIT_FAILURE );
	}

	stack< int, VertexSet > s;
	bool *seen = new bool[n];
	VertexSet *visited = new VertexSet;

	for ( int i = 0; i < n; i++ )
		seen[i] = false;

	s.push( source );
	
	while ( !s.empty() ) {
		int u = s.top();
		s.pop();

		if ( !seen[u] ) {

			seen[u] = true;
			visited->push_back( u );

			for ( int i = 0; i < n; i++ ) {
				if ( ( Gcap[n*u+i] || Gcap[n*i+u] ) && !seen[i] )
					/* (u,i) or (i,u) exists and i not yet visited */
					s.push( i );
			}
		}
	}

	delete [] seen;

	/* Now return visited set */
	return visited;
}

/** Get complement of vertex set S */
VertexSet *VrpData::getComplement( VertexSet *S )
{
	bool *inS = new bool[N];

	for ( int i = 0; i < N; i++ )
		inS[i] = false;

	for ( VertexSet::size_type i = 0; i < S->size(); i++ ) {
		if ( (*S)[i] >= N ) {
			cerr << "getComplement: found a node out of range" << endl;
			exit( EXIT_FAILURE );
		}
		else
			inS[ (*S)[i] ] = true;
	}

	/* Update size */
	int complementSize = N - S->size();

	if ( complementSize == 0 ) {
		delete [] inS;
		return NULL;
	}

	VertexSet *cpl = new VertexSet( complementSize );
	int i, j;
	for ( i = j = 0; i < N; i++ )
		if ( !inS[i] )
			cpl->at(j++) = i;

	delete [] inS;

	return cpl; 
}

/** Row lb = ceiling( ( sum of demands in S ) / ( capacity ) ) */
int VrpData::getSetLb( VertexSet *S )
{
	int weights = 0, volumes = 0;
	for ( VertexSet::size_type i = 0; i < S->size(); i++ ) {
		weights += W[ (*S)[i] ];
		volumes += V[ (*S)[i] ];
	}
	
	double fracLbWeight = double( weights ) / double( Cw );
	double fracLbVolume = double( volumes ) / double( Cv );
	return int( ceil( max( fracLbWeight, fracLbVolume ) ) );
}

/** This must be called only if the solution is connected and integer */
RouteSet VrpData::getRoutes( double eps )
{
	if ( G == NULL ) {
		cerr << "getRoutes: no solution to examine, returning empty set" << endl;
		return RouteSet();
	}

	bool *taken = new bool[N];
	for ( int i = 0; i < N; i++ )
		taken[i] = false;

	// each node should have only two inciding arcs
	RouteSet R;

	for ( int i = 0; i < N; i++ ) {
		// scan forward star of the depot
		if ( ( abs(G[N*0+i]-2.0)<eps || abs(G[N*i+0]-2.0)<eps ) && !taken[i] ) {
			// only one node in this route
			Route ri;
			ri.push_back( i );
			taken[i] = true;

			R.push_back( ri );
		}
		else if ( ( abs(G[N*0+i]-1.0)<eps || abs(G[N*i+0]-1.0)<eps ) && !taken[i] ) {

			Route ri;
			int d = i;
			int from = 0;

			while ( d != 0 ) {
				// add d
				ri.push_back( d );
				taken[d] = true;

				// find next node
				for ( int j = 0; j < N; j++ ) {
					if ( j != from && ( abs(G[N*d+j]-1.0)<eps || abs(G[N*j+d]-1.0)<eps ) ) {
						from = d;
						d = j;
						break;
					}
				}
			}

			R.push_back( ri );
		}
	} // for each (0,i) in sol

	delete [] taken;

	return R;
}

/** Tell if a given set of routes is feasible */
bool VrpData::isFeasible( const RouteSet &S )
{
	if ( S.empty() || S.size() > RouteSet::size_type( M ) )
		return false;

	bool feasible = true;
	bool all = true;
	vector<bool> taken(N);
	taken[0] = true;
	for ( int i = 1; i < N; i++ )
		taken[i] = false;

	for ( RouteSet::size_type i = 0; i < S.size() && feasible; i++ ) {
		// for each route
		int wLoad = 0, vLoad = 0;

		for ( Route::size_type j = 0; j < S[i].size(); j++ ) {
			wLoad += W[ S[i][j] ];
			vLoad += V[ S[i][j] ];
			taken[ S[i][j] ] = true;
		}

		feasible = ( wLoad <= Cw ) && ( vLoad <= Cv );
	}

	if ( feasible )
		for ( int i = 1; i < N && all; i++ )
			all = taken[i];

	return feasible && all;
}
