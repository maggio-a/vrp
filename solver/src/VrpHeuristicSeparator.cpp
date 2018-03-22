/************************************************************
  Autore: Andrea Maggiordomo - mggndr89@gmail.com
  Università di Pisa - Informatica Applicata - Logistica
  2012
*************************************************************/

#include "VrpHeuristicSeparator.hpp"

#include "OsiSolverInterface.hpp"
#include "CglTreeInfo.hpp"
#include "CglCutGenerator.hpp"

#include "CglGomory.hpp"
#include "CoinTime.hpp"

#include "VrpData.hpp"
#include "VrpTypes.hpp"

#include <stack>
#include <vector>

using std::vector;
using std::min;
using std::max;
using std::sort;

/** Constructor */
VrpSep::VrpSep( const VrpData &g, bool gm ) :
	CglCutGenerator(),
	graph( g ),
	gomoryTime( 0.0L ),
	gomoryTotal( 0 ),
	gomoryEffective( 0 ),
	useGomory( gm )
{
	numCols = (graph.N)*(graph.N-1)/2;
}

/* Copy constructor */
VrpSep::VrpSep( const VrpSep &src ) :
	CglCutGenerator(),
	graph( src.graph ),
	gomoryTime( src.gomoryTime ),
	gomoryTotal( src.gomoryTotal ),
	gomoryEffective( src.gomoryEffective ),
	useGomory( src.useGomory ),
	numCols( src.numCols )
{
}

/** Clone */
VrpSep *VrpSep::clone() const
{
	return new VrpSep( *this );
}

/** Assignment operator */ 
VrpSep & VrpSep::operator=( const VrpSep& rhs )
{
	if ( this != &rhs ) {
		CglCutGenerator::operator=( rhs );
		graph = rhs.graph;
	}

	return *this;
}

/** Returns true if needs optimal basis to do cuts */
bool VrpSep::needsOptimalBasis() const
{
  return false;
}

/** Remove duplicate cuts */
void VrpSep::cleanCutSet( vector<VertexSet*> *cuts ) const
{
	for ( vector<VertexSet*>::size_type i = 0; i < cuts->size(); i++ )
		sort( (*cuts)[i]->begin(), (*cuts)[i]->end() );
	
	VertexSet::size_type ilast = 0;
	do {
		bool found = false;
		vector<VertexSet*>::size_type i, ind;
		for ( i = ilast+1; i < cuts->size() && !found; i++ ) {
			if ( (*cuts)[ilast]->size() == (*cuts)[i]->size() 
				&& *((*cuts)[ilast]) == *((*cuts)[i]) ) {
					found = true;
					ind = i;
					break;
			}
		}
	
		if ( found ) {
			delete (*cuts)[ind];
			cuts->erase( cuts->begin() + ind );
		}
		else
			ilast++;
	
	} while ( ilast < cuts->size() );
}

/** First of all identify if there are any connected components detached
	from the depot, in which case cut them and resolve */
void VrpSep::detachedComponents( vector<VertexSet*> *cuts ) const
{
	/* Depth first search from the depot */
	VertexSet *visited = graph.dfs( 0 );

	VertexSet *detached = graph.getComplement( visited );

	delete visited;
	
	if ( detached != NULL ) {
		/* There are some connected components detached from the depot */
		for ( VertexSet::size_type i = 0; i < detached->size(); i++ ) {

			if ( (*detached)[i] == -1 )
				continue;

			VertexSet *c = graph.dfs( (*detached)[i] );
			sort( c->begin(), c->end() );

			cuts->push_back( c );

			/* Remove found component from detached */
			for ( VertexSet::size_type j = 0; j < c->size(); j++ )
				for ( VertexSet::size_type h = 0; h < detached->size(); h++ )
					if ( (*detached)[h] == (*c)[j] )
						(*detached)[h] = -1;
		}

		delete detached;
	}
}

void VrpSep::ccHeuristic( vector<VertexSet*> *cuts, double eps, bool findAll ) const
{
	vector<VertexSet*>::size_type nCuts = cuts->size();

	/* New graph without the depot */;
	double *Gcap = new double[(graph.N-1)*(graph.N-1)];
	for ( int i = 0; i < graph.N-1; i++ )
		for ( int j = 0; j < graph.N-1; j++ )
			Gcap[(graph.N-1)*i + j] = graph.G[graph.N*(i+1) + (j+1)];

	/* Unexplored nodes of Gcap */
	bool *seen= new bool[graph.N-1];
	for ( int i = 0; i < graph.N-1; i++ )
		seen[i] = false;

	for ( int i = 0; i < graph.N-1 && ( cuts->size() < 10 || findAll ); i++ ) {
		if ( seen[i] )
			continue;

		/* For each node of the altered residual graph not seen */

		/* Find the connected component */
		VertexSet *S = graph.dfs( Gcap, i, graph.N-1 );

		for ( VertexSet::size_type j = 0; j < S->size(); j++ ) {
			seen[ (*S)[j] ] = true; /* Mark the nodes as seen */
			(*S)[j]++; /* set is retrieved from a dfs on Gcap
			              node j in Gcap is node j+1 in G, so increase vals */
		}

		double cutVal = graph.getCutValue( *S );
		double SLB = 2*graph.getSetLb( S );

		if ( SLB - cutVal > eps ) {
			cuts->push_back( S );
		}
		else {
			/* Start iterations to find a subset of S that is likely to induce a
			   violated capacity constraint */
			bool inserted = false;
			while ( S->size() >= 2 && !inserted ) {
				int bestIndex;
				bool found = false; // if a candidate node is found
				for ( VertexSet::size_type j = 0; j < S->size(); j++ ) {
					VertexSet Snew( *S );
					Snew.erase( Snew.begin() + j );
					double nv = graph.getCutValue( Snew );					
					if ( nv < cutVal && 2*graph.getSetLb( &Snew ) == SLB ) {
						cutVal = nv;
						bestIndex = j;
						found = true;
					}
				}

				if ( !found ) {
					/* other sub components of S are unlikely to violate constraints */
					//delete S;
					break;
				}
				else {
					S->erase( S->begin() + bestIndex );
					// check if now S does violate smth
					if ( SLB - cutVal > eps ) {
						cuts->push_back( S );
						inserted = true;
						break; // break the while loop
					}
					/* else do nothing, keep iterating until there are
					   no more candidates or a violated cut is found */
				}
			} // while

			if ( !inserted ) // S not in cutSet, we can delete it
				delete S;
		}
	} // for

	delete [] seen;
	delete [] Gcap;

	// if new cuts found, delete duplicates
	if ( cuts->size() > nCuts )
		cleanCutSet( cuts );
}

void VrpSep::shrinkingHeuristic( vector<VertexSet*> *cuts, double eps, bool findAll ) const
{
	vector<VertexSet*>::size_type nCuts = cuts->size();

	/* Build G_i */
	int N = graph.N;
	const double *G = graph.G;

	int nNodes = N;
	double *G_i = new double[N*N];
	for ( int i = 0; i < N*N; i++ )
		G_i[i] = G[i];

	/* Initially each supernode contains only its primitive node */
	vector< VertexSet > superNode( N );
	for ( int i = 0; i < N; i++ )
		superNode[i].push_back( i );


	while ( cuts->size() < 10 || findAll ) {
		/* Step 1 
		   Find connected couples of supernodes that induce a violated capacity
		   constraint */
		for ( int i = 1; i < nNodes; i++ ) {
			for ( int j = i+1; j < nNodes; j++ ) {
				if ( G_i[nNodes*i+j] + G_i[nNodes*j+i] > 0 ) {
					/* Supernodes i and j are connected, find out if
					   {i,j} induces a violated capacity constraint */
					VertexSet *S = new VertexSet;

					for ( VertexSet::size_type k = 0; k < superNode[i].size(); k++ ) {
						/* Add primitive nodes of i */
						S->push_back( superNode[i][k] );
					}

					for ( VertexSet::size_type k = 0; k < superNode[j].size(); k++ )
						/* Add primitive nodes of j */
						S->push_back( superNode[j][k] );

					double cutVal = graph.getCutValue( *S );
					double SLB = 2*graph.getSetLb( S );

					if ( SLB - cutVal > eps ) {
						// new cut
						cuts->push_back( S );
					}
					else
						delete S;
				}
			}
		}

		/* Step 2
		   Shrink an edge e not adjacent to the depot such that xe >= 1
		   choose i,j such that thay yield max x(i,j) */
		int best_i = -1, best_j = -1;
		double best_xij = 0.0, xe;
		for ( int i = 1; i < nNodes; i++ ) {
			for ( int j = i+1; j < nNodes; j++ ) {
				if ( ( xe = G_i[nNodes*i+j] + G_i[nNodes*j+i] ) >= 1 && xe > best_xij ) {
					best_i = i;
					best_j = j;
					best_xij = xe;
				}
			}
		}

		if ( best_xij > 0 && best_i != -1 && best_j != -1) {
			/* then surely best_xij >= 1
			   note that best_i < best_j, new supernode will be best_i */

			for ( VertexSet::size_type k = 0; k < superNode[best_j].size(); k++ )
				superNode[best_i].push_back( superNode[best_j][k] );

			/* remove supernode j */
			superNode.erase( superNode.begin() + best_j );

			/* Now update G_i 
			   The new contracted graph has 1 node less */
			delete [] G_i;
			nNodes--;
			G_i = new double[nNodes*nNodes];

			for ( int h = 0; h < nNodes; h++ ) {
				for ( int k = 0; k < nNodes; k++ ) {
					/* calculate h->k contribution to omega(h,k)
					   in the contracted graph */
					double htok = 0.0;
							
					for ( VertexSet::size_type p = 0; p < superNode[h].size(); p++ )
						/* For each p primitive node of h */
						for ( VertexSet::size_type q = 0; q < superNode[k].size(); q++ )
							/* For each q primitive node of k
							   add x(p,q) of the residual graph G */
							htok += G[N*superNode[h][p] + superNode[k][q]];

					G_i[nNodes*h+k] = htok;
				}
			}
		}
		else {
			delete [] G_i;
			G_i = NULL;
			break;
		}
	} // while ( true )

	if ( G_i != NULL )
		delete [] G_i;

	// if new cuts found, delete duplicates
	if ( cuts->size() > nCuts )
		cleanCutSet( cuts );
}

void VrpSep::greedyRoundCap( vector<VertexSet*> *cuts, double eps, bool findAll ) const
{
	vector<VertexSet*>::size_type nCuts = cuts->size();

	int N = graph.N;
	const double *G = graph.G;

	bool *taken = new bool[N];
	for ( int i = 0; i < N; i++ ) {
		taken[i] = false;
	}
	taken[0] = true ;

	for ( int i = 1; i < N && ( cuts->size() < 10 || findAll ); i++ ) {
		// for each client node
		VertexSet S( 1, i );
		taken[i] = true;

		double sCut = 2;

		bool increased;
		do {
			// try to extend S
			increased = false;

			double minCut = DBL_MAX;
			int bestClient = -1;

			// find v such that omega(delta(S+{v})) is minimized
			for ( int j = 1; j < N; j++ ) {
				if ( taken[j] ) // client already in S
					continue;

				double incidentVal = 0.0L;
				for ( VertexSet::size_type h = 0; h < S.size(); h++ ) {
					// sum arcs incident on S and j
					incidentVal += G[N*S[h]+j] + G[N*j+S[h]];
				}

				/* Calculate omega(delta(S U {j}))
				  remove from sCut the amount of value incident to j and add
				  ( 2-that value ), which is the 'net' contribution of j */
				double cutVal = sCut - incidentVal + (2 - incidentVal);

				if ( cutVal < minCut ) {
					minCut = cutVal;
					bestClient = j;
					increased = true;
				}
			}

			if ( increased ) {
				// add new client
				S.push_back( bestClient );
				taken[bestClient] = true;

				sCut = minCut; // update sCut

				// check if now S violates a constraint, in case add it to the cutset
				double SLB = 2*graph.getSetLb( &S );

				if ( SLB - minCut > eps ) {
					cuts->push_back( new VertexSet( S ) );
				}
			}
		} while ( increased && ( cuts->size() < 10 || findAll ) );

		// clear taken marks
		for ( int j = 1; j < N; j++ )
			taken[j] = false;
	}

	delete [] taken;

	// if new cuts found, delete duplicates
	if ( cuts->size() > nCuts )
		cleanCutSet( cuts );
}

/** Generate cuts for the solver si
The cuts are inserted in the collection cs */
void VrpSep::generateCuts( const OsiSolverInterface &si, OsiCuts &cs,
		const CglTreeInfo info ) const
{
	double eps;
	si.getDblParam( OsiPrimalTolerance, eps );

	graph.buildResidualGraph( si.getColSolution() );
	vector<VertexSet*> *cuts = new vector<VertexSet*>;
	
	assert ( numCols == si.getNumCols() );
	//std::cout << "*******************************************\n";

	detachedComponents( cuts );

	if ( cuts->size() == 0 && info.level == 0 ) {
		// apply all heuristics
		ccHeuristic( cuts, eps );
		shrinkingHeuristic( cuts, eps );
		greedyRoundCap( cuts, eps );
	}
	else if ( cuts->size() == 0 ) {
		ccHeuristic( cuts, eps );

		if ( cuts->size() == 0 )
			shrinkingHeuristic( cuts, eps );

		if ( cuts->size() == 0 )
			greedyRoundCap( cuts, eps );
	}

	// no longer useful
	delete [] graph.G;
	graph.G = NULL;

	if ( cuts->size() > 0 ) {
		addCuts( cuts, cs );
	} else if ( useGomory && info.pass < 250 ) {
		// resort to Gomory's cuts if not too many passes have been performed
		// and heuristics have failed to separate the problem
		
		CoinTimer timer;
		timer.restart();
		
		CglGomory gomory;
		
		// a larger limit seems to counter-balance the gains of finding
		// more cuts by losing too much time in each node
		if ( numCols < 500 )
			gomory.setLimit( 50 );
		else
			gomory.setLimit( 100 );

		gomoryTotal++;
		
		OsiCuts gcs;
		gomory.generateCuts( si, gcs, info );
		
		//double violation = 0.0L;
		double best_v = 0.0L;		
		
		if ( gcs.sizeRowCuts() ) {
			for ( int i = 0; i < gcs.sizeRowCuts(); i++ ) {
				double v;
				if ( (v=gcs.rowCutPtr( i )->violated( si.getColSolution() )) > best_v ) {
					best_v = v;
				}
			}
			
			if ( info.level != 0 || best_v > 5.0L ) {
				// add cuts if either at the root with enough violation
				// or in a sub-problem
				cs.insert( gcs );
				gomoryEffective++;
			}
			
		}
		
		gomoryTime += timer.timeElapsed();
		
	}

	// deallocate sets
	for ( vector<VertexSet*>::size_type i = 0; i < cuts->size(); i++ )
		delete (*cuts)[i];
	delete cuts;
}

void VrpSep::addCuts( vector<VertexSet*> *cuts, OsiCuts &cs ) const
{
	vector<VertexSet*>::size_type k = 0;
	VertexSet *S, *CS; /* Set and complement of the set */
	
	/*VertexSet *err = new VertexSet(9);
	VertexSet& strange = *err;
	strange[0] =3;
	strange[1] =8;
	strange[2] =16;
	strange[3] =18;
	strange[4] =19;
	strange[5] =24;
	strange[6] =26;
	strange[7] =27;
	strange[8] =29;	
	cuts->push_back( err );*/

	while ( k < cuts->size() ) {
		S = (*cuts)[k];
		
		/*if ( S->size() == strange.size()) {
			bool equal = true;
			for ( int i = 0; i < strange.size() && equal; i++ )
				if ( (*S)[i] != strange[i] )
					equal = false;
			
			if ( equal ) {
				k++;
				continue;
			}
		}*/
		
		CS = graph.getComplement( S );

		OsiRowCut rc;

		/* We have exactly (N-ScplSize)*ScplSize variables involved in the row
		   For each node in S we consider the acrs joining it to each node in
		   S_ => |S|*|V\S| arcs involved */
		int nVars = (graph.N - CS->size()) * CS->size();

		int *idxs = new int[nVars];
		double *vals = new double[nVars];
		int h = 0; /* Current index added */
		
		/*for ( VertexSet::size_type i = 0; i < S->size(); i++ ) {
			std::cout << (*S)[i] << " ";
		}
		std::cout << std::endl;
		for ( VertexSet::size_type i = 0; i < CS->size(); i++ ) {
			std::cout << (*CS)[i] << " ";
		}
		std::cout << std::endl;*/

		for ( VertexSet::size_type i = 0; i < S->size(); i++ ) {

			for ( VertexSet::size_type j = 0; j < CS->size(); j++ ) {
				assert( (*S)[i] != (*CS)[j] );
				int p = min( (*S)[i], (*CS)[j] ), q = max( (*S)[i], (*CS)[j] );

				int offset = q - p - 1;
				idxs[h] = graph.vBase[p] + offset;
				assert( idxs[h] <= numCols );
				vals[h] = 1.0L;
				h++;
			}
		}

		assert( h == nVars );
		double rowLb = double( 2*graph.getSetLb( S ) );

		// add cut
		rc.setRow( nVars, idxs, vals );
		rc.setLb( rowLb );
		rc.setUb( COIN_DBL_MAX );
		//rc.setUb( double( 2*graph.M ) );
		rc.setGloballyValid( true );

		cs.insert( rc );

		/* Clean up */
		delete CS;
		delete [] idxs;
		delete [] vals;

		k++;
	}
}
