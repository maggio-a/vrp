/************************************************************
  Autore: Andrea Maggiordomo - mggndr89@gmail.com
  Università di Pisa - Informatica Applicata - Logistica
  2012
*************************************************************/

#ifndef VrpSolver_H
#define VrpSolver_H

#include "OsiClpSolverInterface.hpp"
#include "VrpTypes.hpp"

#include <fstream>
#include <vector>

using std::vector;
using std::ofstream;

class VrpSolver {
public:
	/** Constructor */
	VrpSolver( const char *inputFile, bool standardVrp=false );

	/** Destructor */
	~VrpSolver();

	/** Call CbcModel solver to perform branch and cut */
	int cbcSolve( int logLevel=1, bool doCutOff=false, int whichCuts=0, double maxTime=1e100, double fgt=0.0L );

	/** Call Internal algorithm to perform branch and cut - slower */
	int branchAndCut( double gapTol=0.0L, double branchTarget=0.5L, double TimeOut=600.0L );

	/** Cutting planes - processes only the root node
	    if strongBound then solves the combinatorial relaxation strenghtened with
		the cuts found using integer variables */
	int cuttingPlanes( bool doLog=true, bool strongBound=true );

	/** Call Tabu Search procedure - based on some of the idea discussed by
	    Gendrau, Hertz and Laporte in "A Tabu Search heuristic for
		the Vehicle Routing Problem (1994) */
	int tabuRoute( int iFact, bool doLog=true, int seed=12345, bool silent=false, long lim=100000 );

private:
	// Variables
	// -----------------------------------------------------------------------
	int M;                          /* Max number of trucks allowed */
	int Cw;                         /* Capacity (weight) of the single truck */
	int Cv;                         /* Capacity (volume) of the single truck */
	int N;                          /* Dimension ( number of nodes ) */
	int *W;                         /* Clients demand (weight) */
	int *V;                         /* Clients demand (volume) */
	int *C;                         /* Distances matrix */

	int zlb;
	int zub;
	vector< vector<int> > x_bar;

	// MIP solver specific
	OsiClpSolverInterface *solver_; /* Base model solver */
	char **colNames;                /* Variable names to be used with
								       writeLpNative */
	int *vBase;					    /* Base indices for variables*/

	// Misc
	ofstream log;                   /* Log file */

	enum Dimension {
		weightCapacity,
		volumeCapacity
	};

	// Procedures
	// -----------------------------------------------------------------------

	/** Rounding - suggested in the documentation of TSPLIB */
	int nint( double x );

	/** Gap */
	double solutionGap();

	void logSolution();

	/** local search */
	void ts( int tmax, int p, int q, RouteSet &bestSol, RouteSet &bestFsblSol, long lim );

	/** feasibility */
	bool isFeasible( const RouteSet &S );

	/** dimension-specific feasibility */
	bool feasibleLoad( const RouteSet &S, VrpSolver::Dimension type );

	/** cost evaluations */
	int F1( const RouteSet &solution, bool ignoreUnf=false );
	double F2( const RouteSet &solution, double alpha, double beta );

	/** insertion */
	RouteSet putInRoute( int node, int routeIndex, const RouteSet &sol );

	/** trivial optimizations */
	void two_opt( double *tour, int *perm, int &tlen );
	void postOpt( RouteSet &solution );

	/* Comparison operator to sort clients from nearest to farthest
	   from the source client */
	struct Nearest {
		
		Nearest( int o, int n, const int *c ) {
			Nearest::o = o;
			Nearest::n = n;
			Nearest::c = c;
		}

		bool operator()( const int i, const int j ) const {
			return ( c[n*o + i] < c[n*o + j] );
		}

		int o;
		int n;
		const int *c;
	};
};

inline int VrpSolver::nint( double x ) { return (int) ( x + 0.5L ); }

inline double VrpSolver::solutionGap() { return (zlb>zub) ? 0.0L :
	( double(zub-zlb)/abs( ( zlb!=0 ) ? double(zlb) : 1e-50 ) ); }

/* Column bound to define the value at wich x_idx was fixed
   during branching */
struct ColBound {
	
	ColBound() {
		lower = 0.0L;
		upper = 0.0L;
	}

	ColBound( double l, double u ) {
		lower = l;
		upper = u;
	}

	double lower;
	double upper;
};

#endif
