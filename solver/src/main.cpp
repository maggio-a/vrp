/************************************************************
  Autore: Andrea Maggiordomo - mggndr89@gmail.com
  Università di Pisa - Informatica Applicata - Logistica
  2012
*************************************************************/

#include "VrpSolver.hpp"
#include "utils.hpp"

#include <iostream>

using std::cerr;
using std::endl;

void printOptions() {
	
	cerr << "Algorithms" << endl
		<< "  :: Tabu Search meta-heuristic" << endl
		<< "  -h [tMax(int) seed(int)]"<< endl
		<< "     tMax: maximum non-improving iterations (default 50)" << endl
		<< "     seed: (default 12345)" << endl
		<< endl
		<< "  :: Cutting planes (root only)" << endl
		<< "  -p [doLog(0,1) strongBound(0,1)]" << endl
		<< "     doLog: log to file the procedure's outcome (false/true default 1)" << endl
		<< "     strongBound: perform branch and bound on the extended model (false/true default 1)" << endl
		<< endl
		<< "  :: Branch and cut (own)" << endl
		<< "  -b [timeOut(seconds) gapTolerance(float) branchTarget(float) ]" << endl
		<< "     timeout: defines maximum running time (default unlimited)" << endl
		<< "     fractionGapTolerance: interrupt search when relative gap is <= this value (default 0)" << endl
		<< "     branchTarget: branching is performed on variable closest to this fraction (default 0.5=.5)" << endl
		<< endl
		<< "  :: Branch and cut (CbcModel)" << endl
		<< "  -c [logLevel(int) doCutOff(0,1) whichCuts(int) timeout(seconds) fractionGapTolerance(double)] " << endl
		<< "     logLevel: solver's verbosity [ 0 - silent, 1 - reduced, >1 verbose ] (default 1)" << endl
		<< "     doCutOff: impose a heuristic cut off value on the objective function (default 0)" << endl
		<< "     whichCuts: choose solver's separators [ 0 - vrp heuristics, >0 try Gomory's cuts if heuristics fail, "
		         << "2 also use CglProbing ] (default 0)" << endl
		<< "     timeout: set maximum running time (default unlimited)" << endl
		<< "     fractionGapTolerance: interrupt search when relative gap is <= this value (default 0)" << endl;
}


int main( int argc, char **argv )
{
	int tMax = 50;
	int seed = 12345;

	int doLog = 1;
	int strongBound = 1;

	double frGapTol = 0.0L;
	double bt = 0.5L;
	double to = 1e100;
	int logLevel = 1;
	bool doCutOff = false;
	int whichCuts = 0;

	if ( argc < 3 ) {
		cerr << "Usage: " << argv[0] << " filename algorithm" << endl;
		printOptions(); 
		exit( EXIT_FAILURE );
	}

	std::string arg = argv[2];

	if ( arg == "-h" ) {
		switch ( argc ) {
		case 5:
			get_arg( argv[4], seed );
		case 4:
			get_arg( argv[3], tMax );
		}
	}
	else if ( arg == "-p" ) {
		switch ( argc ) {
		case 5:
			get_arg( argv[4], strongBound );
		case 4:
			get_arg( argv[3], doLog );
		}
	}
	else if ( arg == "-b" ) {
		switch( argc ) {
		case 6:
			get_arg( argv[5], bt );
		case 5:
			get_arg( argv[4], frGapTol );
		case 4:
			get_arg( argv[3], to );
		}
	}
	else if ( arg == "-c" ) {
		switch ( argc ) {
		case 8:
			get_arg( argv[7], frGapTol );
		case 7:
			get_arg( argv[6], to );
		case 6:
			get_arg( argv[5], whichCuts );
		case 5:
			get_arg( argv[4], doCutOff );
		case 4:
			get_arg( argv[3], logLevel );
		}
	}
	else {
		cerr << "Usage: " << argv[0] << " filename algorithm" << endl;
		printOptions(); 
		exit( EXIT_FAILURE );
	}

	VrpSolver problem( argv[1] );

	if ( arg == "-h" )
		problem.tabuRoute( tMax, true, seed );
	else if ( arg == "-p" )
		problem.cuttingPlanes( doLog != 0, strongBound != 0 );
	else if ( arg == "-b" )
		problem.branchAndCut( frGapTol, bt, to );
	else if ( arg == "-c" )
		problem.cbcSolve( logLevel, doCutOff, whichCuts, to, frGapTol );

	return 0;
}
