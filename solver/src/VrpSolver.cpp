/************************************************************
  Autore: Andrea Maggiordomo - mggndr89@gmail.com
  Università di Pisa - Informatica Applicata - Logistica
  2012
*************************************************************/

#include "VrpSolver.hpp"

#include "CoinModel.hpp"
#include "OsiClpSolverInterface.hpp"

#include "VrpTypes.hpp"
#include "VrpHeuristicSeparator.hpp"
#include "VrpData.hpp"

#include "utils.hpp"
#include "CoinTime.hpp"

#include "CbcModel.hpp"
#include "OsiAuxInfo.hpp" // OsiBabSolver
#include "CbcCutGenerator.hpp"
#include "CbcFeasibilityNoStrong.hpp"

#include "CglProbing.hpp"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <queue>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <climits>

using std::cout;
using std::cerr;
using std::endl;
using std::setprecision;
using std::ifstream;
using std::ostringstream;
using std::ios_base;
using std::deque;
using std::abs;
using std::min;
using std::max;
using std::sort;
using std::set;
using std::swap;

/** Constructor */
VrpSolver::VrpSolver( const char *inputFile, bool standardVrp )
{

	/* Read instance */
	ifstream instance( inputFile );

	/* Problem constraints */
	if ( !instance.is_open() ) {
		cerr << "Cannot open file " << inputFile << "\n";
		exit( EXIT_FAILURE );
	}

	instance >> M;
	if ( instance.fail() ) {
		cerr << "Error reading M from instance file\n";
		exit( EXIT_FAILURE );
	}

	if ( M <= 0 ) {
		cerr << "The number of trucks available must be positive\n";
		exit( EXIT_FAILURE );
	}

	instance >> Cw;
	if ( instance.fail() ) {
		cerr << "Error reading Cw from instance file\n";
		exit( EXIT_FAILURE );
	}

	if ( Cw <= 0 ) {
		cerr << "Weight capacity of the trucks must be positive\n";
		exit( EXIT_FAILURE );
	}

	instance >> Cv;
	if ( instance.fail() ) {
		cerr << "Error reading Cv from instance file\n";
		exit( EXIT_FAILURE );
	}

	if ( Cv <= 0 ) {
		cerr << "Volume capacity of the trucks must be positive\n";
		exit( EXIT_FAILURE );
	}

	/* Nodes info */
	instance >> N;
	if ( instance.fail() ) {
		cerr << "Error reading N from instance file\n";
		exit( EXIT_FAILURE );
	}

	if ( N <= 0 ) {
		cerr << "No clients specified by the instance\n";
		exit( EXIT_FAILURE );
	}

	/* Read coordinates */
	double *xpos = new double[N];
	double *ypos = new double[N];

	int id; /* Useless */

	for ( int i = 0; i < N; i++ ) {
		instance >> id; /* Reads the node id */
		if ( instance.fail() ) {
			cerr << "Error reading node info from instance file\n";
			exit( EXIT_FAILURE );
		}

		instance >> xpos[i];
		if ( instance.fail() ) {
			cerr << "Error reading node info from instance file\n";
			exit( EXIT_FAILURE );
		}

		instance >> ypos[i];
		if ( instance.fail() ) {
			cerr << "Error reading node info from instance file\n";
			exit( EXIT_FAILURE );
		}
	}

	/* Now generate the distances matrix
	   C[(N*i) + j] specifies the distance from i to j 
	   note that because the distances are euclidean the matrix
	   is actually simmetric and so C[(N*i)+j] = C[(N*j)+i] for
	   any given i and j */
	C = new int[N*N];

	double xd, yd;
	for ( int i = 0; i < N; i++ ) {
		for ( int j = 0; j < N; j++ ) {
			xd = xpos[i] - xpos[j];
			yd = ypos[i] - ypos[j];

			C[(N*i) + j] = nint( sqrt( xd*xd + yd*yd ) );
		}
	}

	delete [] xpos;
	delete [] ypos;

	/* Read clients demand */
	W = new int[N];
	V = new int[N];
	for ( int i = 0; i < N; i++ ) {
		instance >> id; /* Read the node id */
		if ( instance.fail() ) {
			cerr << "Error reading demand info from instance file\n";
			exit( EXIT_FAILURE );
		}

		instance >> W[i];
		if ( instance.fail() ) {
			cerr << "Error reading demand info from instance file\n";
			exit( EXIT_FAILURE );
		}

		int v;
		instance >> v;
		if ( instance.fail() ) {
			cerr << "Error reading demand info from instance file\n";
			exit( EXIT_FAILURE );
		}
		V[i] = standardVrp ? 0 : v;
	}

	/* Create indices to easily access the variables
	   Since we operate with an undirected graph, we associate a variable to
	   each edge - suppose we have 6 edges, then we associate:

	     x0 -> edge 0-1 => vBase[0] = 0
		 x1 -> edge 0-2
		 x2 -> edge 0-3
		 x3 -> edge 1-2 => vBase[1] = 3
		 x4 -> edge 1-3
		 x5 -> edge 2-3 => vBase[2] = 5

	   so its easy to retrieve the variable associated with a given edge */
	int nVars = (N*(N-1))/2;
	vBase = new int[N-1];
	vBase[0] = 0;
	for ( int i = 1; i < N-1; i++ )
		vBase[i] = vBase[i-1] + N - 1 - (i-1);

	assert( vBase[N-2] == nVars-1 );

	/* Now build the model
	   The approach used here is to solve a combinatorial relaxation of the
	   problem (in which we disregard the capacity constraints) with continuous
	   variables. What we do is then find some violated
	   constraints, add them to the problem, and resolve with the new
	   constraints added in a branch and cut fashion, loking for feasible
	   integer solutions along the way.
	   The master model, stored in solver_, will hold the initial relaxation */

	CoinModel *model = new CoinModel;

	model->setOptimizationDirection( 1.0 ); /* Minimize */
	
	/* sum( j in STAR(0), x_0j ) = 2*M

	   The first N-2 represent arcs inciding on the depot, the sum of
	   these variables must be equal to 2*M */
	int *idxs = new int[N-1];
	double *vals = new double[N-1];
	for ( int i = 0; i < N-1; i++ ) {
		idxs[i] = i;
		vals[i] = 1.0L;
	}

	model->addRow( N-1, idxs, vals, double( 2*M ), double (2*M ) );
	
	/* For each customer, impose a degree of 2 on the incident arcs
	   in solution */
	for ( int customer = 1; customer < N; customer++ ) {

		int i = 0;
		for ( int k = 0; k < customer; k++, i++ ) {
			int offset = customer - k - 1;
			idxs[i] = vBase[k] + offset;
			vals[i] = 1.0L;
		}

		if ( customer < N-1 ) {
			for ( int k = vBase[customer]; k < vBase[customer] + N - 1 - customer; k++, i++ ) {
				idxs[i] = k;
				vals[i] = 1.0L;
			}
		}

		model->addRow( N-1, idxs, vals, 2.0L, 2.0L );
	}

	delete [] idxs;
	delete [] vals;

	colNames = new char*[nVars];

	/* Set column constraints */
	for ( int i = 0; i < N-1; i++ ) {
		for( int j = i+1; j < N; j++ ) {
			// edges i-j
			int offset = j - i - 1;
			int colId = vBase[i] + offset;
			model->setColumnBounds( colId, 0.0L, ( i == 0 ) ? 2.0L : 1.0L );
			model->setColumnObjective( colId, double( C[N*i+j] ) );
			model->setInteger( colId );

			// naming
			ostringstream os;
			os << "x(" << i << "," << j << ")";
			colNames[colId] = os2ch( os );
		}
	}

	solver_ = new OsiClpSolverInterface;
	solver_->loadFromCoinModel( *model );

	delete model;

	// Debug print 
	// solver->writeLpNative( tmp, NULL, colNames );

	zlb = - INT_MAX;
	zub = INT_MAX;

	// open log file
	ostringstream logstream;
	log.open( "solution.log", ios_base::out|ios_base::trunc );
	log << ":: " << inputFile << " ::" << endl;

}

/** Destructor */
VrpSolver::~VrpSolver()
{
	delete solver_;
	delete [] W;
	delete [] V;
	delete [] C;
	delete [] vBase;

	int nVars = (N*(N-1))/2;
	for ( int i = 0; i < nVars; i++ )
		delete [] colNames[i];
	delete [] colNames;

}

/** Branch and cut provided by Cbc */
int VrpSolver::cbcSolve( int logLevel, bool doCutOff, int whichCuts, double maxTime, double fgt )
{
	CoinTimer wallclock;
	// Load solver in the model
	CbcModel model( *solver_ );
	log << endl << "Cbc branch and cut" << endl;
	log << "- doCutOff " << doCutOff << endl;
	log << "- whichCuts " << whichCuts << endl;
	log << "- maxTime " << maxTime << endl;
	log << "- gap tolerance " << fgt*100.0L << "%" << endl;

	VrpData G( M, N, Cw, Cv, W, V, vBase );
	VrpSep separator( G, (whichCuts != 0) );
	
	if ( whichCuts == 2 ) {
		CglProbing generator1;
		model.addCutGenerator( &generator1, -1, "Probing" );
	}

	model.addCutGenerator( &separator, 1, "VrpHeuristics", true, false, false, -100, 1, -1 );
	// get statistics
	int numberGenerators = model.numberCutGenerators();
	int iGenerator;
	for ( iGenerator = 0; iGenerator < numberGenerators; iGenerator++ ) {
    	CbcCutGenerator *generator = model.cutGenerator( iGenerator );
    	generator->setTiming( true );
	}
	
	// First continuous solve
	model.initialSolve();

	/* Before we start the branch and cut we need to tell the solver
	 - it has to call the cut generator as long as cuts are found
	 - an integer solution may be cut out by the cut generator
	 - to ignore solutions found after strong branching as those are not passed
	   through the cut generator and might therefore be integer but unfeasible */

	// Subsequent calls of the generator
	iGenerator = numberGenerators-1;
	model.cutGenerator( iGenerator )->setMustCallAgain(true);
	//model.cutGenerator( 0 )->setMustCallAgain( true );

	// Try to cut integer solutions
	OsiBabSolver solverCharacteristics;
	solverCharacteristics.setSolverType( 4 );
	model.solver()->setAuxiliaryInfo( &solverCharacteristics );

	// Filter solutions found during strong branching
	CbcFeasibilityNoStrong ns( G );
	model.setProblemFeasibility( ns );

	// Also say don't recompute solution - example allCuts.cpp
	// I think this prevents double checks when continuous integer...
	model.setSpecialOptions( 4 );

	// Parameters
	if ( logLevel == 0 )
		model.setLogLevel( logLevel );
	else if ( logLevel == 1 )
		model.solver()->setHintParam( OsiDoReducePrint, true );
	else
		model.setLogLevel( 2 );

	model.setNumberStrong( 20 );
	
	double cutoffTime = 0.0L;
 
	if ( doCutOff ) {
		CoinTimer timer;
		// set cutoff to heuristic solution value + 1 to avoid cutting the optimal value
		// solve multiple times with different seeds to shake things up
		cout << "Calculating cut off on objective value" << std::flush;
		for ( int seed = 0; seed < 5; ++seed ) {
			cout << "." << std::flush;
			tabuRoute( 30, false, seed, true, 5000 ); //suppress all messages and limit iterations
		}
		cout << endl;
		double cutoff = double( ( zub == INT_MAX ) ? zub : zub+1 );
		model.setCutoff( cutoff );
		cutoffTime = timer.timeElapsed();
		cout << ":: Setting cut off on objective value to " << cutoff << endl << endl;
		log << "Setting cut off on objective value to " << cutoff << endl << endl;
	}

	model.setMaximumSeconds( maxTime );
	model.setAllowableFractionGap( fgt );

	// start search
	try {
		model.branchAndBound();
	} catch( CoinError &e ) {
		e.print();
		exit( EXIT_FAILURE );
	}

	cout << endl;

	if ( model.isAbandoned() ) {
		cerr << endl << "Search was abandoned!" << endl;
		exit( EXIT_FAILURE );
	}
	else {
		bool unf = false;
		const double *best_sol = NULL;
		if ( model.isProvenInfeasible() ) {
			cout << "Problem proven to be unfeasible" << endl;
			log << "Problem proven to be unfeasible" << endl;
			unf = true;
		}
		else if ( model.isProvenOptimal() ) {
			cout << endl << "Problem solved to optimality" << endl;
			log << endl << "Problem solved to optimality" << endl << endl;
			best_sol = model.getCbcColSolution();
		}
		else {
			cout << endl << "Search interrupted, retrieving best solution found" << endl;
			log << endl << "Search interrupted, retrieving best solution found" << endl;
			best_sol = model.bestSolution();
		}
		
		if ( !unf && best_sol != NULL ) {
			
			// optimal or interrupted - output a solution if available
			//G.buildResidualGraph( model.getCbcColSolution() );
			G.buildResidualGraph( best_sol );
			RouteSet r = G.getRoutes( model.getDblParam( CbcModel::CbcIntegerTolerance ) );
			//delete [] sol;

			if ( G.isFeasible( r ) ) {

				cout << "========================================" << endl;
				for ( RouteSet::size_type i = 0; i < r.size(); i++ ) {
					cout << "Route #" << i << ": ";
					log << "Route #" << i << ": ";
					int w = 0, v = 0;
					for ( Route::size_type j = 0; j < r[i].size(); j++ ) {
						cout << r[i][j] << " ";
						log << r[i][j] << " ";
						v += V[ r[i][j] ];
						w += W[ r[i][j] ];
					}
					cout << "[ W=" << w << ", V=" << v <<" ]" << endl;
					log << "[ W=" << w << ", V=" << v <<" ]" << endl;
				}
				cout << "Cost " << model.getObjValue() << endl;
				log << "Cost " << model.getObjValue() << endl;
				if ( model.status() == 1 || model.secondaryStatus() == 2 ) { // search interrupted
					double gap = (model.getObjValue() - model.getBestPossibleObjValue())/model.getBestPossibleObjValue();
					cout << "Best possible " << model.getBestPossibleObjValue() << ", gap " << gap*100 << "%" << endl;
					log << "Best possible " << model.getBestPossibleObjValue() << ", gap " << gap*100 << "%" << endl;
				}
				cout << "========================================" << endl;
			}
		} // if !unf
		else if ( !unf ) {
			// no fesible solution available
			cout << "No feasible solution available" << endl;
			log << "No feasible solution available" << endl;
		}
		
		if ( cutoffTime ) {
			cout << endl << "Cutoff stage took " << cutoffTime << " seconds" << endl;
			log << endl << "Cutoff stage took " << cutoffTime << " seconds" << endl;
		}
		
		log << endl << "Execution time: " << model.getCurrentSeconds() << " seconds" << endl;
		log << model.getNodeCount() << " nodes processed" << endl;
	}
	
	cout << endl;
	
	log << endl << "Wallclock: " << wallclock.timeElapsed() << " seconds." << endl << endl;
	
	cout << "Cuts at root node changed objective from " << model.getContinuousObjective()
		<< " to " << model.rootObjectiveAfterCuts() << endl;
	
	for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
		CbcCutGenerator *generator = model.cutGenerator( iGenerator );
		cout << generator->cutGeneratorName() << " - Called "
			<< generator->numberTimesEntered() << " times, "
			<< generator->numberCutsInTotal() << " row cuts and "
			<< generator->numberColumnCuts() << " column cuts ";
		if ( generator->timing() )
			cout << "(" << generator->timeInCutGenerator() << " seconds)" << endl;
	    else
			cout << endl;
	}
	
	if ( whichCuts ) {
		iGenerator = numberGenerators-1;
		VrpSep *sep = dynamic_cast<VrpSep*>( model.cutGenerator( iGenerator )->generator() );
	
		cout << "Separation resorted to Gomory's cuts "<< sep->gomoryTotal
			<< " times, finding additional cuts " << sep->gomoryEffective << " times("
			<< sep->gomoryTime << " seconds)"<< endl;
	}

	return 0;
}

/** Branch and cut - own implementation */ 
int VrpSolver::branchAndCut( double gapTol, double branchTarget, double timeOut )
{
	if ( branchTarget < 0.0L || branchTarget > 1.0L ) {
		cerr << "Error: branch target must be in [0,1]" << endl;
		exit( EXIT_FAILURE );
	}

	log << "== Branch and cut ";
	for ( int i = 0; i < 61; i++ )
		log << "=";
	log << endl;

	log << "Gap tolerance " << gapTol << endl;
	log << "Branch on the variable closest to: " << branchTarget << endl;
	log << "Time limit: " << timeOut << " seconds" << endl;
	for ( int i = 0; i < 79; i++ )
		log << "=";
	log << endl << endl;

	// call taburoute to find a valid upper bound zub
	cout << "Calculating upper bound on objective value" << std::flush;
	for ( int seed = 0; seed < 5; ++seed ) {
		cout << "." << std::flush;
		tabuRoute( 30, false, seed, true, 5000 ); //suppress all messages and limit iterations
	}

	if ( zub == INT_MAX ) {
		cout << "Notice: starting branch and cut procedure with no upper bound" << endl;
		log << "Notice: starting branch and cut procedure with no upper bound" << endl;

		int trivial_zub = 0; // all arcs in solution
		for ( int i = 0; i < N; i++ )
			for ( int j = 0; j < N; j++ )
				trivial_zub += C[N*i+j];
		zub = trivial_zub;
	}

	// build solver from base model
	OsiClpSolverInterface *solver = new OsiClpSolverInterface( *solver_ );

	double eps;
	solver->getDblParam( OsiPrimalTolerance, eps );

	int ctot = 0;

	CoinTimer timer;
	timer.restart();

	logSolution();

	// load separator
	VrpData info( M, N, Cw, Cv, W, V, vBase );
	VrpSep separator( info );

	// a tree node is just a vector of bounds on variables
	typedef vector<ColBound> TreeNode;

	// each sub-polyhedron is defined by a list of column bounds
	// use a deque for the tree levels
	deque<TreeNode> parentLevel, childLevel; // levels of the tree

	// push first node, no restrictions
	TreeNode root( solver->getNumCols() );
	const double *colLower = solver->getColLower();
	const double *colUpper = solver->getColUpper();
	for ( int i = 0; i < solver->getNumCols(); i++ ) {
		ColBound b( colLower[i], colUpper[i] );
		root[i] = b;
	}
	// assert( (N*(N-1))/2 == root.size() );
	// push root node in the parent level
	parentLevel.push_front( root );

	int levelLowerBound = INT_MAX;
	
	// free all variables
	for ( int i = 0; i < solver->getNumCols(); i++ )
		solver->setColBounds( i, 0.0L, ( i < N-1 ) ? 2.0L : 1.0L );
	
	cout << "Performing initial solve" << endl;
	
	// initial solve
	solver->initialSolve();

	if ( !solver->isProvenOptimal() ) {
		log << "Preliminary solve failed!" << endl;
		cout << "Preliminary solve failed!" << endl;
		exit( EXIT_FAILURE );
	}
	else {
		// possibly update the lower bound
		int lb = int( solver->getObjValue() );
		lb = ( abs(lb - int(lb)) <= eps ) ? lb : lb+1; // rounds val to the next integer
		if ( lb > zlb ) {
			log << "Updating lower bound: " << zlb << " -> " << lb << endl;
			zlb = lb;
		}
	}
	
	solver->setHintParam( OsiDoReducePrint, true, OsiHintDo );
	
	cout << endl << "Branch and cut started..." << endl;
	int it1 = 0;
	// start branch and cut procedure
	do {
		it1++;
		
		if ( it1%500 == 0 ) {
			// print some info to reassure the user :>
			cout << "After " << it1 << " nodes processed - best objective value is "
			     << zub << ", best possible is " << zlb << ", gap = " << solutionGap()
			     << endl; 
			cout << "Procedure was started " << timer.timeElapsed()
			     << " seconds ago" << endl;
		}
	
		if ( parentLevel.empty() ) {
			// done sweeping current tree level, child level becomes parent level
			parentLevel = childLevel;

			// update lower bound
			if ( levelLowerBound > zlb ) {
				log << "Updating lower bound: " << zlb << " -> "
					<< levelLowerBound << endl;
				zlb = levelLowerBound;
			}
			// reset level lower bound
			levelLowerBound = INT_MAX;
			// new child level is empty
			childLevel.clear();
		}

		// check halting conditions
		double gap = solutionGap();
		if (  gap == 0 ) {
			log << endl << "Optimal solution found, stopping" << endl ;
			cout << endl << "Optimal solution found, stopping" << endl ;
			break;
		}
		else if ( gap <= gapTol ) {
			log << endl << "Solution gap less than or equal to gap tolerance, stopping" << endl;
			cout << endl << "Solution gap less than or equal to gap tolerance, stopping" << endl;
			break;
		}
		else if ( timer.timeElapsed() >= timeOut ) {
			log << "TimeOut reached" << endl;
			cout << "TimeOut reached" << endl;
			break;
		}

		// fetch next node
		TreeNode Xc = parentLevel.front();
		parentLevel.pop_front();
	
		// bound variables according to Xc
		for ( TreeNode::size_type i = 0; i < Xc.size(); i++ )
			solver->setColBounds( i, Xc[i].lower, Xc[i].upper );
	
		// start cutting planes
		bool pruned;
		int it2 = 0;

		do {
			it2++;
	
			pruned = false;

			// solve
			solver->resolve();
	
			if ( solver->isProvenPrimalInfeasible() ) {
				// subproblem is empty, prune
				//log << "Node pruned - unfeasible" << endl;
				pruned = true;
			}
			else if ( solver->isProvenOptimal() ) {
	
				// the problem is solved
				const double *x = solver->getColSolution();
				int subProbLb = int( solver->getObjValue() );
				subProbLb = ( abs(subProbLb - int(subProbLb)) <= eps ) ? subProbLb : subProbLb+1; // rounds val to the next integer
	

				if ( subProbLb >= zub ) {
					// lower bound of the subtree higher than the current
					// solution value, no business in searching here
					pruned = true;
					//log << "Node pruned - upper bound " << zub << ", lb "
					//	<< solver->getObjValue() << endl;
				}
				else {
					// try to cut
					// if inequalities are found apply them and resolve
					// if not, check if the solution is feasible and eventually take it
					OsiCuts cs;
					if ( it1 == 1 ) {
						// inform the generator that this is the root node
						CglTreeInfo root;
						root.level = 1;
						separator.generateCuts( *solver, cs, root );
					}
					else
						separator.generateCuts( *solver, cs );

					if ( cs.sizeRowCuts() == 0 ) {
						// no cuts found, if the solution is feasible we take it and prune the node
						// otherwise we have to branch
						bool integer = true;
	
						for ( int i = 0; i < solver->getNumCols() && integer; i++ ) {
							if ( ( x[i] > 0.0L && x[i] < 1.0L ) || ( x[i] > 1.0L && x[i] < 2.0L ) )
								integer = false;
						}
		
						if ( integer ) {
							// we know solution is connected since no cuts were found
							VrpData solGraph = info;
							solGraph.buildResidualGraph( x );
							RouteSet routes = solGraph.getRoutes();
							// check if it is feasible
							if ( solGraph.isFeasible( routes ) ) {
								// update best solution
								x_bar = routes;
								zub = subProbLb; //F1( routes ); //solver->getObjValue();
								pruned = true;
								log << "Node pruned - optimal solution found" << endl;
								logSolution();
								log << "[Time elapsed: " << timer.timeElapsed() << "]" << endl;
								
								cout << "Node " << it1 <<": found an integer solution of value "
								     << subProbLb << endl;
								
							}
						}

						// if not pruned we failed to find cuts and we have to branch so break
						if ( !pruned ) {
							//cout << "Interrupting the cutting plane iterations (2)" << endl;
							break;
						}
					}
					else {
						// apply cuts to the solver
						solver->applyCuts( cs );
						//for ( int i = 0; i < cs.sizeRowCuts(); i++ ) {
						//	OsiRowCut *rc = cs.rowCutPtr( i );
						//	solver->addRow( rc->row(), rc->lb(), rc->ub() );
						//}
						ctot += cs.sizeRowCuts();
					}
				}
			} // isProvenOptimal
			else {
				cerr << "Something bad happened, exiting" << endl;
				exit( EXIT_FAILURE );
			}
	
		} while ( !pruned );

		// bounding procedure done, possibly update level lower bound
		int lb = int( solver->getObjValue() );
		lb = ( abs(lb - int(lb)) <= eps ) ? lb : lb+1; // rounds val to the next integer
		levelLowerBound = min( levelLowerBound, lb );
	
		if ( !pruned ) {
			// Branch on Xi
			int min_id;
			double min_slack = 1.0L;
			bool found = false;
			const double *x = solver->getColSolution();

			// find the "most fractional" variable
			for ( int i = 0; i < solver->getNumCols(); i++ ) {
				double target = ( x[i] <= 1 ) ? branchTarget : 1 + branchTarget;
				double slack = abs( target - x[i] );
				if ( slack < min_slack ) {
					found = true;
					min_slack = slack;
					min_id = i;
				}
			}
	
			if ( found ) {
				double i_ub = floor( x[min_id] );
				double i_lb = ceil( x[min_id] );
				// new nodes
				TreeNode n1 = Xc;
				n1[min_id].upper = i_ub; // x[min_id] <= floor
	
				TreeNode n2 = Xc;
				n2[min_id].lower = i_lb; //x[min_id] > floor

				// childLevel.push_front( n1 );
				// childLevel.push_front( n2 );

				childLevel.push_back( n1 );
				childLevel.push_back( n2 );

				/*int lb = int( solver->getObjValue() );
				lb = ( abs(lb - int(lb)) <= eps ) ? lb : lb+1; // rounds val to the next integer
				levelLowerBound = min( levelLowerBound, lb );*/
			}
			else {
				cerr << "Error: could not branch on Xc" << endl;
				log << "Error: could not branch on Xc" << endl;
			}
		}
	
	} while ( !parentLevel.empty() || !childLevel.empty() );

	bool opt = false;
	if ( parentLevel.empty() && childLevel.empty() ) {
		log << endl << " ** Enumeration completed ** " << endl;
		cout << endl << " ** Enumeration completed ** " << endl;
		opt = true;
	}

	// print info
	if ( !x_bar.empty() ) {
		cout << "========================================" << endl;
		for ( RouteSet::size_type i = 0; i < x_bar.size(); i++ ) {
			cout << "Route #" << i << ": ";
			int w = 0, v = 0;
			for ( Route::size_type j = 0; j < x_bar[i].size(); j++ ) {
				cout << x_bar[i][j] << " ";
				v += V[ x_bar[i][j] ];
				w += W[ x_bar[i][j] ];
			}
			cout << "[ W=" << w << ", V=" << v <<" ]" << endl;
		}
		cout << "Cost " << zub << ", gap " << ( (opt) ? 0.0L : solutionGap() ) << endl;
		cout << "========================================" << endl;
	}
	cout << endl << "[Time elapsed: " << timer.timeElapsed() << "]" << endl;
	cout << "Stopped after " << it1 << " nodes processed and " << ctot
	     << " RCIs found" << endl;
	
	
	logSolution();
	log << endl << "[Time elapsed: " << timer.timeElapsed() << "]" << endl;
	log << "Stopped after " << it1 << " nodes processed and " << ctot
	    << " RCIs found" << endl;

	delete solver;


	return 0;
}

/** Log solution to file */
void VrpSolver::logSolution() {

	log << "========================================" << endl;

	log << "Best solution" << endl;
	if ( x_bar.empty() ) 
		log << "No feasible solution found" << endl;
	else {
		for ( RouteSet::size_type i = 0; i < x_bar.size(); i++ ) {
			log << "Route #" << i << ": ";
			int w = 0, v = 0;
			for ( Route::size_type j = 0; j < x_bar[i].size(); j++ ) {
				log << x_bar[i][j] << " ";
				v += V[ x_bar[i][j] ];
				w += W[ x_bar[i][j] ];
			}
			log << "[ W=" << w << ", V=" << v <<" ]" << endl;
		}
		log << "Cost " << zub << ", gap " << solutionGap() << endl;
	}

	log << "========================================" << endl;

}

/** Cutting planes - processing done only at root */
int VrpSolver::cuttingPlanes( bool doLog, bool strongBound )
{	
	OsiClpSolverInterface *solver = new OsiClpSolverInterface( *solver_ );

	double eps;
	solver->getDblParam( OsiPrimalTolerance, eps );

	if ( doLog ) {
		log << "== Cutting planes ";
		for ( int i = 0; i < 61; i++ )
			log << "=";
		log << endl;
		for ( int i = 0; i < 79; i++ )
			log << "=";
		log << endl;
	}

	/* Solve initial problem */
	cout << "Initial solve:" << endl;
	solver->initialSolve();

	if ( solver->isProvenOptimal() ) {
		cout << "Objective value: " << setprecision( 8 )
		     << solver->getObjValue() << endl;

		int val = int( solver->getObjValue() );
		val = ( abs(val - int(val)) <= eps ) ? val : val+1; // rounds val to the next integer

		zlb = max( zlb, val );
	}
	else
		cout << "Could not solve the problem" << endl;

	VrpData info( M, N, Cw, Cv, W, V, vBase );
	VrpSep separator( info );

	/* Now start the row generation
	   The solution is a non directed graph ( basically a disjointed vertex
	   cycle cover ) on which we operate to find as many violated capacity
	   inequalities as possible. The procedure returns a set of client
	   subsets from which we derive the inequalities to add to the problem */

	// trick the separator into thinking the node is always the root
	// so it applies all heuristics
	CglTreeInfo trick;
	trick.level = 1;

	int nIter = 0;
	//bool sameObj;
	//double legacyObj = solver->getObjValue();
	
	solver->setHintParam( OsiDoReducePrint, true, OsiHintDo );
	do {
		cout << "Iteration " << nIter << ":" << endl;
		if ( doLog )
			log << "Iteration " << nIter << ":" << endl;

		OsiCuts cs;
		separator.generateCuts( *solver, cs, trick );

		// no cuts found - break
		if ( cs.sizeRowCuts() == 0 ) {
			cout << "Could not generate another row, " 
				 << "interrupting the cutting planes iterations" << endl;

			if ( doLog )
				log << "Could not generate another row, " 
					<< "interrupting the cutting planes iterations" << endl;
			break;
		}
		
		// apply cuts and resolve
		solver->applyCuts( cs );

		/* Re-solve the problem with the new constraints */
		solver->resolve();

		/* Print obj value */
		if ( solver->isProvenOptimal() ) {

			int val = int( solver->getObjValue() );
			val = ( abs(val - int(val)) <= eps ) ? val : val+1; // rounds val to the next integer

			zlb = max( zlb, val );

			cout << "Objective value " << setprecision( 8 )
				 << solver->getObjValue() << endl;
		
			if ( doLog ) {
				log << "Objective value " << setprecision( 8 )
					<< solver->getObjValue() << endl;

				log << "Number of constraints " << solver->getNumRows() << endl;
			}
		}
		else {
			cout << "Could not solve the problem, exiting" << endl;
			if ( doLog )
				log << "Could not solve the problem, exiting" << endl;
			exit( EXIT_FAILURE );
		}

		//sameObj = abs( legacyObj - solver->getObjValue() ) < eps;
		//legacyObj = solver->getObjValue();

	} while ( nIter++ < 100 /*&& !sameObj*/ ); // until iter cap or obj function soen't change

	if ( nIter == 100 )
		cout << "Cutting plane: iterations limit reached" << endl;

	// if "strong" bound is request, solve with integer variables
	if ( strongBound ) {
		cout << "Starting branch and bound on the extended model" << endl;

		solver->branchAndBound();

		if ( solver->isProvenOptimal() ) {
			cout << "B&B Objective value " << setprecision( 8 )
				 << solver->getObjValue() << endl;
			if ( doLog )
				log << "B&B Objective value " << setprecision( 8 )
					<< solver->getObjValue() << endl;

			int val = int( solver->getObjValue() );
			val = ( abs(val - int(val)) <= eps ) ? val : val+1; // rounds val to the next integer

			zlb = max( zlb, val );

		}
		else {
			cout << "Could not solve the problem, exiting" << endl;
			if ( doLog )
				log << "Could not solve the problem, exiting" << endl;
			exit( EXIT_FAILURE );
		}
	}

	delete solver;

	return 0;
}

/** Tabu search meta-heuristic */
int VrpSolver::tabuRoute( int iFact, bool doLog, int seed, bool silent, long lim )
{
	srand( seed );
	
	if ( silent )
		doLog = false;
	else
		lim = lim*5; // raise limit

	if ( doLog ) {
		log << "== TabuRoute Heuristic ";           
		for ( int i = 0; i < 56; i++ )
			log << "=";
		log << endl;
		for ( int i = 0; i < 79; i++ )
			log << "=";
		log << endl << endl;
	}

	RouteSet bestFsblSol, bestSol( M );

	/* Initialization phase, try to find a decent solution */
	CoinTimer timer;

	bool *taken = new bool[N]; // marks
	int *ind = new int[N]; // sorting
	double *tour = new double[N*N]; // tsp based heuristic

	// Trivial starting solution: closest node next
	for ( int i = 0; i < N; i++ ) {
		taken[i] = false;
		ind[i] = i;
	}
	taken[0] = true;

	int lastNode;
	bool added;
	for ( int i = 0; i < M-1; i++ ) {
		lastNode = 0; // start from depot
		int cleft = Cw; // residual capacity = Q
		int vleft = Cv;
		do {
			added = false;
			sort(ind, ind+N, Nearest( lastNode, N, C ) );

			for ( int j = 0; j < N && !added; j++ ) {
				if ( !taken[ ind[j] ] && W[ ind[j] ] <= cleft && V[ ind[j] ] <= vleft ) {
					added = true;

					// update marks
					taken[ ind[j] ] = true;
					// update cleft vleft
					cleft -= W[ ind[j] ];
					vleft -= V[ ind[j] ];
					// add ind[j] to solution
					bestSol[i].push_back( ind[j] );
					// update last node added
					lastNode = ind[j];
				}
			}
		} while ( added );
	}

	// the last route takes all that is left disregarding capacity,
	// this may lead to an unfeasible solution
	lastNode = 0;
	do {
		added = false;
		sort( ind, ind+N, Nearest( lastNode, N, C ) );

		for ( int j = 0; j < N; j++ ) {
			if ( !taken[ ind[j] ] ) {
				added = true;
				taken[ ind[j] ] = true;
				bestSol[M-1].push_back( ind[j] );
				lastNode = ind[j];
			}
		}
	} while ( added );

	if ( isFeasible( bestSol ) ) {
		postOpt( bestSol );
		bestFsblSol = bestSol;
	}

	if ( doLog ) {
		log << "F1 = " << F1( bestSol ) << " ";
		log << "F2(1.0,1.0) = "<< F2( bestSol, 1.0, 1.0 ) << endl << endl;
	}


	// trivial solution 2: based on a 2-optimized NearestNeighbour tsp
	// pick lambda different starting nodes and keep the best solutions
	bool *picked = new bool[N]; // remember starting nodes already picked
	for ( int i = 0; i < N; i++ )
		picked[i] = false;

	for ( int lambda = 0; lambda < sqrt( (double) N )/2; lambda++ ) {
		for ( int i = 0; i < N; i++ ) {
			taken[i] = false;
			ind[i] = i;
		}
		// clear tour matrix
		for ( int i = 0; i < N; i++ )
			for ( int j = 0; j < N; j++ ) 
				tour[N*i+j] = tour[N*j+i] = 0;

		taken[0] = true;

		int d;
		do {
			d = rand()%(N-1) + 1; // random starting node
		} while ( picked[d] );

		taken[d] = true;
		tour[N*0+d] = tour[N*d+0] = 1;

		picked[d] = true; // avoid picking twice the same node

		int count = 2, tlen = C[N*0+d];
		while ( count < N ) {
			// pick nearest neighbour
			int cost = INT_MAX, n = 0;

			for ( int i = 1; i < N; i++ ) { // 0 surely taken
				if ( C[N*d+i] < cost && !taken[i] ) {
					// taken[d] == true, so i != d
					cost = C[N*d+i];
					n = i;
				}
			}
			
			if ( n == 0 ) {
				cerr << "Error: TabuRoute" << endl;
				exit( EXIT_FAILURE );
			}

			taken[n] = true;
			tour[N*n+d] = tour[N*d+n] = 1;
			tlen += cost;
	
			d = n;
			count++;
		}

		// close the tour
		tour[N*d+0] = tour[N*0+d] = 1;
		tlen += C[N*d+0];

		// now 2-optimize the tour
		int v1 = tlen;
		int perm[4];
		for ( int i = 0; i < N; i++ ) {
			perm[0] = i;
			for ( int j = 0; j < N; j++ ) {
				if ( i == j )
					continue;

				perm[1] = j;

				for ( int k = 0; k < N; k++ ) {
					if ( j == k || i == k )
						continue;

					perm[2] = k;

					for ( int l = 0; l < N; l++ ) {
						if ( k == l || j == l || i == l  )
							continue;

						perm[3] = l;
						two_opt( tour, perm, tlen );
					}
				}
			}
		}

		int v2 = tlen;
		if ( doLog ) {
			log << "Lambda" << lambda << " tlen delta: " << v1-v2 << endl;
		}

		// and build routes
		// reset marks
		for ( int i = 0; i < N; i++ ) {
			taken[i] = false;
		}
		taken[0] = true;

		vector< vector<int> > tmpSol(M);
		d = 0;
		for ( int i = 0; i < M-1; i++ ) {
			int cleft = Cw;
			int vleft = Cv;
			do {
				added = false;
				for ( int n = 0; n < N; n++ ) {
					if ( tour[N*d+n] && !taken[n] && W[n] <= cleft && V[n] <= vleft ) {
						taken[n] = true;
						added = true;
						cleft -= W[n];
						vleft -= V[n];
						tmpSol[i].push_back( n );
						d = n;
						break;
					}
				}
			} while ( added );
		}

		// d is the last node inserted
		// last route again disregards capacity and may be unfeasible
		do {
			added = false;
			for ( int n = 0; n < N; n++ ) {
				if ( tour[N*d+n] && !taken[n] ) {
					taken[n] = true;
					added = true;
					tmpSol[M-1].push_back( n );
					d = n;
					break;
				}
			}
		} while ( added );

		if ( isFeasible( tmpSol ) )
			postOpt( tmpSol );

		if ( doLog ) {
			log << "F1 = " << F1( tmpSol ) << " ";
			log << "F2(1.0,1.0) = "<< F2( tmpSol, 1.0, 1.0 ) << endl << endl;
		}

		if ( F2( tmpSol, 1.0, 1.0 ) < F2( bestSol, 1.0, 1.0 ) )
			bestSol = tmpSol;
		if ( isFeasible( tmpSol ) && F1( tmpSol ) < F1( bestFsblSol ) )
			bestFsblSol = tmpSol;
	} // for lambda iterations


	// delete data structures used in the initialization phase
	delete [] ind;
	delete [] taken;
	delete [] tour;
	delete [] picked;

	if ( doLog ) {
		log << "== Phase 0 ==" << endl << endl;
		log << "Phase 0 took " << timer.timeElapsed() << " seconds" << endl << endl;
	
		if ( !bestSol.empty() ) {
			log << "Current best feasible solution:" << endl;
			for ( RouteSet::size_type i = 0; i < bestSol.size(); i++ ) {
				log << "Route #" << i << ": ";

				for ( Route::size_type j = 0; j < bestSol[i].size(); j++ )
					log << bestSol[i][j] << " ";
				log << endl;
			}
		}
		else
			log << "No feasible solution found" << endl;
		
		log << endl << "Best objective value: " << F1( bestFsblSol ) << endl;
		log << "Best objective2 value (alpha = beta = 1.0): " << F2( bestSol, 1.0, 1.0 ) << endl << endl;
	}

	if ( !silent ) {
		cout << "== Phase 0 ==" << endl << endl;

		cout << "Best objective value: " << F1( bestFsblSol ) << endl;
		cout << "Best objective2 value (alpha = beta = 1.0): " << F2( bestSol, 1.0, 1.0 ) << endl << endl;

		cout << "== Phase 1 ==" << endl << endl;
	}

	int p, q;
	
	// diversifying iterations
	//p = min( N-2, 10 ), q = min( N-1, 5*M );
	//p = N-2, q = min( N-1, 5*M );
	//p = min( N-2, 25 ), q = N-1;
	p = N-2, q = min( N-1, 5 );

	timer.restart();
	if ( !silent )
		cout << "Searching..." << endl;
	ts( iFact*(N-1), p, q, bestSol, bestFsblSol, lim );
	
	// get upper bound
	int bestVal = F1( bestFsblSol );

	if ( doLog ) {
		log << "== Phase 1 ==" << endl << endl;
		log << "Phase 1 took " << timer.timeElapsed() << " seconds" << endl << endl;

		// log phase 1 best solution
		for ( RouteSet::size_type i = 0; i < bestFsblSol.size(); i++ ) {
			log << "Route #" << i << ": ";

			for ( Route::size_type j = 0; j < bestFsblSol[i].size(); j++ ) {
				log << bestFsblSol[i][j] << " ";
			}
			log << endl;
		}
		log << "Objective value: " << bestVal << endl << endl;
	}
	
	if ( !silent ) {	
		cout << "Objective value: " << bestVal << endl << endl;
		cout << "== Phase 2 ==" << endl << endl;
	}

	// intensifying iterations
	//p = min( N-1, 10 ), q = N-1;
	p = N-2, q = N-1;

	timer.restart();
	if ( !silent )
		cout << "Searching..." << endl;
	ts( (N-1)*M, p, q, bestSol, bestFsblSol, lim );

	if ( doLog ) {
		log << "== Phase 2 ==" << endl << endl;
		log << "Phase 2 took " << timer.timeElapsed() << " seconds" << endl << endl;

		// log phase 2 best solution
		for ( RouteSet::size_type i = 0; i < bestFsblSol.size(); i++ ) {
			log << "Route #" << i << ": ";

			for ( Route::size_type j = 0; j < bestFsblSol[i].size(); j++ ) {
				log << bestFsblSol[i][j] << " ";
			}
			log << endl;
		}
	}

	bestVal = F1( bestFsblSol );

	if ( doLog ) {
		log << "Objective value: " << bestVal << endl << endl;
		log << "Trying to post-optimize the best feasbile solution found"
		<< endl << endl;
	}

	if ( !silent ) {
		cout << "Objective value: " << bestVal << endl << endl;
		cout << "Trying to post-optimize the best feasbile solution found"
			 << endl << endl;
	}

	postOpt( bestFsblSol );

	bestVal = F1( bestFsblSol );

	// Update problem variables
	if ( !bestFsblSol.empty() && bestVal < zub ) {
		zub = bestVal;
		x_bar = bestFsblSol;
	}

	for ( RouteSet::size_type i = 0; i < bestFsblSol.size(); i++ ) {
		if ( doLog )
			log << "Route #" << i << ": ";
		if ( !silent )
			cout << "Route #" << i << ": ";

		for ( Route::size_type j = 0; j < bestFsblSol[i].size(); j++ ) {
			if ( doLog )
				log << bestFsblSol[i][j] << " ";
			if ( !silent )
				cout << bestFsblSol[i][j] << " ";
		}
		if ( doLog ) 
			log << endl;
		if ( !silent )
			cout << endl;
	}

	if ( doLog ) 
		log << "Objective value: " << bestVal << endl << endl;
	if ( !silent )
		cout << "Objective value: " << bestVal << endl << endl;

	return 0;
}

/** TabuSearch - actual local search procedure */
void VrpSolver::ts( int tmax, int p, int q, RouteSet &bestSol, RouteSet &bestFsblSol, long lim )
{
	// tabu list tl[i][r] defines for how many iterations (i,r) is tabu
	int **tl = new int*[N];
	for ( int i = 0; i < N; i++ )
		tl[i] = new int[M];

	for ( int i = 0; i < N; i++ )
		for ( int j = 0; j < M; j++ )
			tl[i][j] = 0;

	double alpha = 1.0;
	double beta = 1.0;

	int lastWeightFsbl, lastWeightUnfsbl;
	int lastVolumeFsbl, lastVolumeUnfsbl;

	// taken marks to select q unique customers
	bool *taken = new bool[N-1];
	// for sorting customers
	int *ind = new int[N-1];

	lastWeightFsbl = lastWeightUnfsbl = 0;
	lastVolumeFsbl = lastVolumeUnfsbl = 0;

	// iteration solution
	RouteSet iSol = bestSol;

	int notImproved = 0, t = 0;
	//for ( int t = 0; t < tmax; t++ ) {
	while ( notImproved < tmax ) {

		//cout << "  iteration " << t << endl;
		t++;
		if ( t >= lim )	{
			//cout << "Iterations cap reached" << endl;
			break;
		}

		if ( t%10 == 0 ) {
			// weights penalties
			if ( lastWeightFsbl >= 10 )
				alpha = alpha / 2;
			else if ( lastWeightUnfsbl >= 10 )
				alpha = alpha * 2;

			// volume penalties
			if ( lastVolumeFsbl >= 10 )
				beta = beta / 2;
			else if ( lastVolumeUnfsbl >= 10 )
				beta = beta * 2;
		}

		// clear marks
		for ( int i = 0; i < N-1; i++ )
			taken[i] = false;

		vector<int> Z; /* Clients Set to be examined */
		while ( Z.size() < vector<int>::size_type( q ) ) {
			//int d = rand() % N;
			int d = rand() % (N-1);
			if ( !taken[d] ) {
				taken[d] = true;
				Z.push_back( d+1 );
			}
		}

		int val_f1 = F1( bestFsblSol );
		double val_f2 = F2( bestSol, alpha, beta );

		while ( !Z.empty() ) {
			// For each i in Z

			// Get i
			int i = Z.back();
			Z.pop_back();

			// Find route to which i belongs
			int ir = -1;
			for ( RouteSet::size_type j = 0; j < iSol.size() && ir == -1; j++ )
				for ( Route::size_type h = 0; h < iSol[j].size(); h++ ) 
					if ( iSol[j][h] == i )
						ir = j;

			if ( ir == -1 ) {
				cerr << "TabuSearch: node " << i
				     << " doesn't belong to any route" << endl;
				exit( EXIT_FAILURE );
			}

			// Sort other clients from closest to most distant
			
			for ( int j = 0; j < N-1; j++ )
				ind[j] = j+1; // do not consider 0

			sort( ind, ind+N-1, Nearest( i, N, C ) );
			// Obviously now ind[0] == i, since C[N*i+i] == 0

			// select candidate routes
			set<int> routes;
			for ( int j = 1; j <= p; j++ ) {
				// for each neighbour
				bool found = false;

				for ( RouteSet::size_type h = 0; h < iSol.size() && !found; h++ ) {
					// for each route in sol
					for ( Route::size_type k = 0; k < iSol[h].size() && !found; k++ ) {
						// for each node in the route
						if ( ind[j] == iSol[h][k] ) {
							found = true;
							routes.insert( h );
						}
							
					}
				}
			}
			// pick empty routes as candidates too
			for ( RouteSet::size_type j = 0; j < iSol.size(); j++ )
				if ( iSol[j].empty() )
					routes.insert( j );

			/* Now routes contains all candidate routes to host node i
			   we pick the best yielded solution (i, r*) */
			vector< vector<int> > bestNgh; // best neighbour
			bool updated = false;

			do {
				// put i in the best new route
				int r = *( --routes.end() );
				routes.erase( --routes.end() );

				// current neighbour solution
				vector< vector<int> > currNgh = putInRoute( i, r, iSol );

				// check if (i, r) was forbidden
				if ( tl[i][r] > 0 ) {
					// Move found in the tabu list
					if ( isFeasible( currNgh ) && F1( currNgh ) > F1( bestFsblSol ) )
						// is feasible but not better than what we have, skip
						continue;
					if ( !isFeasible( currNgh ) && F2( currNgh, alpha, beta ) > F2( bestSol, alpha, beta ) )
						// is unfeasible and not better than what we already have
						continue;
				}

				// we proceed with the local search, either (i,r) is not taboo
				// or it is but yields a "better" solution

				if ( F2( bestNgh, alpha, beta ) > F2( currNgh, alpha, beta ) ) {
					bestNgh = currNgh;
					updated = true;
				}

			} while ( !routes.empty() );

			// bestNgh contains the best neighbour of iSol

			if ( updated == false )
				// could not place i in any route, nothing has changed
				continue;
			else {
				// update solutions and counters
				if ( F2( bestNgh, alpha, beta ) < F2( bestSol, alpha, beta ) )
					bestSol = bestNgh;

				if ( isFeasible( bestNgh ) ) {
					lastWeightFsbl++;
					lastVolumeFsbl++;

					lastWeightUnfsbl = lastVolumeUnfsbl = 0;

					if ( F1( bestNgh ) < F1( bestFsblSol ) )
						bestFsblSol = bestNgh;
				}
				else {
					if ( feasibleLoad( bestNgh, weightCapacity ) ) {
						lastWeightFsbl++;
						lastWeightUnfsbl = 0;
					}
					else {
						lastWeightFsbl = 0;
						lastWeightUnfsbl++;
					}

					if ( feasibleLoad( bestNgh, volumeCapacity ) ) {
						lastVolumeFsbl++;
						lastVolumeUnfsbl = 0;
					}
					else {
						lastVolumeFsbl = 0;
						lastVolumeUnfsbl++;
					}
				}

				iSol = bestNgh;

				// declare (i, ir) tabu, ir is the route to which i belonged
				// before changing iSol
				int theta = 5 + ( rand()%6 );
				tl[i][ir] = theta;
			}

		} // while ( !Z.empty() )

		for ( int i = 0; i < N; i++ )
			for ( int j = 0; j < M; j++ )
				if ( tl[i][j] > 0 )
					tl[i][j]--;

		if ( val_f1 == F1( bestFsblSol ) && val_f2 == F2( bestSol, alpha, beta ) )
		//if ( val_f1 == F1( bestFsblSol ) )
			notImproved += 1; // nothing has changed
		else
			notImproved = 0;

	}

	delete [] taken;
	delete [] ind;

	for ( int i = 0; i < N; i++ )
		delete [] tl[i];
	delete [] tl;
}

/** Decide if a given set of routes is feasible */
bool VrpSolver::isFeasible( const RouteSet &S )
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

/** Decides if a given route is feasible with respect to the specified dimension */
bool VrpSolver::feasibleLoad( const RouteSet &S, VrpSolver::Dimension type )
{
	if ( S.empty() )
		return false;

	bool feasible = true;
	int *dem = ( type == weightCapacity ) ? W : V;
	int cap = ( type == weightCapacity ) ? Cw : Cv;

	for ( RouteSet::size_type i = 0; i < S.size() && feasible; i++ ) {
		// for each route
		int load = 0;

		for ( Route::size_type j = 0; j < S[i].size(); j++ ) {
			load += dem[ S[i][j] ];
		}

		feasible = ( load <= cap );
	}

	return feasible;
}

/** TabuSearch - auxiliary objective function */
double VrpSolver::F2( const RouteSet &solution, double alpha, double beta )
{
	if ( solution.empty() )
		return DBL_MAX;

	int sumWeight = 0;
	int sumVolume = 0;

	for ( RouteSet::size_type i = 0; i < solution.size(); i++ ) {
		int wLoad = 0;
		int vLoad = 0;

		for ( Route::size_type j = 0; j < solution[i].size(); j++ ) {
			// Sum demands of current route
			wLoad += W[solution[i][j]];
			vLoad += V[solution[i][j]];
		}
		
		sumWeight += max( 0, wLoad-Cw );
		sumVolume += max( 0, vLoad-Cv );
	}

	return F1( solution, true ) + alpha*sumWeight + beta*sumVolume;
}

/** TabuSearch - objective function */
int VrpSolver::F1( const vector< vector<int> > &solution, bool ignoreUnf )
{
	if ( solution.empty() || ( !ignoreUnf && !isFeasible( solution ) ) )
		return INT_MAX;

	int sum = 0;

	for ( RouteSet::size_type i = 0; i < solution.size(); i++ ) {
		if ( !solution[i].empty() ) {
			// sum 0 -> solution[i][0]
			sum += C[N*0 + solution[i][0]];

			for ( Route::size_type j = 1; j < solution[i].size(); j++ ) 
				//sum s[i][j-1] -> s[i][j]
				sum += C[N*solution[i][j-1] + solution[i][j]];

			sum += C[N*solution[i].back() + 0];
		}
	}

	return sum;
}

/** Insert a node in a route */
RouteSet VrpSolver::putInRoute( int nd, int ind, const RouteSet &sol )
{
	RouteSet s = sol; // Shadow solution

	/* erase node nd from the current solution */
	Route::iterator target;
	RouteSet::size_type inR;
	bool found = false;
	for ( inR = 0; inR < s.size(); inR++ ) {
		for ( target = s[inR].begin(); target != s[inR].end(); target++ )
			if ( *target == nd ) {
				found = true;
				break;
			}
		if ( found )
			break;
	}
	if ( found )
		s[inR].erase( target );

	// if s is empty, simply put nd in s[ind]
	if ( s[ind].empty() ) {
		s[ind].push_back( nd );
		return s;
	}

	/* Find best position to insert node nd 

	   extraMileage: suppose to insert node k in position i, then the extra
	   cost is found by calculating C(i-1,k) + C(k,i) - C(i-1,i) */
	int extraMileage;
	int bestEm; // best insertion cost
	int bestPos; // best insertion index

	// if first in the route
	bestPos = 0;
	bestEm = extraMileage = C[N*0 + nd] + C[N*nd + s[ind][0]] - C[N*0 + s[ind][0]];
	//bestPosPrice = basePrice + extraMileage;

	// in between other nodes
	for ( Route::size_type i = 1; i < s[ind].size(); i++ ) {
		extraMileage = C[N*s[ind][i-1] + nd] + C[N*nd + s[ind][i]] -
			C[N*s[ind][i-1] + s[ind][i]];
		if ( extraMileage < bestEm ) {
			bestPos = i;
			bestEm = extraMileage;
		}
	}

	// if last in the route
	extraMileage = C[N*s[ind].back() + nd] + C[N*nd + 0] - 
		C[N*s[ind].back() + 0];
	// Update route
	if ( extraMileage < bestEm )
		s[ind].push_back( nd );
	else 
		s[ind].insert(  s[ind].begin() + bestPos, nd );

	return s;
}

void VrpSolver::two_opt( double *tour, int *perm, int &tlen )
{
	int i, j, k, l;
	i = perm[0];
	j = perm[1];
	k = perm[2];
	l = perm[3];
	if ( !( i != j && i != k && i != l && j != k && j != l && k != l ) ) {
		return;
	}
	if ( !tour[N*i+j] || !tour[N*k+l] ) {
		// at least one of the 2 arcs is not in the graph
		return;
	}

	tour[N*i+j] = tour[N*j+i] = tour[N*k+l] = tour[N*l+k] = 0;

	VertexSet *component = VrpData::dfs( tour, i, N );
	// attach i to k and j to l, if i and k are in the same component
	// swap k and l
	for ( VertexSet::size_type h = 0; h < component->size(); h++ ) {
		if ( (*component)[h] == k ) {
			swap( k, l );
			break;
		}
	}

	delete component;

	int tnew = tlen - C[N*i+j] - C[N*k+l] + C[N*i+k] + C[N*j+l];

	if ( tnew < tlen ) {
		// apply the move and update tlen
		tour[N*i+k] = tour[N*k+i] = tour[N*j+l] = tour[N*l+j] = 1;
		tlen = tnew;
	}
	else
		// restore tour
		tour[N*i+j] = tour[N*j+i] = tour[N*k+l] = tour[N*l+k] = 1;
}

void VrpSolver::postOpt( RouteSet &solution )
{
	/* 2-Optimize each tour in the solution */
	double *tour = new double[N*N];

	for ( RouteSet::size_type h = 0; h < solution.size(); h++ ) {

		if ( solution[h].empty() )
			continue;

		// build tour for the current solution
		for ( int i = 0; i < N*N; i++ )
			tour[i] = 0.0;

		int d = 0; // last node inserted
		int n = solution[h].front(); // next node to insert
		int tlen = 0; // tour length

		// depot->solution[h].front
		tour[N*d+n] = tour[N*n+d] = 1.0;
		tlen += C[N*d+n];

		// solution[h][t-1]->solution[h][t]
		d = n;
		for ( Route::size_type t = 1; t < solution[h].size(); t++ ) {
			// n is the next node to insert in the tour
			n = solution[h][t];

			tour[N*d+n] = tour [N*n+d] = 1.0;
			tlen += C[N*d+n];

			d = n;
		}

		//solution[h].back->depot
		d = solution[h].back();
		n = 0;
		tour[N*d+n] = tour[N*n+d] = 1.0;
		tlen += C[N*d+n];

		// 2 optimize the tour
		int perm[4];
		for ( int i = 0; i < N; i++ ) {
			perm[0] = i;
			for ( int j = 0; j < N; j++ ) {
				if ( i == j )
					continue;

				perm[1] = j;

				for ( int k = 0; k < N; k++ ) {
					if ( j == k || i == k )
						continue;

					perm[2] = k;

					for ( int l = 0; l < N; l++ ) {
						if ( k == l || j == l || i == l  )
							continue;

						perm[3] = l;
						two_opt( tour, perm, tlen );
					}
				}
			}
		}

		// update solution[h]
		vector<int> *path = VrpData::dfs( tour, 0, N );
		solution[h] = *path;
		solution[h].erase( solution[h].begin() ); // remove 0

		delete path;
	}

	delete [] tour;
}
