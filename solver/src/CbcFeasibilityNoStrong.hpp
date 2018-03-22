#ifndef NoStrong_H
#define NoStrong_H

#include "CbcFeasibilityBase.hpp"
#include "CbcModel.hpp"

#include "VrpData.hpp"

#define CHECK_SOLUTION 0

/** Class to prevent integer solutions found during strong branching to be
    considered feasible, since an integer solution found during strong branching
    may very well be cut by the cut generator */
class CbcFeasibilityNoStrong : public CbcFeasibilityBase {
public:
	VrpData graph;
	// Default Constructor
	CbcFeasibilityNoStrong( const VrpData &g ) : graph( g ) { }

	virtual ~CbcFeasibilityNoStrong() { }
	// Copy constructor
	CbcFeasibilityNoStrong( const CbcFeasibilityNoStrong &rhs ) : graph( rhs.graph ) { }

	// Assignment operator
	CbcFeasibilityNoStrong & operator=( const CbcFeasibilityNoStrong& rhs ) {
		return *this;
	}

	/// Clone
	virtual CbcFeasibilityBase *clone() const {
		return new CbcFeasibilityNoStrong( graph );
	}

	/** On input mode:
	 0 - called after a solve but before any cuts
	-1 - called after strong branching
	Returns :
	 0 - no opinion
	-1 - pretend infeasible
	 1 - pretend integer solution */
  virtual int feasible( CbcModel * model, int mode );
};

int CbcFeasibilityNoStrong::feasible(CbcModel * model, int mode)
{
#if CHECK_SOLUTION
	if ( mode == -1 ) {

	const double *x = model->getColSolution();
	//const double *x = solution;
	int numberColumns = model->getNumCols();

	// check if the solution is connected
	VrpData solGraph( graph );
	solGraph.buildResidualGraph( x );

	vector<int> *connected = solGraph.dfs( 0 );
	int size = connected->size();
	delete connected;

	if ( size < solGraph.N ) {
		// solution is not connected, pretend infeasible
		return mode;
	}

	// check if the solution is integer
	// this should not be necessary
	bool integer = true;
	for ( int i = 0; i < numberColumns && integer; i++ ) {
		if ( ( x[i] > 0.0L && x[i] < 1.0L ) || ( x[i] > 1.0L && x[i] < 2.0L ) )
			integer = false;
	}

	if ( integer ) {
		// check if the routes are feasible
		vector< vector<int> > routes = solGraph.getRoutes();

		if ( !solGraph.isFeasible( routes ) )
			return mode;
		else {
			// solution is feasible
			return 1;
		}

	}
	else
		return mode;

	}
	else
#endif
		return mode;
}

#endif
