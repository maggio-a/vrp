/************************************************************
  Autore: Andrea Maggiordomo - mggndr89@gmail.com
  Università di Pisa - Informatica Applicata - Logistica
  2012
*************************************************************/

#ifndef VrpHeuristicSeparator_H
#define VrpHeuristicSeparator_H

#include "CglCutGenerator.hpp"

#include "VrpData.hpp"
#include "VrpTypes.hpp"

#include <vector>

using std::vector;

class VrpSep : public CglCutGenerator {

public:
	mutable VrpData graph; // problem info
	mutable double gomoryTime;     // time spent by CglGomory
	mutable long gomoryTotal;      // total calls to Cglgomory
	mutable long gomoryEffective;  // effective calls to CglGomory

	/** Generate cuts for the solver si
	The cuts are inserted in the collection cs */
	virtual void generateCuts( const OsiSolverInterface &si, OsiCuts &cs,
		const CglTreeInfo info=CglTreeInfo() ) const;

	/** Constructor */
	VrpSep( const VrpData &, bool gomory=false );

	/* Copy constructor */
	VrpSep( const VrpSep & );

	/** Clone */
	VrpSep *clone() const;

	/** Assignment operator */ 
	VrpSep & operator=( const VrpSep& rhs );

	/** Destructor */
	~VrpSep() { }

	/// Return true if needs optimal basis to do cuts
	bool needsOptimalBasis() const;

private:

	void cleanCutSet( vector<VertexSet*> *cuts ) const;

	void detachedComponents( vector<VertexSet*> *cuts ) const;

	/* Separation heuristics 

	   ccHeuristic - Connected components heuristic described in
	   T.K. Ralphs, L. Kopman, W.R. Pulleyblank, and L.E. Trotter Jr.
	   ''On the Capacitated Vehicle Routing Problem'' */
	void ccHeuristic( vector<VertexSet*> *cuts, double eps, bool findAll=false ) const;

	/* Shrinking heuristic described in 	
	   T.K. Ralphs, L. Kopman, W.R. Pulleyblank, and L.E. Trotter Jr.
	   ''On the Capacitated Vehicle Routing Problem'' */
	void shrinkingHeuristic( vector<VertexSet*> *cuts, double eps, bool findAll=false ) const;
	
	/* Greedy rounded capacity heuristic 
	   Something similar is described in P. Augerat et al.
	   ''Separating capacity constraints in the CVRP using tabu search'' */
	void greedyRoundCap( vector<VertexSet*> *cuts, double eps, bool findAll=false ) const;

	/* Converts sets of customers in cuts and adds them in cs */
	void addCuts( vector<VertexSet*> *cuts, OsiCuts &cs ) const;
	
	bool useGomory;
	int numCols;

};

#endif
