/************************************************************
  Autore: Andrea Maggiordomo - mggndr89@gmail.com
  Università di Pisa - Informatica Applicata - Logistica
  2012
*************************************************************/

#include "utils.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

using namespace std;

inline double rand_val() {
	return ( double( rand() ) / RAND_MAX );
}

/* Parameter file structure:
   <VolumeCapacity>
   <VolumeDemandRatio> // demands-capacity tightness
   <RangeCapacity> // vf_i = W[i]/wCap
					  V[i] ~ vf_i*vCap +|- range*( vf_i*vCap )
					  with range in (0, RangeCapacity)

   */

int main( int argc, char **argv )
{
	int vehiclesNumber;
	int weightCap, volumeCap;
	double volumeDemandRatio;
	double rangeCap;
	int N;

	int random = 0;
	int seed = 12345;

	switch( argc ) {
	case 5:
		get_arg( argv[4], seed );
	case 4:
		get_arg( argv[3], random );
	case 3:
		break;
	default:
		cerr << "Usage: " << argv[0] << " param inst [random] [seed]" << endl
			<< "param  :: input parameters file" << endl
			<< "inst   :: vrp input instance" << endl
			<< "random :: 0 volume demands according to parameters (default), !0 random proportions" << endl
			<< "seed   :: default " << seed << endl; 
		exit( EXIT_FAILURE );
	}

	srand( seed );

	// open parameters file
	ifstream param( argv[1] );
	if ( !param.is_open() ) {
		cerr << "Cannot open file " << argv[1] << endl;
		exit( EXIT_FAILURE );
	}

	// open instance
	ifstream instance( argv[2] );

	if ( !instance.is_open() ) {
		cerr << "Cannot open file " << argv[2] << "\n";
		exit( EXIT_FAILURE );
	}

	// read parameters
	instance >> vehiclesNumber;
	if ( param.fail() ) {
		cerr << "Error reading vehiclesNumber" << endl;
		exit( EXIT_FAILURE );
	}
	if ( vehiclesNumber <= 0 ) {
		cerr << "Error: vehicles number must be positive" << endl;
		exit( EXIT_FAILURE );
	}

	instance >> weightCap;
	if ( param.fail() ) {
		cerr << "Error reading weightCap" << endl;
		exit( EXIT_FAILURE );
	}
	if ( weightCap <= 0 ) {
		cerr << "Error: weight capacity must be positive" << endl;
		exit( EXIT_FAILURE );
	}

	// from param
	param >> volumeCap;
	if ( param.fail() ) {
		cerr << "Error reading volumeCap" << endl;
		exit( EXIT_FAILURE );
	}
	if ( volumeCap <= 0 ) {
		cerr << "Error: volume capacity must be positive" << endl;
		exit( EXIT_FAILURE );
	}

	param >> volumeDemandRatio;
	if ( param.fail() ) {
		cerr << "Error reading volumeDemandRatio" << endl;
		exit( EXIT_FAILURE );
	}
	if ( volumeDemandRatio <= 0 || volumeDemandRatio >= 1 ) {
		cerr << "Error: demand/volume ratio capacity must be in ( 0 , 1 )" << endl;
		exit( EXIT_FAILURE );
	}

	param >> rangeCap;
	if ( param.fail() ) {
		cerr << "Error reading rangeCap" << endl;
		exit( EXIT_FAILURE );
	}
	if ( rangeCap < 0 || rangeCap >= 0.99 ) {
		cerr << "Error: rangeCap must be in [ 0 , 0.99 )" << endl;
		exit( EXIT_FAILURE );
	}

	instance >> N;
	if ( param.fail() ) {
		cerr << "Error reading n" << endl;
		exit( EXIT_FAILURE );
	}
	if ( N <= 0 ) {
		cerr << "Error: number of clients must be positive" << endl;
		exit( EXIT_FAILURE );
	}

	// get coordinates
	int *xpos = new int[N];
	int *ypos = new int[N];
	int id;

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

	// read weight demands
	vector<int> W( N );
	for ( int i = 0; i < N; i++ ) {
		instance >> id; /* Reads the node id */
		if ( instance.fail() ) {
			cerr << "Error reading node info from instance file\n";
			exit( EXIT_FAILURE );
		}

		instance >> W[i];
		if ( instance.fail() ) {
			cerr << "Error reading node info from instance file\n";
			exit( EXIT_FAILURE );
		}
	}

	// now assign volume demands
	vector<int> V( N );
	V[0] = 0;

	// calculate total demand
	int totVolumeDem = int( vehiclesNumber*volumeCap*volumeDemandRatio );
	int vMean = totVolumeDem / (N-1);

	// get fractions for each customer
	vector<double> vf(N);
	vf[0] = 0;
	for ( int i = 1; i < N; i++ ) {
		vf[i] = double( W[i] )/weightCap;
	}

	for ( int i = 1; i < N; i++ ) {

		int sign = ( rand_val() >= 0.5 ) ? +1 : -1;

		if ( random ) {
			do {
				V[i] = vMean + sign*int( rand_val() * vMean );
			} while ( V[i] <= 0 );
			totVolumeDem -= V[i];
		}
		else {
			// get base assignment
			int v_i = max( int( vf[i]*volumeCap ), 1 );

			double range;
			do {
				range = rand_val();
			} while ( range > rangeCap );

			// assign volume demand to i
			do {
				V[i] = v_i + sign*int( range*v_i );
			} while ( V[i] <= 0 );
			
			totVolumeDem -= V[i];
		}
	}

	// Fix volume demands
	int unchanged = 0; // if too many iterations without adjusting, adjust anyway
	while ( totVolumeDem != 0 ) {
		//int i = int( randVar() * (N-1) ) + 1; // sometimes generates N
		int i;
		while( ( i = int( rand_val() * (N-1) ) + 1 ) >= N );

		if ( random ) {
			if ( totVolumeDem > 0 ) {
				V[i]++;
				totVolumeDem--;
			}
		if ( totVolumeDem < 0 && V[i] > ( double( 1 ) / (N-1) )*volumeCap ) {
				// avoid decreasing too much too many
				V[i]--;
				totVolumeDem++;
			}
		}
		else {
			int v_i = max( int( vf[i]*volumeCap ), 1 );

			// determine bounds on V[i]
			int ViUpper = v_i + int( rangeCap*v_i );
			int ViLower = v_i - int( rangeCap*v_i );

			if ( totVolumeDem > 0 && (V[i] < ViUpper || unchanged > 2*N) ) {
				// add only if V[i] doesn't exceed range
				V[i]++;
				totVolumeDem--;
				
				unchanged = 0;
			}
			else if ( totVolumeDem < 0 && V[i] > ViLower ) {
				// decrease only if V[i] doesn't exceed range
				V[i]--;
				totVolumeDem++;
				
				unchanged = 0;
			}
			else {
				unchanged++;
			}
		}
	}

	assert( totVolumeDem == 0 );


	// output file
	//ostringstream os;
	//os << "n" << N << "-k" << vehiclesNumber << ".data";
	//char *outName = os2ch( os );

	ofstream out( /*outName*/"out.txt", ios_base::out|ios_base::trunc );

	if ( !out.is_open() ) {
		cerr << "Cannot open output file " << /*outName*/"out.txt" << endl;
		exit( EXIT_FAILURE );
	}

	//delete [] outName;

	// depot info section
	out << vehiclesNumber << endl << endl;
	out << weightCap << " " << volumeCap << endl << endl;

	// graph section
	// dimension
	out << N << endl << endl;
	// output coordinates
	for ( int i = 0; i < N; i++ )
		out << " " << i+1 << " " << xpos[i] << " " << ypos[i] << endl;
	
	// output demands
	out << endl;	
	for ( int i = 0; i < N; i++ ) {
		out << i+1 << " " << W[i] << " " << V[i] << endl;
	}

	//print memo of the parameters used:
	out << endl << "#PARAMS: VOLUMECAP=" << volumeCap
		<< " VOLUMEDEMANDRATIO=" << volumeDemandRatio
		<< " RANGECAP=" << rangeCap
		<< " RANDOM=" << (random != 0)
		<< " SEED=" << seed << endl;

	return 0;
}
