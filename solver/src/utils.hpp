/************************************************************
  Autore: Andrea Maggiordomo - mggndr89@gmail.com
  Università di Pisa - Informatica Applicata - Logistica
  2012
*************************************************************/

#ifndef Misc_Utils_H
#define Misc_Utils_H

#include <iostream>
#include <sstream>
#include <cstdlib>

template<class T> void get_arg( const char *arg, T &var )
{
	std::istringstream is( arg );
	is >> var;
	if ( is.fail() ) {
		std::cerr << "Error reading input parameters" << std::endl;
		exit( EXIT_FAILURE );
	}
}

char *os2ch( std::ostringstream& os );

#endif
