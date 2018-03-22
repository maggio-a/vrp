/************************************************************
  Autore: Andrea Maggiordomo - mggndr89@gmail.com
  Università di Pisa - Informatica Applicata - Logistica
  2012
*************************************************************/

#include "utils.hpp"

char *os2ch( std::ostringstream& os )
{
	std::string from = os.str();

	int len = from.length();
	char *to = new char[len+1];

	for ( int i = 0; i < len; i++ )
		to[i] = from[i];

	to[len] = '\0';

	return to;
}
