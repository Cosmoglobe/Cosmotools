#include <handy.h>
using namespace std;

// Read an ascii powspec and output it as fits
int main(int argn, char ** argv)
{
	try {
		vector<char*> args;
		char * psfile = 0, * ofile = 0;
		int ordering = DIAGONAL_FIRST, skipl = 0;
		bool scaled = true;
		for(char ** i = argv+1; *i; i++)
			if(!strcmp(*i,"-s")) scaled ^= true;
			else if(!strcmp(*i,"-l")) skipl ^= true;
			else args.push_back(*i);
		if(args.size() < 1) serror("Usage: ps2fits [options] ps ofile");
		psfile = args[0];
		ofile  = args[1];

		// Generate our starting map, using the TT column of a file
		vector<PowSpec> p = fits2powspec(psfile, 0, -1);
		powspec2ascii(p, ofile, skipl, scaled);
	} catch(const Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); }
}
