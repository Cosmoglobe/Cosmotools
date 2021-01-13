#include <cstdlib>
#include <cstdio>
using namespace std;
int main(int argn, char ** args)
{
	double mjd = atof(args[1]);
	printf("%ld\n", (long int)((mjd-40587)*24*60*60));
}
