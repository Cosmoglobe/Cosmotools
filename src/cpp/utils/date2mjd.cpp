#include <cstdio>
#include <cstdlib>
#include <ctime>

using namespace std;

void memzero(char * a, int n) { for(int i = 0; i < n; i++) a[i] = 0; }
time_t toffset()
{
	static time_t foo = -1;
	if(foo != -1) return foo;
	tm tim; memzero((char*)&tim, sizeof(tim));
	strptime("1970-1-1 00:00:00", "%Y-%m-%d %H:%M:%S", &tim);
	foo = mktime(&tim);
	return foo;
}
double date2mjd(char * date, char * fmt = "%Y-%m-%d %H:%M:%S")
{
	tm tim; memzero((char*)&tim, sizeof(tim));
	if(!strptime(date, fmt, &tim))
		fprintf(stderr, "Could not parse %s according to %s", date, fmt);
	time_t t = mktime(&tim) - toffset();
	return t/86400.0 + 40587.0;
}

int main(int argn, char ** args)
{
	if(argn > 2) printf("%lf\n", date2mjd(args[1], args[2]));
	else printf("%lf\n", date2mjd(args[1]));
};
