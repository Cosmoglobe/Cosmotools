/* Given a CES list, a list of mjds where the phase switch status changes,
 * and a delay in seconds wanted after each change. Adjusts the starting
 * times of each ces as necessary, printing a warning if more than 1% of
 * the CES is lost.
 */

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cmath>
using namespace std;

static const double inf = 1.0/0.0;
int min(int a, int b) { return a < b ? a : b; }

struct CES { double mjd[2]; string rest; };
istream & operator>>(istream & in, CES & ces)
{
	string s;
	if(getline(in,s))
	{
		stringstream ss(s);
		ss >> ces.mjd[0] >> ces.mjd[1];
		getline(ss, ces.rest);
	}
	return in;
}
bool operator<(const CES & a, const CES & b) { return a.mjd[0] < b.mjd[0]; }

// Returns -1 if a is completely before b, 1 if a is completely after,
// and 0 if they overlap.
int compare(const CES & a, const CES & b)
{ return a.mjd[1] < b.mjd[0] ? -1 : a.mjd[0] > b.mjd[1] ? 1 : 0; }

void output(int res, const CES & a, const CES & b)
{
	if(res < 0) printf("%13.7lf %13.7lf %s\n", a.mjd[0], a.mjd[1], a.rest.c_str());
	else if(res > 0) return;
	else printf("%13.7lf %13.7lf %s\n", b.mjd[1], a.mjd[1], a.rest.c_str());
}

int main(int argc, char ** argv)
{
	// Prepare the data
	char * fn1 = argv[1], * fn2 = argv[2];
	double delay = atof(argv[3])/24/60/60;
	ifstream f1(fn1), f2(fn2);
	CES c; vector<CES> v1, v2;
	double d;
	while(f1 >> c) v1.push_back(c);
	while(f2 >> d) {
		CES ces;
		ces.mjd[0] = d;
		ces.mjd[1] = d + delay;
		v2.push_back(ces);
	}
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());

	// Add fake stop-ces at the end
	CES stop = { { inf, inf }, "stop" };
	v1.push_back(stop);
	v2.push_back(stop);

	// Sweep through the files in parallel, keeping
	// time synced.
	int i1 = 0, i2 = 0;
	for(;;)
	{
		if(i1 == v1.size()-1 && i2 == v2.size()-1) break;
		int res = compare(v1[i1],v2[i2]);
		output(res, v1[i1], v2[i2]);
		if(res < 0) i1++;
		else if(res > 0) i2++;
		else
		{
			int r1 = compare(v1[i1+1],v2[i2]),
				r2 = compare(v1[i1],v2[i2+1]);
			if(r1 != 0) i2++;
			if(r2 != 0) i1++;
		}
	}
}
