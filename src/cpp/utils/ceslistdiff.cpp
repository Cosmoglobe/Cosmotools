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
static const double s2mjd = 1.0/(24*3600);
int min(int a, int b) { return a < b ? a : b; }

struct CES { double mjd[2]; string name; };
istream & operator>>(istream & in, CES & ces)
{
	string s;
	if(getline(in, s))
	{
		stringstream ss(s);
		ss >> ces.mjd[0] >> ces.mjd[1];
		while(ss >> ces.name);
	}
	return in;
}
bool operator<(const CES & a, const CES & b) { return a.mjd[0] < b.mjd[0]; }

// Returns -1 if a is completely before b, 1 if a is completely after,
// and 0 if they overlap.
int compare(const CES & a, const CES & b)
{ return a.mjd[1] < b.mjd[0] ? -1 : a.mjd[0] > b.mjd[1] ? 1 : 0; }

void print_post(int res, const CES & a, const CES & b)
{
	if(res != 0) printf("___");
	else
	{
		string sa, sb;
		for(int i = 0; i < a.name.size(); i++)
			if(a.name[i] != '_') sa += a.name[i];
		for(int i = 0; i < b.name.size(); i++)
			if(b.name[i] != '_') sb += b.name[i];
		char x = sa == sb ? '_' : 'o';
		double df = fabs(a.mjd[0]-b.mjd[0]);
		double dt = fabs(a.mjd[1]-b.mjd[1]);
		char y = df > 50*s2mjd ? 'F' : df > 10*s2mjd ? 'f' : '_';
		char z = dt > 50*s2mjd ? 'T' : dt > 10*s2mjd ? 't' : '_';
		printf("%c%c%c", x, y, z);
	}
}
enum { NORMAL, EQUAL, EMPTY };
void print_ces(const CES & a, int type)
{
	switch(type)
	{
	case(NORMAL): printf("%13.7lf %13.7lf %8s ", a.mjd[0], a.mjd[1], a.name.c_str()); break;
	case(EQUAL):  printf("%36s ", "===================================="); break;
	case(EMPTY):  printf("%36s ", ""); break;
	}
}

void output(int res, const CES & a, const CES & b, bool aeq, bool beq)
{
	print_ces(a, res > 0 ? EMPTY : aeq ? EQUAL : NORMAL);
	print_ces(b, res < 0 ? EMPTY : beq ? EQUAL : NORMAL);
	print_post(res, a, b);
	printf("\n");
}

int main(int argc, char ** argv)
{
	// Prepare the data
	char * fn1 = argv[1], * fn2 = argv[2];
	ifstream f1(fn1), f2(fn2);
	CES c; vector<CES> v1, v2;
	while(f1 >> c) v1.push_back(c);
	while(f2 >> c) v2.push_back(c);
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());

	// Add fake stop-ces at the end
	CES stop = { { inf, inf }, "stop" };
	v1.push_back(stop);
	v2.push_back(stop);

	// Sweep through the files in parallel, keeping
	// time synced.
	int i1 = 0, i2 = 0, i1p = -1, i2p = -1;
	for(;;)
	{
		if(i1 == v1.size()-1 && i2 == v2.size()-1) break;
		int res = compare(v1[i1],v2[i2]);
		output(res, v1[i1], v2[i2], i1==i1p, i2==i2p);
		if(res <= 0) i1p = i1;
		if(res >= 0) i2p = i2;
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
