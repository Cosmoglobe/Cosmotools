// Given a l2-database, the chicago ces database and a ces number, output
// the corresponding oslo ces number, or -1 -1 if none apply.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

struct CES { CES(int r = -1, int s = -1):run(r),seg(s) {} int run, seg; };
struct Range { Range(double r1 = 0, double r2 = 0):from(r1),to(r2) {} double from,to; };
bool operator<(const CES & a, const CES & b) { return a.run == b.run ? a.seg < b.seg : a.run < b.run; }
bool operator<(const Range & a, const Range & b) { return a.from < b.from; }

// We assume that no ranges in the table overlap.
template <class T> map<Range,T> overlapping_elements(const map<Range,T> & table, Range r)
{
	map<Range,T> res;
	typedef typename map<Range,T>::const_iterator it;
	it fstart = table.lower_bound(r);
	it bstart = fstart; --bstart;
	if(fstart != table.begin() && bstart->first.to > r.from)
		res[bstart->first] = bstart->second;
	for(it i = fstart; i != table.end() && i->first.from < r.to; i++)
		res[i->first] = i->second;
	return res;
}

void translate(const CES & ces, map<CES,Range> & from,
	map<Range,CES> & to, bool verbose = true)
{
	// Is the input ces defined at all?
	map<CES,Range>::const_iterator ci = from.find(ces);
	if(ci == from.end())
	{
		if(verbose) printf("undefined\n");
		else printf("-1 0\n");
		return;
	}
	const Range & r = from[ces];
	map<Range, CES> res = overlapping_elements(to, r);
	typedef map<Range, CES>::iterator it;
	if(!verbose)
	{
		if(res.empty()) printf("-1 1\n");
		else
		{
			for(it i = res.begin(); i != res.end(); i++)
				printf("%d %d  ", i->second.run, i->second.seg);
			printf("\n");
		}
	}
	else
	{
		printf("%d %d %13.7lf %13.7lf => %d\n", ces.run, ces.seg,
			r.from, r.to, res.size());
		int j = 1;
		for(it i = res.begin(); i != res.end(); i++, j++)
			printf("  %d %13.7lf %13.7lf: %d %d\n", j, i->first.from,
				i->first.to, i->second.run, i->second.seg);
	}
}

void help()
{
	printf("Syntax: ces_translate [options] [run seg]\n"
		"Translates from chicago run and seg to oslo\n"
		"run and seg. If run and seg are not specified\n"
		"run and seg pairs are read from standard input.\n"
		"Note: In general, there is no 1-1 mapping between\n"
		"chicago and oslo ceses, so more than one output\n"
		"is possible for each input. Also, if the input\n"
		"is not defined, this is indicated by a run number\n"
		"of -1.\n"
		"Options:\n"
		" -r: Reverse the lookup direction.\n"
		" -c: Specify a chicago ces database manually.\n"
		" -o: Specify an oslo level2 database manually.\n"
		" -v: More verbose output. Prints the mjd ranges matched against.\n"
		" -h: Display this message and exit.\n");
}

int main(int argc, char ** argv)
{
	vector<char*> args;
	char * l2filename = "/projects/quiet/runlist.txt",
		* chicagofilename = "/data4/quiet/calib_data/chicago/ces_defs.txt";
	bool reverse = false;
	bool verbose = false;
	for(char ** i = argv+1; *i; i++)
		if(!strcmp("-o",*i)) l2filename = *++i;
		else if(!strcmp("-c",*i)) chicagofilename = *++i;
		else if(!strcmp("-h",*i)) { help(); return 1; }
		else if(!strcmp("-r",*i)) reverse ^= true;
		else if(!strcmp("-v",*i)) verbose ^= true;
		else args.push_back(*i);

	// Read chicago db
	map<CES,Range> c2r_chicago;
	map<Range,CES> r2c_chicago;
	int c, s; double r1, r2;
	FILE * cfile = fopen(chicagofilename,"r");
	while(fscanf(cfile, "%d %d %lf %lf", &c, &s, &r1, &r2) == 4)
	{
		c2r_chicago[CES(c,s)] = Range(r1,r2);
		r2c_chicago[Range(r1,r2)] = CES(c,s);
	}
	fclose(cfile);

	// Read oslo db
	map<CES,Range> c2r_oslo;
	map<Range,CES> r2c_oslo;
	ifstream oslo(l2filename);
	string line;
	while(getline(oslo,line))
	{
		stringstream ss; ss << line;
		if(line.substr(0,2) == "  " && line[4] != ' ') ss >> c >> r1 >> r2;
		c2r_oslo[CES(c,0)] = Range(r1,r2);
		r2c_oslo[Range(r1,r2)] = CES(c,0);
	}
	oslo.close();

	map<CES,Range> & from_c2r = !reverse ? c2r_chicago : c2r_oslo;
	map<Range,CES> & from_r2c = !reverse ? r2c_chicago : r2c_oslo;
	map<CES,Range> & to_c2r   = !reverse ? c2r_oslo    : c2r_chicago;
	map<Range,CES> & to_r2c   = !reverse ? r2c_oslo    : r2c_chicago;

	// Perform the lookup
	if(args.empty())
	{
		while(getline(cin, line))
		{
			CES ces(0,0);
			int n = sscanf(line.c_str(), "%d %d", &ces.run, &ces.seg);
			if(n >= 1) translate(ces, from_c2r, to_r2c, verbose);
		}
	}
	else
	{
		CES ces(atoi(args[0]), args.size() > 1 ? atoi(args[1]) : 0);
		translate(ces, from_c2r, to_r2c, verbose);
	}
}
