#include <cstdlib>
#include <cstdio>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
using namespace std;

static const double inf = 1.0/0.0;

vector<string> tokenize(const string & line, const string & sep)
{
	vector<string> result;
	int i = line.find_first_not_of(sep);
	int j = line.find_first_of(sep, i);
	while(i != string::npos)
	{
		result.push_back(line.substr(i,j-i));
		i = line.find_first_not_of(sep, j);
		j = line.find_first_of(sep,i);
	}
	return result;
}

int indent(const string & s, const string & white = " ")
{ return s.find_first_not_of(white); }

typedef pair<int,int> diode;

set<diode> all_diodes(int nmod, int ndi)
{
	set<diode> res;
	for(int i = 0; i < nmod; i++)
		for(int j = 0; j < ndi; j++)
			res.insert(diode(i,j));
	return res;
}

static const int any = -1;
static const double dany = inf;

struct Cut { double from, to; set<int> run, seg, mod, di; };
double   get_double(const string & s) { return s == "all" ? dany : atof(s.c_str()); }
int      get_int   (const string & s) { return s == "all" ? any  : atoi(s.c_str()); }
set<int> get_set   (const string & s) {
	set<int> res;
	vector<string> tokens = tokenize(s,",");
	for(int i = 0; i < tokens.size(); i++)
		res.insert(get_int(tokens[i]));
	return res;
}
bool isin(const set<int> & s, int i) { return s.find(i) != s.end(); }

int main(int argn, char ** argv)
{
	char * infile = "-", * outfile = "-";
	char * cutfile = 0;
	int ndiode = 4;
	vector<char*> args;
	for(char ** i = argv+1; *i; i++)
		if(false);
		else args.push_back(*i);
	if(args.empty()) { fprintf(stderr, "Expected [options] cutfile [runlist [outfile]]\n"); return 1; }
	cutfile = args[0];
	if(args.size() > 1) infile  = args[1];
	if(args.size() > 2) outfile = args[2];

	ifstream ifile; if(strcmp(infile,"-")) ifile.open(infile);
	istream & in = strcmp(infile,"-") ? ifile : cin;
	FILE * out = strcmp(outfile,"-") ? fopen(outfile,"w") : stdout;

	string line;
	// Now process the cut file
	ifstream cut(cutfile);
	string type, fromw, tow, runw, segw, modw, diw;
	vector<Cut> cuts;
	while(getline(cut,line))
	{
		stringstream ss(line);
		Cut c;
		ss >> type >> fromw >> tow >> runw >> segw >> modw >> diw;
		if(type[0] == '#') continue;
		c.from = get_double(fromw);
		if(c.from == dany) c.from = -inf;
		c.to   = get_double(tow);
		if(c.to == dany) c.from = inf;;
		c.run  = get_set(runw);
		c.seg  = get_set(segw);
		c.mod  = get_set(modw);
		c.di   = get_set(diw);
		cuts.push_back(c);
	}

	int run, seg, foo, nmod; double from, to;
	// Quick and dirty: We simply look at the indentation
	while(getline(in,line))
	{
		stringstream ss(line);
		switch(indent(line))
		{
			case 4: ss >> run >> from >> to >> foo >> nmod; break;
			case 6:
			{
				ss >> seg >> from >> to;
				set<diode> rejected;
				for(int i = 0; i < cuts.size(); i++)
				{
					if(cuts[i].to < from || cuts[i].from > to) continue;
					if(!isin(cuts[i].run,any) && !isin(cuts[i].run,run)) continue;
					if(!isin(cuts[i].seg,any) && !isin(cuts[i].seg,seg)) continue;
					// We passed the tests, so apply this cut to the segment
					for(set<int>::iterator m = cuts[i].mod.begin(); m != cuts[i].mod.end(); m++)
					for(int m1 = *m == any ? 0 : *m, m2 = *m == any ? nmod : *m+1; m1 < m2; m1++)
					for(set<int>::iterator d = cuts[i].di.begin();  d != cuts[i].di.end();  d++)
					for(int d1 = *d == any ? 0 : *d, d2 = *d == any ? ndiode : *d+1; d1 < d2; d1++)
						rejected.insert(diode(m1,d1));
				}
				if(rejected == all_diodes(nmod, ndiode)) continue;
				// OK. Some diodes survived. Output the run.
				fprintf(out, "%4d %2d %3d", run, seg, rejected.size());
				for(set<diode>::iterator i = rejected.begin(); i != rejected.end(); i++)
					fprintf(out, " %d %d", i->first, i->second);
				fprintf(out,"\n");
			}
		}
	}
}
