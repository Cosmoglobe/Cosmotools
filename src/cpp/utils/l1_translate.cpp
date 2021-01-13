// Translate between l1 files and run+seg.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iostream>

using namespace std;

struct CES { CES(int r = -1, int s = -1):run(r),seg(s) {} int run, seg; };
struct Range { Range(double r1 = 0, double r2 = 0):from(r1),to(r2) {} double from,to; };
bool operator<(const CES & a, const CES & b) { return a.run == b.run ? a.seg < b.seg : a.run < b.run; }
bool operator<(const Range & a, const Range & b) { return a.from < b.from; }

// strip out directories, prune space
string get_basename(const string & line)
{
	int a = line.find_last_of("/"),
		b = line.find_last_not_of(" ");
	if(a != string::npos) return line.substr(a+1, b-a);
	else return line.substr(0, b+1);
}

// We assume that no ranges in the table overlap.
template <class T> map<Range,T> overlapping_elements(const map<Range,T> & table, Range r)
{
	map<Range,T> res;
	typedef typename map<Range,T>::const_iterator it;
	it fstart = table.lower_bound(r);
	it bstart = fstart; --bstart;
	if(fstart != table.begin() && bstart->first.to > r.from)
		res[bstart->first] = bstart->second;
	for(it i = fstart; i->first.from < r.to; i++)
		res[i->first] = i->second;
	return res;
}

void translate(const CES & ces, map<CES,Range> & from,
	map<Range,string> & to)
{
	const Range & r = from[ces];
	map<Range, string> res = overlapping_elements(to, r);
	printf("%d %d %13.7lf %13.7lf => %d\n", ces.run, ces.seg,
		r.from, r.to, res.size());
	int i = 1;
	for(map<Range,string>::const_iterator it = res.begin(); it != res.end(); it++, i++)
		printf("  %d %13.7lf %13.7lf %s\n", i, it->first.from, it->first.to,
			it->second.c_str());
}

void translate(const string & s, map<string,Range> & from,
	map<Range,CES> & to)
{
	Range & r = from[s];
	map<Range, CES> res = overlapping_elements(to, r);
	printf("%s %13.7lf %13.7lf => %d\n", s.c_str(),
		r.from, r.to, res.size());
	int i = 1;
	for(map<Range,CES>::const_iterator it = res.begin(); it != res.end(); it++, i++)
		printf("  %d %13.7lf %13.7lf %d %d\n", i, it->first.from, it->first.to,
			it->second.run, it->second.seg);
}

void help()
{
	printf("Syntax: l1_translate [options] [run seg/l1 file]\n"
		"Options:\n"
		" -1: Specify an l1 database manually.\n"
		" -2: Specify an l2 database manually.\n"
		" -h: Display this message and exit.\n");
}

int get_indent(const string & s, char ind = ' ')
{ return s.find_first_not_of(ind); }

int main(int argc, char ** argv)
{
	vector<char*> args;
	char * l1filename = "/data4/quiet/runlist.txt",
		* l2filename = "/data4/quiet/runlist_l2.txt";
	bool reverse = false;
	bool verbose = false;
	for(char ** i = argv+1; *i; i++)
		if(!strcmp("-1",*i)) l1filename = *++i;
		else if(!strcmp("-2",*i)) l2filename = *++i;
		else if(!strcmp("-h",*i)) { help(); return 1; }
		else if(!strcmp("-v",*i)) verbose ^= true;
		else args.push_back(*i);

	// Read l1 db:
	map<string,Range> from_l1;
	map<Range,string> to_l1;
	ifstream l1file(l1filename);
	string line;
	while(getline(l1file,line))
	{
		stringstream ss; ss << line;
		if(get_indent(line) != 6) continue;
		int seq; Range r; string s;
		ss >> seq >> r.from >> r.to >> s;
		s = get_basename(s);
		from_l1[s] = r;
		to_l1[r] = s;
	}
	l1file.close();

	// Read oslo db
	map<CES,Range> from_l2;
	map<Range,CES> to_l2;
	ifstream l2file(l2filename);
	CES ces;
	while(getline(l2file,line))
	{
		stringstream ss; ss << line;
		int indent = get_indent(line);
		Range r;
		if(indent == 4) ss >> ces.run;
		else if(indent == 6)
		{
			ss >> ces.seg >> r.from >> r.to;
			from_l2[ces] = r;
			to_l2[r] = ces;
		}
	}
	l2file.close();

	if(args.empty())
	{
		string line;
		while(getline(cin,line))
		{
			stringstream ss; ss << line;
			CES ces;
			if(ss >> ces.run >> ces.seg) translate(ces, from_l2, to_l1);
			else translate(get_basename(line), from_l1, to_l2);
		}
	}
	else
	{
		if(args.size() == 1) translate(get_basename(args[0]), from_l1, to_l2);
		else if(args.size() == 2)
			translate(CES(atoi(args[0]), atoi(args[1])), from_l2, to_l1);
		else { help(); return 1; }
	}
}
