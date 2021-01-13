#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <serror.h>
#include <map>

using namespace std;
using namespace skn;

/* Given a set of CES lists and corresponding file lists, merge them by matching
 * up the CESes, and for ones that are in both, output the longest one. */

static const double inf = 1.0/0.0;
struct File { double mjd[2]; string name; };
struct CES
{
	double mjd[2];
	string object;
	vector<File> files;
	double az, el, dk, lon, lat;
	int cid;
};
istream & operator>>(istream & in, File & f) { return in >> f.mjd[0] >> f.mjd[1] >> f.name; }
istream & operator>>(istream & in, CES & c)
{ return in >> c.mjd[0] >> c.mjd[1] >> c.az >> c.el >> c.dk >> c.lon >> c.lat >> c.object;}
bool operator<(const File & a, const File & b) { return a.mjd[0] < b.mjd[0]; }
bool operator<(const CES & a,  const CES & b)  { return a.mjd[0] < b.mjd[0]; }

/* Take in sorted ces lists, and iterate together. These lists are
 * supposed to represent sections of the same, true CES range, but cut off
 * at different points. This means that a single CES only may overlap with
 * one single CES on the other side. */
vector<CES> merge(const vector<vector<CES> > & ceses)
{
	vector<CES> res;
	int n = ceses.size(), done = 0;
	vector<int> i(n);
	while(done < n)
	{
		// Find the earliest start from those still active.
		int b = 0; double v = inf;
		for(int j = 0; j < n; j++)
			if(i[j] < ceses[j].size() && ceses[j][i[j]].mjd[0] < v)
				{ b = j; v = ceses[j][i[j]].mjd[0]; }
		// Find overlaps
		vector<int> overlap;
		for(int j = 0; j < n; j++)
			if(i[j] < ceses[j].size() && ceses[j][i[j]].mjd[0] < ceses[b][i[b]].mjd[1])
				overlap.push_back(j);
		// Merge overlapping
		CES c = ceses[b][i[b]];
		for(int j = 0; j < overlap.size(); j++)
		{
			int k = overlap[j];
			if(ceses[k][i[k]].mjd[1] > c.mjd[1]) c.mjd[1] = ceses[k][i[k]].mjd[1];
		}
		res.push_back(c);
		// Advance all that overlapped
		for(int j = 0; j < overlap.size(); j++)
			i[overlap[j]]++;
		// Count done
		done = 0;
		for(int j = 0; j < n; j++) done += i[j] >= ceses[j].size();
	}
	return res;
}

vector<File> merge(const vector<vector<File> > & files)
{
	map<string,File> m;
	for(int i = 0; i < files.size(); i++)
	for(int j = 0; j < files[i].size(); j++)
		m[files[i][j].name] = files[i][j];
	vector<File> res;
	for(map<string,File>::iterator i = m.begin(); i != m.end(); i++)
		res.push_back(i->second);
	sort(res.begin(), res.end());
	return res;
}

int main(int argn, char ** argv)
{
	try {
		vector<char*> args;
		for(char ** i = argv+1; *i; i++)
			if(false);
			else args.push_back(*i);

		if(args.size() % 2 != 0) serror("ceslist_merge [ceses1.txt files1.txt [ceses2.txt files2.txt [...]]] ceslist_out.txt filelist_out.txt");

		char * ocesfile = args[args.size()-2], * ofilefile = args[args.size()-1];
		vector<char*> cesfiles, filefiles;
		for(int i = 0; i < args.size()-2; i+=2)
		{
			cesfiles.push_back(args[i]);
			filefiles.push_back(args[i+1]);
		}

		// Read all the CES lists and all the file lists
		vector<vector<CES> > ceses;
		vector<vector<File> > files;
		for(int i = 0; i < cesfiles.size(); i++)
		{
			vector<CES> clist; CES c;
			vector<File> flist; File f;
			ifstream cfile(cesfiles[i]), ffile(filefiles[i]);
			while(cfile >> c) clist.push_back(c);
			while(ffile >> f) flist.push_back(f);
			sort(clist.begin(), clist.end());
			ceses.push_back(clist);
			files.push_back(flist);
		}

		vector<CES> oces = merge(ceses);
		vector<File> ofiles = merge(files);

		FILE * out = fopen(ocesfile, "w");
		for(int i = 0; i < oces.size(); i++)
		{
			const CES & c = oces[i];
			fprintf(out, "%13.7lf %13.7lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %s\n",
				c.mjd[0], c.mjd[1], c.az, c.el, c.dk, c.lon, c.lat, c.object.c_str());
		}
		fclose(out);

		out = fopen(ofilefile, "w");
		for(int i = 0; i < ofiles.size(); i++)
			fprintf(out, "%13.7lf %13.7lf %s\n", ofiles[i].mjd[0], ofiles[i].mjd[1], ofiles[i].name.c_str());
		fclose(out);
	} catch(Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); }
}
