#include <seph.h>
#include <algorithm>
#include <map>
#include <fstream>
#include <serror.h>
#include <cstring>

using namespace std;
using namespace skn;

/* Given a ces list and file range list (both output by ces_detect), produce a
 * properly formatted level1-runlist. */

struct File { Vector<double,2> mjd; string name; };
struct CES
{
	Vector<double,2> mjd;
	string object;
	vector<File> files;
	double az, el, dk, lon, lat;
	int cid;
};
bool operator<(const File & a, const File & b) { return a.mjd(0) < b.mjd(0); }
bool operator<(const CES & a,  const CES & b)  { return a.mjd(0) < b.mjd(0); }

int main(int argn, char ** argv)
{
	try {
		vector<char*> args;
		char * database = "/projects/quiet/external_data/unix.405",
			* patches = "/projects/quiet/auxilliary/patches_galactic.txt";
		bool own_object = false;
		for(char ** i = argv+1; *i; i++)
			if(!strcmp("-db",*i)) database = *++i;
			else if(!strcmp("-patches",*i)) patches = *++i;
			else if(!strcmp("-o",*i)) own_object ^= true;
			else args.push_back(*i);

		vector<Object> objdb = read_objs(patches);
		Ephfile eph(database);

		if(args.size() != 2) serror("Usage: ces2runlist ces_list.txt range_list.txt");
		char * ces_file = args[0], * range_file = args[1];

		// Read in the basic CES information
		vector<CES> ceses;
		ifstream ices(ces_file);
		CES ces;
		while(ices >> ces.mjd(0) >> ces.mjd(1) >> ces.az >> ces.el >>
			ces.dk >> ces.lon >> ces.lat >> ces.object)
			ceses.push_back(ces);
		sort(ceses.begin(), ceses.end());

		// Read in the file range information
		vector<File> files_tmp;
		ifstream irange(range_file);
		File f;
		while(irange >> f.mjd(0) >> f.mjd(1) >> f.name) files_tmp.push_back(f);
		sort(files_tmp.begin(), files_tmp.end());
		// Only keep unique files
		vector<File> files;
		if(!files_tmp.empty()) files.push_back(files_tmp[0]);
		for(int i = 1; i < files_tmp.size(); i++)
			if(files_tmp[i].name != files_tmp[i-1].name)
				files.push_back(files_tmp[i]);

		// Will need these for binary search
		vector<double> file_from, file_to;
		for(int i = 0; i < files.size(); i++)
		{
			file_from.push_back(files[i].mjd(0));
			file_to.push_back(files[i].mjd(1));
		}

		// Now loop through each CES, determine which object it belongs to,
		// and which files contain it.
		for(int i = 0; i < ceses.size(); i++)
		{
			CES & c = ceses[i];
			if(own_object) c.object = identify(eph, objdb, (c.mjd(0)+c.mjd(1))/2, c.lon*M_PI/180, c.lat*M_PI/180);
			vector<double>::iterator
				from = upper_bound(file_to  .begin(), file_to  .end(), c.mjd(0)),
				to   = upper_bound(file_from.begin(), file_from.end(), c.mjd(1));
			int fi = from-file_to.begin(), ft = to-file_from.begin();
			for(int j = fi; j < ft; j++)
				c.files.push_back(files[j]);
			c.cid = i+1;
		}

		// Group into classes
		map<string,vector<CES> > ces_db;
		for(int i = 0; i < ceses.size(); i++)
			ces_db[ceses[i].object].push_back(ceses[i]);

		// And finally output in runlist-format
		printf("%d\n", ces_db.size());
		for(map<string,vector<CES> >::iterator i = ces_db.begin(); i != ces_db.end(); i++)
		{
			vector<CES> & v = i->second;
			printf("%s %d\n", i->first.c_str(), v.size());
			for(int j = 0; j < v.size(); j++)
			{
				CES & c = v[j];
				printf("  %d %13.7lf %13.7lf %2d %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf\n",
					c.cid, c.mjd(0), c.mjd(1), c.files.size(), c.az, c.el, c.dk, c.lon, c.lat);
				for(int k = 0; k < c.files.size(); k++)
					printf("    %13.7lf %13.7lf %s\n", c.files[k].mjd(0), c.files[k].mjd(1),
						c.files[k].name.c_str());
			}
		}
	} catch(Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); }
}
