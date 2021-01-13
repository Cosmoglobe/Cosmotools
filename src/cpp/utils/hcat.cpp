/* Concatenate slices of hdf and fits maps. Unlike mapcat, this program
 * will not try to be a swiss army knife - it will just slice and
 * concatenate */

#include <cstdlib>
#include <cstdio>
#include <stimer.h>
#include <limits>
#include <hmap.h>
#include <hmap_io.h>
#include <cstring>
using namespace std;
using namespace skn;

struct Ifile
{
	string filename;
	Slice slice;
};

typedef HMap<double> HdfMap;

void help()
{
	fprintf(stderr,
"hcat: Slice and concatenate healpix maps in fits and hdf format.\n"
"Usage: hcat {options} ifile1{[slice]} {ifile2{[slice]} {...}} ofile\n"
"       Here {} indicates optional sections.\n"
"Options:\n"
"  -v: Increase verbosity.\n"
"  -q: Decrease verbosity.\n"
"  -d: Dimension to concatenate along.\n"
"      Dimensions are 0=map, 1=component and 2=pixel.\n"
"      Default: 2\n"
"  -float: Write output in single precision.\n"
"  -double: Write output in double precision (default).\n"
"  -h: Display this help message.\n"
"Slice syntax:\n"
"  [{from}{:{num}{:{step}}}]\n"
"  I.e. like fortran, or scipy arrays, but the second element\n"
"  is the number of elements, not the end position.\n"
"  [5:5:2] = 5,7,9,11,13\n"
"Examples:\n"
"  Convert from fits to hdf:\n"
"    hcat a.fits b.hdf\n"
"  Concatenate several fits files into one hdf file:\n"
"    hcat foo*.fits bar.hdf\n"
"  Extract the first map in a hdf map to another hdf file:\n"
"    hcat bar.hdf[0] baz.hdf\n"
"  Extract the fifth map and second component into one fits file:\n"
"    hcat bar.hdf[4,1] single.fits\n"
"  Create TQU fits map from separate T, Q and U maps:\n"
"    hcat -d 1 T.fits Q.fits U.fits TQU.fits\n"
"  Extract all maps in hdf file to fits files:\n"
"    hcat bar.hdf map%%02d.fits\n"
"    This produces map00.fits, map01.fits and so on.\n"
"Remarks:\n"
"  Operations in the pixel direction are nonintuitive and\n"
"    probably best avoided. I.e. don't specify -d 2.\n"
"  Maps being concatenated must be compatible.\n"
"");
	exit(1);
}

int main(int argc, char ** argv)
{
	try { 
		vector<char*> args;
		vector<int> hdu;
		double t1, t2;
		int verbosity = 0;
		int dimension = 0;
		int dhdu = 2, chdu = dhdu;
		char * type = 0;
		for(char ** i = argv+1; *i; i++)
			if(!strcmp(*i,"-float")) type = "float";
			else if(!strcmp(*i,"-double")) type = "double";
			else if(!strcmp(*i,"-v")) verbosity++;
			else if(!strcmp(*i,"-q")) verbosity--;
			else if(!strcmp(*i,"-hdu")) chdu = atoi(*++i);
			else if(!strcmp(*i,"-d")) dimension = atoi(*++i);
			else if(!strcmp(*i,"-h")) help();
			else if(**i == '-') serror("Unrecognized option: %s", *i);
			else {
				args.push_back(*i);
				hdu.push_back(chdu);
				chdu = dhdu;
			}
		if(args.size() < 2) help();

		string ofname = args.back();
		args.pop_back();
		vector<char*> ifnames = args;

		// Parse all the input args, turning them into Ifiles
		vector<Ifile> ifiles;
		for(int i = 0; i < ifnames.size(); i++)
		{
			Ifile ifile;
			string str = ifnames[i];
			if(str[str.size()-1] != ']') // No slice
				ifile.filename = str;
			else
			{
				int j = str.find_last_of('[');
				if(j == string::npos) serror("Unmatched [] in slice!");
				ifile.filename = str.substr(0,j);
				ifile.slice = parse_slice(str.substr(j+1,str.size()-j-2));
			}
			ifiles.push_back(ifile);
		}

		// Read in everything, making a HdfMap for each input argument
		t1 = wall_time();
		vector<HdfMap> maps(ifiles.size());
		for(int i = 0; i < ifiles.size(); i++)
		{
			double t3 = wall_time();
			HdfMap map = read_hmap<double>(ifiles[i].filename,0,hdu[i]);
			maps[i] = map(ifiles[i].slice);
			double t4 = wall_time();
			if(verbosity > 1) fprintf(stderr, "Read %s: %7.3lf\n", ifnames[i], t4-t3);
		}
		t2 = wall_time();
		if(verbosity > 0) fprintf(stderr, "Read maps: %7.3lf\n", t2-t1);

		// Concatenate. May support other directions later.
		t1 = wall_time();
		HdfMap full = concatenate(maps, dimension);
		t2 = wall_time();
		if(verbosity > 0) fprintf(stderr, "Concatenated maps: %7.3lf\n", t2-t1);

		// And output
		t1 = wall_time();
		write_hmap(ofname, full, 0, type);
		t2 = wall_time();
		if(verbosity > 0) fprintf(stderr, "Wrote %s: %7.3lf\n", ofname.c_str(), t2-t1);
	} catch(Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); }
}
