#include <handy.h>
#include <hmap.h>
#include <hmap_io.h>

using namespace std;
using namespace skn;

int main(int argc, char ** argv)
{
	try {
		vector<char*> args;
		int ncomp = 0, which = 0, lmax = 0;
		for(char ** i = argv+1; *i; i++)
			if(!strcmp("-map",*i)) which = atoi(*i++)-1;
			else if(!strcmp("-lmax",*i)) lmax = atoi(*i++);
			else if(!strcmp("-ncomp",*i)) ncomp = atoi(*i++);
			else if(**i == '-' && strlen(*i) > 1) serror("Unknown option %s!", *i);
			else args.push_back(*i);
		if(args.size() < 2) serror("Usage: map2ps [options] map ofile");

		HMap<double> imaps = read_hmap<double>(args[0]);
		int nside = imaps.get_nside(), npix = 12*nside*nside;
		if(!ncomp) ncomp = imaps.size(1) < 3 ? 1 : 3;
		if(!lmax)  lmax  = 2*nside;
		vector<Map> maps(ncomp, Map(nside, imaps.get_ordering() ? RING : NEST, nside_dummy()));
		for(int c = 0; c < ncomp; c++)
			for(int i = 0; i < nside; i++)
				maps[c][i] = imaps.pix(which, c, i);
		vector<Almc> alms = map2alm(maps, lmax);
		vector<PowSpec> ps = alm2powspec(alms);
		powspec2file(ps, args[1]);
	} catch(const Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); }
}
