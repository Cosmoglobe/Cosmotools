#include <handy.h>
#include <hmap.h>
#include <hmap_io.h>
#include <stimer.h>
using namespace std;

int main(int argn, char ** argv)
{
	try {
		vector<char*> args;
		char * psfile = 0, * ofile = 0;
		int ns = 0x100, lmax = 0, ordering = DIAGONAL_FIRST, skipl = -1, nc;
		bool scaled = true, pixwin = false, cauchy = false;
		int verbosity = 0;
		char * beamarg = 0;
		double t1, t2;
		for(char ** i = argv+1; *i; i++)
			if(!strcmp("-n",*i) || !strcmp("-nside",*i)) ns = atoi(*++i);
			else if(!strcmp("-p",*i) || !strcmp("-pixwin",*i)) pixwin ^= true;
			else if(!strcmp("-b",*i)) beamarg = *++i;
			else if(!strcmp("-s",*i)) scaled ^= true;
			else if(!strcmp("-cauchy",*i)) cauchy ^= true;
			else if(!strcmp("-seed",*i)) rng = planck_rng(atoi(*++i));
			else if(!strcmp("-v",*i)) verbosity++;
			else if(!strcmp("-q",*i)) verbosity--;
			else if(**i == '-' && strlen(*i) > 1) serror("Unknown option %s!", *i);
			else args.push_back(*i);
		if(args.size() < 2) {
			fprintf(stderr, "ps2map [options] ps ofile\n"
				" The input power spectrum can be either .txt of .fits\n"
				" The spectrum is usually assumed to be scaled by l(l+1)/2pi\n"
				" To avoid this, specify the -s option.\n"
				" Options:\n"
				"  -n -nside: Specify nside of output map\n"
				"  -b: Specify beam as fwhm or beam file\n"
				"  -s: Toggle l-scaling of power spectrum.\n"
				"  -seed: Specify random number seed.\n"
				"  -v: Increase verbosity.\n"
				"  -q: Decrease verbosity.\n");
			exit(1);
		}
		psfile = args[0]; if(!strcmp(psfile,"-")) psfile = "/dev/stdin";
		ofile  = args[1];
		if(!lmax) lmax = ns * 3;

		string fname(psfile);
		bool fits = fname.size() >= 5 && fname.substr(fname.size()-5) == ".fits";

		// Generate our starting map, using the TT column of a file
		if(verbosity > 0) fprintf(stderr, "Reading power spectrum ");
		t1 = wall_time();
		vector<PowSpec> p = fits ?
			fits2powspec(psfile, lmax, skipl) :
			ascii2powspec(psfile,lmax, skipl, scaled);
		t2 = wall_time();
		if(verbosity > 0) fprintf(stderr, "%.6lf\n", t2-t1);

		if(verbosity > 0) fprintf(stderr, "Processing spectrum ");
		t1 = wall_time();
		// Smooth with a beam of the ordero f the reduced pixel size
		if(beamarg)
		{
			double beam;
			if(sscanf(beamarg, "%lf", &beam) == 1)
			{
				double sigma = beam*M_PI/180/60/sqrt(8*log(2));
				for(int i = 0; i < p.size(); i++)
					for(int l = 0; l < p[i].Lmax(); l++)
						p[i].tt(l) *= exp(-l*(l+1)*sigma*sigma);
			}
			else
			{
				fname = string(beamarg);
				bool fits = fname.size() >= 5 && fname.substr(fname.size()-5) == ".fits";
				vector<PowSpec> b = fits ?
					fits2powspec(beamarg, lmax, skipl) :
					ascii2powspec(beamarg,lmax, skipl, 0);
				for(int i = 0; i < p.size(); i++)
					for(int l = 0; l < p[i].Lmax(); l++)
					{
						double scale = l > b[0].Lmax() ? 0 : b[0].tt(l);
						p[i].tt(l) *= scale*scale;
					}
			}
		}
		if(pixwin)
		{
			const arr<double> pw = get_window(ns);
			for(int i = 0; i < p.size(); i++)
				for(int l = 0; l < p[i].Lmax(); l++)
					p[i].tt(l) *= l < pw.size() ? pw[l]*pw[l] : 0;
		}
		t2 = wall_time();
		if(verbosity > 0) fprintf(stderr, "%.6lf\n", t2-t1);

		// All the output made this 5 times as long as it was!
		// Might remove it later.
		if(verbosity > 0) fprintf(stderr, "Making alms ");
		t1 = wall_time();
		vector<Almc>    a = cauchy ? powspec2alm_cauchy(p, ordering, rng) : powspec2alm(p,ordering, rng);
		nc = a.size();
		t2 = wall_time();
		if(verbosity > 0) fprintf(stderr, "%.6lf\n", t2-t1);

		if(verbosity > 0) fprintf(stderr, "Making maps ");
		t1 = wall_time();
		vector<Map>     m = alm2map(a,ns,RING);
		a.clear();
		t2 = wall_time();
		if(verbosity > 0) fprintf(stderr, "%.6lf\n", t2-t1);

		if(verbosity > 0) fprintf(stderr, "Writing output ");
		HMap<double> om(ns, 1, m.size(), 1);
		for(int i = 0; i < m.size(); i++)
		{
			for(int j = 0; j < 12*ns*ns; j++)
				om(0,i,j) = m[i][j];
			m[i] = Map(1, RING); // free memory
		}

		t1 = wall_time();
		write_hmap(ofile, om);
		t2 = wall_time();
		if(verbosity > 0) fprintf(stderr, "%.6lf\n", t2-t1);
	} catch(const Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); }
}
