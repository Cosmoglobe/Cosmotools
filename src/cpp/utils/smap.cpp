#include <handy.h>
#include <stoken.h>
#include <cstring>
using namespace std;
using namespace skn;

bool starts_with(const string & s, const string & target, int n)
{
	if(s.size() < n) return false;
	for(int i = 0; i < target.size() && i < s.size(); i++)
		if(target[i] != s[i]) return false;
	return true;
}

Pixset prune_edge(const SMap & smap, double radius)
{
	int n = int(radius/smap.pixsize());
	map<int,Pixset> to_remove = smap.layers(n);
	Pixset sparsity = smap.get_sparsity();
	for(int i = 0; i <= n; i++)
		sparsity -= to_remove[i];
	return sparsity;
}

Pixset grow_edge(const SMap & smap, double radius)
{
	int n = int(radius/smap.pixsize());
	map<int,Pixset> to_add = smap.layers(0,-n);
	Pixset sparsity = smap.get_sparsity();
	for(int i = 0; i >= -n; i--)
		sparsity += to_add[i];
	return sparsity;
}

int main(int argn, char ** argv)
{
	try
	{
		vector<char*> args;
		bool smooth = false;
		for(char ** i = argv+1; *i; i++)
			if(!strcmp(*i,"-i")) smooth ^= true;
			else args.push_back(*i);
		if(args.size() < 2) serror(
			"Syntax: smap sparse   file [ofile]: Make file sparse.\n"
			"        smap full     file [ofile]: Make file dense.\n"
			"        smap prune    file [val|range [ofile]]: Remove pixels by value.\n"
			"        smap transfer file ofile: Give ofile file's pixel set.\n"
			"        smap intersection file1 file2 [ofile1 [ofile2]. If\n"
			"          outfiles not given, overwrite file1 and file2 with\n"
			"          their intersection. Otherwise, file1's part of the\n"
			"          intersection goes into ofile1, and so on..\n"
			"        smap crop     file arcmins [ofile]: Remove arcmins\n"
			"        smap nside    file nside [ofile]: Resize map.\n"
			"If ofile is not specified, the result is generally written back\n"
			"to the input file."
			);

		char * command = args[0], * filename = args[1];

		vector<SMap> maps = fits2smap<double>(filename);
		if(starts_with(command, "sparse", 1))
		{
			prune(maps);
			unify(maps);
			char * ofilename = args.size() > 2 ? args[2] : filename;
			smap2fits(ofilename, maps);
		}
		else if(starts_with(command, "fullsky", 1))
		{
			char * ofilename = args.size() > 2 ? args[2] : filename;
			map2fits(ofilename, smap2map(maps));
		}
		else if(starts_with(command, "prune", 1))
		{
			vector<string> range = tokenize(args[2],": ");
			if(range.size() == 1) prune(maps, atof(range[0].c_str()));
			else if(range.size() == 2) prune(maps, atof(range[0].c_str()), atof(range[1].c_str()));
			else serror("Prune expected a range of 1 or 2 elements, but got %d!", range.size());
			unify(maps);
			char * ofilename = args.size() > 3 ? args[3] : filename;
			smap2fits(ofilename, maps);
		}
		else if(starts_with(command, "transfer", 1))
		{
			if(args.size() < 3) serror("Output file needed for transfer");
			char * ofilename = args[2];
			vector<SMap> omaps = fits2smap<double>(ofilename);
			//prune(maps);
			//unify(maps);
			for(int i = 0; i < omaps.size(); i++)
			{
				if(omaps[i].Scheme() != maps[i].Scheme())
					maps[i].swap_scheme();
				omaps[i].set_sparsity(maps[i].get_sparsity());
			}
			smap2fits(ofilename, omaps);
		}
		else if(starts_with(command, "intersection", 1))
		{
			if(args.size() < 3) serror("Intersection: map1 map2 [omap1 [omap2]].\nInput maps are overwritten with intersection if no omaps are present.\nOtherwise, map1's part of the intersection goes into omap1, and similarly for omap2\nif present");
			char * filename2 = args[2];
			vector<SMap> maps2 = fits2smap<double>(filename2);
			//prune(maps); prune(maps2);
			//unify(maps); unify(maps2);
			Pixset common = getpix(maps) & getpix(maps2);
			setpix(maps,common); setpix(maps2,common);
			if(args.size() == 3)
			{
				smap2fits(filename,maps);
				smap2fits(filename2,maps2);
			}
			else
			{
				smap2fits(args[3],maps);
				if(args.size() > 4) smap2fits(args[4],maps2);
			}
		}
		else if(starts_with(command, "remove", 1))
		{
			if(args.size() < 3) serror("remove: map1 map2 [omap1].\nRemove the pixels in map2 from those in map1. map1 is overwritten if omap is absent.");
			char * filename2 = args[2];
			vector<SMap> maps2 = fits2smap<double>(filename2);
			//prune(maps); prune(maps2);
			//unify(maps); unify(maps2);
			Pixset res = getpix(maps); res -= getpix(maps2);
fprintf(stderr, "%d %d %d\n", getpix(maps).size(), getpix(maps2).size(), res.size());
			setpix(maps,res);
			smap2fits(args[args.size() > 3 ? 3 : 1],maps);
		}
		else if(starts_with(command, "crop", 1))
		{
			if(args.size() < 3) serror("smap crop file arcmins [ofile]");
			double radius = atof(args[2])*M_PI/180/60;
			for(int i = 0; i < maps.size(); i++)
			{
				Pixset pix = radius >= 0 ? prune_edge(maps[i], radius) : grow_edge(maps[i], -radius);
				maps[i].set_sparsity(pix);
			}
			char * ofilename = args.size() > 3 ? args[3] : filename;
			smap2fits(ofilename, maps);
		}
		else if(starts_with(command, "nside", 1))
		{
			if(args.size() < 3) serror("smap nside file nside [ofile]");
			//prune(maps);
			vector<SMap> out = maps;
			int nside = atoi(args[2]);
			if(!smooth)
				for(int i = 0; i < out.size(); i++)
					out[i].SetNside(nside);
			else
			{
				Healpix_Ordering_Scheme scheme = out[0].Scheme();
				for(int i = 0; i < out.size(); i++)
					out[i].SetNside(nside, RING);
				for(int i = 0; i < out.size(); i++)
				for(int j = 0; j < out[i].size(); j++)
				{
					int k = out[i].index_to_pixel(j);
					pointing p = out[i].pix2ang(k);
					out[i](j) = maps[i].interpolated_value(p);
				}
				for(int i = 0; i < out.size(); i++)
					out[i].SetNside(nside, scheme);
			}
			unify(out);
			char * ofilename = args.size() > 3 ? args[3] : filename;
			smap2fits(ofilename, out);
		}
		else serror("Unrecognized command %s!\n", command);
	} catch(const Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); }
}
