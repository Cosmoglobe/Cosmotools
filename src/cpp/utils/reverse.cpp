#include <cstdlib>
#include <cstdio>
#include <vector>
#include <map>
#include <spng.h>
#include <serror.h>
#include <hmap_io.h>
#include <pointing.h>
#include <healpix_map.h>

using namespace std;
using namespace skn;

typedef HMap<double> Map;

void help()
{
	fprintf(stderr, "Usage: reverse colorbar cmin cmax image lmin lmax bmin bmax nside omap\n");
}

/* read_png is broken on titan, so use this one */
Image read_ppm(const char * fname)
{
	FILE * f = fopen(fname, "r");
	int w, h, mval;
	fscanf(f, "P6 %d %d %d\n", &w, &h, &mval);
	Image res(w, h);
	for(int y = 0; y < res.size(1); y++)
	for(int x = 0; x < res.size(0); x++)
	{
		for(int j = 0; j < 3; j++)
			fread(&res(x,y)[j],1,1,f);
	}
	return res;
}

/* Simple and slow linear search */
void find_closest(const vector<Pixel> & cols, const Pixel & c, int & bi, double & bd)
{
	bd = 1.0/0.0;
	bi = 0;
	for(int i = 0; i < cols.size(); i++)
	{
		double d = 0, a;
		for(int j = 0; j < 3; j++) { a = (double)c[j]-cols[i][j]; d += a*a; };
		if(d < bd) { bd = d; bi = i; }
	}
	//if(bd < 64)
	//printf("%02x%02x%02x %02x%02x%02x %15.7e\n", c[0],c[1],c[2],cols[bi][0],cols[bi][1],cols[bi][2],bd);
}

int main(int argc, char ** argv)
{
	try {
		vector<char *> args;
		for(char ** i = argv+1; *i; i++)
			if(false);
			else args.push_back(*i);
		/* Arguments: colorbar cmin cmax image lmin lmax bmin bmax nside omap */
		if(args.size() != 10) { help(); return 1; }
		int dlim = 128;

		fprintf(stderr, "Reading..\n");
		Image colorbar = read_ppm(args[0]);
		double cmin = atof(args[1]), cmax = atof(args[2]);
		Image img = read_ppm(args[3]);
		double lmin = atof(args[4]), lmax = atof(args[5]),
			bmin = atof(args[6]), bmax = atof(args[7]);
		int nside = atoi(args[8]);
		char * ofile = args[9];

		write_png("foo.png", img);

		// With that out of the way, let us build the color lookup. Exact lookups
		// are fast with maps, but are fragile. Distance-based lookups are (much)
		// slower, but robust. There is also the possibility of reverse-engineering
		// the color formula. Most likely, a linear formula is in use, in which
		// case the derivative should be linewise constant. I will try direct lookup
		// for now.
		fprintf(stderr, "Building color map..\n");
		map<Pixel,double> color_table;
		map<Pixel,int>    color_counts;
		int cw = colorbar.size(0);
		for(int i = 0; i < cw; i++)
		{
			double v = (cmax-cmin)*i/cw;
			color_table[colorbar[i]] += v;
			color_counts[colorbar[i]]++;
		}
		map<Pixel,double>::iterator cti;
		map<Pixel,int>::iterator cci;
		vector<Pixel> colors;
		for(cti = color_table.begin(), cci = color_counts.begin(); cti != color_table.end(); cti++, cci++)
		{
			cti->second /= cci->second;
			colors.push_back(cti->first);
		}

		// The color table now maps each uniqe color to a value, and the colors vector
		// holds a lienar copy for fast access.

		// Now run through each pixel in the input map, calculate its coordinates, and
		// look up its value based on the color. Accumulate valid colors in healpix
		// pixels.
		fprintf(stderr, "Buliding healpix map..\n");
		Healpix_Base hpix(nside, RING, nside_dummy());
		Map omap(nside, 1, 1, 2);
		for(int x = 0; x < img.size(0); x++)
		{
			double l = lmin + (lmax-lmin)*x/img.size(0);
			for(int y = 0; y < img.size(1); y++)
			{
				int y2 = img.size(1)-y-1;
				double b = bmin + (bmax-bmin)*y2/img.size(1);
				pointing point(halfpi-b*pi/180, l*pi/180);
				int pix = hpix.ang2pix(point);
				map<Pixel,double>::const_iterator it = color_table.find(img(x,y));
				if(it != color_table.end())
				{
					omap.pix(0,0,pix) += it->second;
					omap.pix(1,0,pix) += 1;
				}
				else
				{
					int ind; double dist;
					find_closest(colors, img(x,y), ind, dist);
					if(dist < dlim)
					{
						omap.pix(0,0,pix) += color_table[colors[ind]];
						omap.pix(1,0,pix) += 1;
					}
				}
			}
		}
		for(int i = 0; i < omap.size(2); i++)
		{
			if(omap.pix(1,0,i) > 0) omap.pix(0,0,i) /= omap.pix(1,0,i);
			else omap.pix(0,0,i) = Healpix_undef;
		}

		write_hmap(ofile, omap);

	} catch(Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); return 1; }
	return 0;
}
