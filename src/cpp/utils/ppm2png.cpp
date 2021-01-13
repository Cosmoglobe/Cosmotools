#include <cstdlib>
#include <cstdio>
#include <spng.h>

using namespace std;
using namespace skn;

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

int main(int argc, char ** argv)
{
	Image img = read_ppm(argv[1]);
	write_png(argv[2], img);
	return 0;
}
