/* Dump fiels from level2- and level3-files to ascii format. */
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cstring>
#include <stoken.h>
#include <sgrid.h>
#include <serror.h>
#include <hdf5.h>
using namespace std;
using namespace skn;

typedef vector<vector<int> > Sel;
Sel parse_sel(const string & s)
{
	Sel res;
	vector<string> fields = tokenize(s, ",");
	for(int i = 0; i < fields.size(); i++)
	{
		vector<string> entries = tokenize_full(fields[i], ":");
		if(entries.size() > 3) serror("Unrecognized slice %s!", fields[i].c_str());
		vector<int> k;
		for(int j = 0; j < entries.size(); j++)
			k.push_back(entries[j].empty() ? -1 : atoi(entries[j].c_str()));
		res.push_back(k);
	}
	return res;
}

void help()
{
	fprintf(stderr, ""
		"lxtxt: Dump level2- and level3 files (and other hdf files) as an ascii table.\n"
		"Syntax: lxtxt input.hdf field[slice] field[slice] ...\n"
		"  Any number of fields can be specified. For each field, an optional\n"
		"  slice can be specified, in the format offset:num:step for each dimension.\n"
		"Examples:\n"
		"  lxtxt a.hdf time\n"
		"    Prints the time column of a.hdf.\n"
		"  lxtxt a.hdf time tod[3].\n"
		"    Prints two columns: The time column and the index 3 column in the tod.\n"
		"  lxtxt a.hdf tod[::4,::1000]\n"
		"    Prints every 1000th sample for every fourth diode in tod.\n");
}

int main(int argc, char ** argv)
{
	try {
		vector<char*> args;
		char * format = "%21.14le";
		bool hack_sec = true, decimate = false;
		for(argv++; *argv; argv++)
			if(!strcmp(*argv,"-fmt")) format = *++argv;
			else if(!strcmp(*argv,"-h")) { help(); return 0; }
			else if(!strcmp(*argv,"-d")) decimate ^= true;
			else if(strlen(*argv) > 1 && **argv == '-')
				return fprintf(stderr, "Unsupported option %s!\n", *argv), 1;
			else break;
		char fmt[0x100];
		sprintf(fmt, " %s", format);
		for(; *argv; argv++) args.push_back(*argv);

		if(args.empty())
			return fprintf(stderr, "Syntax: lxtxt lxfile.hdf [fields]\n"), 1;

		char * filename = args[0];
		if(args.size() == 1) return 0; // Nothing to do

		hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

		vector<Grid<double,2> > data;
		for(int i = 1; i < args.size(); i++)
		{
			string s = args[i];
			Sel sel;
			if(s[s.size()-1] == ']')
			{
				int j = s.find_last_of('[');
				sel = parse_sel(s.substr(j+1, s.size()-j-2));
				s = s.substr(0,j);
			}

			/* Hack: Artificial data set "sec" or "seconds". Maps to "time", but
			 * does some postprocessing. */
			bool do_sec = false;
			if(hack_sec && s == "sec" || s == "seconds") { do_sec = true; s = "time"; }

			hid_t set = H5Dopen(file, s.c_str(), H5P_DEFAULT);
			hid_t space = H5Dget_space(set);
			int dim = H5Sget_simple_extent_ndims(space);
			vector<hsize_t> dims(dim);
			H5Sget_simple_extent_dims(space, &dims[0], NULL);
			if(dim > 2 && sel.empty()) serror("Needs slice for %d-dimensional set %s!", dim, s.c_str());
			// If we have a slice, use it on the DataSpace
			vector<hsize_t> efflen = dims;
			vector<hsize_t> offset(dim), len(dim), stride(dim);
			if(!sel.empty())
			{
				sel.resize(dim);

				for(int j = 0; j < sel.size(); j++)
				{
					if(!sel[j].empty())
					{
						if(sel[j].size() == 1)
						{
							offset[j] = sel[j][0] < 0 ? 0 : sel[j][0];
							stride[j] = 1;
							len[j] = 1;
						}
						else
						{
							offset[j] = sel[j][0] < 0 ? 0 : sel[j][0];
							if(sel[j].size() > 2)
								stride[j] = sel[j][2] <= 0 ? 1 : sel[j][2];
							else
								stride[j] = 1;
							len[j] = sel[j][1] < 0 ? (dims[j]-offset[j])/stride[j] : sel[j][1];
						}
					}
					else
					{
						offset[j] = 0;
						stride[j] = 1;
						len[j] = dims[j];
					}
				}
				/* If we want decimating strides, we will read in all the data, and
				 * decimate later. */
				if(decimate)
				{
					vector<hsize_t> stride_(dim,1), len_(dim);
					for(int i = 0; i < dim; i++) len_[i] = len[i]*stride[i];
					H5Sselect_hyperslab(space, H5S_SELECT_SET, &offset[0], &stride_[0], &len_[0], NULL);
				}
				else
					H5Sselect_hyperslab(space, H5S_SELECT_SET, &offset[0], &stride[0], &len[0], NULL);
				efflen = len;
			}
			else
			{
				offset = vector<hsize_t>(dim,0);
				stride = vector<hsize_t>(dim,1);
				len    = dims;
			}
			// Filter out dimensions with length 1
			vector<hsize_t> memlen, memstride;
			for(int i = 0; i < efflen.size(); i++)
				if(efflen[i] > 1) {
					memlen.push_back(efflen[i]);
					memstride.push_back(stride[i]);
				}
			memlen.resize(2);
			memstride.resize(2);
			for(int i = 0; i < memlen.size(); i++)
				if(memlen[i] < 1) {
					memlen[i] = 1;
					memstride[i] = 1;
				}

			// If the result is one-dimensional, make sure the longest
			// dimension is the last one
			if(memlen[1] == 1) {
				swap(memlen[0], memlen[1]);
				swap(memstride[0], memstride[1]);
			}

			/* Now actually read in the data. Decimating strides makes this
			 * uglier. */
			Grid<double,2> g(memlen[0], memlen[1]);
			if(!decimate)
			{
				hid_t mem = H5Screate_simple(2, &memlen[0], NULL);
				H5Dread(set, H5T_NATIVE_DOUBLE, mem, space, H5P_DEFAULT, &g[0]);
			}
			else
			{
				vector<hsize_t> readlen = memlen;
				for(int i = 0; i < 2; i++) readlen[i] *= memstride[i];
				Grid<double,2> full(readlen[0], readlen[1]);
				hid_t mem = H5Screate_simple(2, &readlen[0], NULL);
				H5Dread(set, H5T_NATIVE_DOUBLE, mem, space, H5P_DEFAULT, &full[0]);
				/* The actual decimation */
				for(int i = 0; i < readlen[0]; i++)
				for(int j = 0; j < readlen[1]; j++)
					g(i/memstride[0],j/memstride[1]) += full(i,j);
				g /= memstride[0]*memstride[1];
			}

			/* Hack: Support fake seconds stream */
			if(do_sec) g = (g-g[0])*(24*60*60);

			data.push_back(g);
		}

		// Check the compatibility of the data sets
		int n = data[0].size(1);
		for(int i = 1; i < data.size(); i++)
		{
			//if(data[i].size(1) != n) serror("out[%d] has size %d, but out[0] has size %d!", i, data[i].size(1), n);
			n = min(n, data[i].size(1));
		}

		// Finally, output the data with rows as columns
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < data.size(); j++)
			for(int k = 0; k < data[j].size(0); k++)
				printf(fmt, data[j](k,i));
			printf("\n");
		}
	} catch(Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); return 1; }
}
