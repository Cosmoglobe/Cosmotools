#include <complex>
#include <set>
#include <map>
#include <fftw3.h>
#include "tabcommon.h"

using namespace std;
using namespace skn;

vector<double> powspec(vector<double> & data)
{
	int n = data.size(), m = data.size()/2+1;
	vector<complex<double> > tmp(m);
	fftw_execute(fftw_plan_dft_r2c_1d(n, &data[0], (fftw_complex*)&tmp[0], FFTW_ESTIMATE));
	vector<double> res(m);
	for(int i = 0; i < m; i++)
		res[i] = (tmp[i]*conj(tmp[i])).real()/n;
	return res;
}
vector<double> crosspow(vector<double> & a, vector<double> & b)
{
	int n = a.size(), m = a.size()/2+1;
	vector<complex<double> > ta(m), tb(m);
	fftw_execute(fftw_plan_dft_r2c_1d(n, &a[0], (fftw_complex*)&ta[0], FFTW_ESTIMATE));
	fftw_execute(fftw_plan_dft_r2c_1d(n, &b[0], (fftw_complex*)&tb[0], FFTW_ESTIMATE));
	vector<double> res(m);
	for(int i = 0; i < m; i++)
		res[i] = 0.5*(ta[i]*conj(tb[i])+conj(ta[i])*tb[i]).real()/
			(abs(ta[i])*abs(tb[i]));
	return res;
}
vector<double> powtime(vector<double> & time)
{
	vector<double> res(time.size()/2+1);
	double dt = (time.back()-time.front())/(time.size()-1);
	for(int i = 0; i < res.size(); i++)
		res[i] = i/dt/time.size();
	return res;
}

int main(int argn, char ** argv)
{
	try {
		vector<char*> args;
		int timecol = 0; // which column is considered time
		char *ifilename = 0, *ofilename = 0;
		vector<vector<int> > cols, tcol;
		for(char ** i = argv+1; *i; i++)
			if(!strcmp(*i,"-t")) timecol = atoi(*++i)-1;
			else if(!(tcol = parse_sel_cross(*i)).empty()) cols = tcol;
			else args.push_back(*i);

		if(args.size() > 0) ifilename = args[0];
		if(args.size() > 1) ofilename = args[1];

		// This is a bit clumsy. Determine which columns to read.
		set<int> toread;
		map<int,int> col_mapping;
		for(int i = 0; i < cols.size(); i++)
		for(int j = 0; j < cols[i].size(); j++)
			toread.insert(cols[i][j]);
		vector<int> icols;
		int k = 0;
		for(set<int>::iterator i = toread.begin(); i != toread.end(); i++)
		{
			icols.push_back(*i);
			col_mapping[*i] = k++;
		}

		vector<vector<double> > itable = read_cols(ifilename,icols);

		// Handle empty input list
		if(cols.empty())
		{
			for(int i = 0; i < icols.size(); i++)
			{
				cols.push_back(vector<int>(1,icols[i]));
				col_mapping[icols[i]] = icols[i]-1;
			}
		}

		vector<vector<double> > otable;
		// Normal powspec output
		for(int i = 0; i < cols.size(); i++)
		{
			if(cols[i].front()-1 == timecol)
				otable.push_back(powtime(itable[col_mapping[cols[i][0]]]));
			else if(cols[i].size() == 1)
				otable.push_back(powspec(itable[col_mapping[cols[i][0]]]));
			else
				otable.push_back(crosspow(
					itable[col_mapping[cols[i][0]]],
					itable[col_mapping[cols[i][1]]]));
		}

		// Output
		FILE * ofile = ofilename ? fopen(ofilename,"w") : stdout;
		for(int i = 0; i < otable[0].size(); i++)
		{
			for(int j = 0; j < otable.size(); j++)
				fprintf(ofile, "%22.14le", otable[j][i]);
			fprintf(ofile, "\n");
		}
	} catch(Error & e) { fprintf(stderr, e.msg.c_str()); }
}
