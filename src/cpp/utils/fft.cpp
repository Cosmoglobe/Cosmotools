#include <complex>
#include <fftw3.h>
#include "tabcommon.h"

using namespace std;
using namespace skn;

vector<complex<double> > fft(vector<double> & data)
{
	int n = data.size(), m = data.size()/2+1;
	vector<complex<double> > tmp(m);
	fftw_execute(fftw_plan_dft_r2c_1d(n, &data[0], (fftw_complex*)&tmp[0], FFTW_ESTIMATE));
	for(int i = 0; i < m; i++)
		tmp[i] /= sqrt(n);
	return tmp;
}

vector<complex<double> > fftc(vector<complex<double> > & data, int dir = 1)
{
	int n = data.size();
	vector<complex<double> > tmp(n);
	fftw_execute(fftw_plan_dft_1d(n, (fftw_complex*)&data[0], (fftw_complex*)&tmp[0], dir < 0 ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_ESTIMATE));
	for(int i = 0; i < n; i++)
		tmp[i] /= sqrt(n);
	return tmp;
}

vector<double> ffttime(vector<double> & time)
{
	vector<double> res(time.size()/2+1);
	double dt = (time.back()-time.front())/(time.size()-1);
	for(int i = 0; i < res.size(); i++)
		res[i] = i/dt/time.size();
	return res;
}

vector<double> fft_inv(vector<complex<double> > & data)
{
	int m = data.size(), n = (data.size()-1)*2;
	vector<double> tmp(n);
	fftw_execute(fftw_plan_dft_c2r_1d(n, (fftw_complex*)&data[0],&tmp[0], FFTW_ESTIMATE));
	for(int i = 0; i < n; i++)
		tmp[i] /= sqrt(n);
	return tmp;
}
vector<double> ffttime_inv(vector<double> & freq)
{
	vector<double> res((freq.size()-1)*2);
	// We have to choose whether we want odd or even original array lengths
	// to work. This is because the information about the odd/evenness of the
	// original array is lost after the transformation. I have chosen to prefer
	// even lengths
	double maxtime = 1.0/freq[1]*(res.size()-1)/res.size();
	for(int i = 0; i < res.size(); i++)
		res[i] = i*maxtime/(res.size()-1);
	return res;
}

int main(int argn, char ** argv)
{
	try {
		vector<char*> args;
		int timecol = 0; // which column is considered time
		bool inverse = false;
		char *ifilename = 0, *ofilename = 0;
		vector<int> cols, tcol;
		for(char ** i = argv+1; *i; i++)
			if(!strcmp(*i,"-t")) timecol = atoi(*++i)-1;
			else if(!strcmp(*i,"-i")) inverse ^= 1;
			else if(!(tcol = parse_sel(*i)).empty()) cols = tcol;
			else args.push_back(*i);

		if(args.size() > 0) ifilename = args[0];
		if(args.size() > 1) ofilename = args[1];
		
		vector<vector<double> > table = read_cols(ifilename, cols);
		// Make sure cols has the right size
		if(inverse && table.size() % 2 == !(timecol < table.size())) {
			cols.push_back(0);
			table.push_back(vector<double>(table.front().size()));
		}

		// Fourier transform. A bit complicated due to the expansion of complex numbers to columns
		vector<vector<double> > output;
		if(!inverse)
		{
			for(int i = 0; i < cols.size(); i++)
			{
				output.push_back(vector<double>());
				if(cols[i]-1 != timecol) output.push_back(vector<double>());
			}
			for(int i = 0, j = 0; i < cols.size(); i++, j++)
				if(cols[i]-1 == timecol) output[j] = ffttime(table[i]);
				else {
					vector<complex<double> > res = fft(table[i]);
					for(int k = 0; k < res.size(); k++)
					{
						output[j].push_back(res[k].real());
						output[j+1].push_back(res[k].imag());
					}
					j++;
				}
		}
		else
		{
			for(int i = 0, count = 0; i < cols.size(); i++)
				if(cols[i]-1 == timecol || !(count++%2))output.push_back(vector<double>());
			for(int i = 0, j = 0; i < output.size(); i++, j++)
			{
				if(cols[j]-1 == timecol) output[i] = ffttime_inv(table[j]);
				else
				{
					vector<complex<double> > tmp(table[j].size());
					for(int k = 0; k < tmp.size(); k++)
						tmp[k] = complex<double>(table[j][k],table[j+1][k]);
					output[i] = fft_inv(tmp);
					j++;
				}
			}
		}
		// Output
		FILE * ofile = ofilename && strcmp(ofilename,"-") ? fopen(ofilename,"w") : stdout;
		for(int i = 0; i < output[0].size(); i++)
		{
			for(int j = 0; j < output.size(); j++)
				fprintf(ofile, "%22.14le", output[j][i]);
			fprintf(ofile, "\n");
		}
	} catch(Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); }
}
