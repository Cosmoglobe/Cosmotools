#ifndef HANDYSTUFF
#define HANDYSTUFF

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <xcomplex.h>
//#include <cxxutils.h>
#include <healpix_data_io.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <alm.h>
#include <alm_healpix_tools.h>
#include <alm_powspec_tools.h>
#include <fitshandle.h>
#include <powspec_fitsio.h>
#include <powspec.h>
#include <planck_rng.h>
#include <hpixdata.h>
#include <serror.h>
#include <sparse_healpix_map.h>
#include <fits_types.h>
#include <limits>
#include <fitsio.h>

using namespace skn;

// Common definitions
typedef Healpix_Map<double> Map;
typedef Sparse_Healpix_Map<double> SMap;
typedef Sparse_Healpix_Map<float> SFMap;
typedef xcomplex<double> Complex;
typedef Alm<Complex> Almc;

// Conversions between Map and SMap
template <class T> Healpix_Map<T> smap2map(const Sparse_Healpix_Map<T> & smap)
{
	Healpix_Map<T> res(smap.Nside(), smap.Scheme(), nside_dummy());
	for(int i = 0; i < smap.Npix(); i++)
		res[i] = smap[i];
	return res;
}
template <class T> Sparse_Healpix_Map<T> map2smap(const Healpix_Map<T> & map)
{
	Sparse_Healpix_Map<T> res(map.Nside(), map.Scheme(), nside_dummy());
	for(int i = 0; i < map.Npix(); i++)
		res(i) = map[i];
	return res;
}
template <class T> std::vector<Healpix_Map<T> > smap2map(const std::vector<Sparse_Healpix_Map<T> > & smap) {
	std::vector<Healpix_Map<T> > res;
	for(int i = 0; i < smap.size(); i++) res.push_back(smap2map(smap[i]));
	return res;
}
template <class T> std::vector<Sparse_Healpix_Map<T> > map2smap(const std::vector<Healpix_Map<T> > & map) {
	std::vector<Sparse_Healpix_Map<T> > res;
	for(int i = 0; i < map.size(); i++) res.push_back(map2smap(map[i]));
	return res;
}

static const double NaN = std::numeric_limits<double>::quiet_NaN();

const char * COLNAMES[] = { "signal", "Q-pol", "U-pol", 0 };

// Get a random number from the operating system
unsigned int getseed()
{
	FILE * urandom = fopen("/dev/urandom", "r");
	unsigned int res;
	fread(&res, sizeof(unsigned int), 1, urandom);
	return res;
}
static planck_rng rng(getseed(),getseed(),getseed(),getseed());

int nside2order(int nside) { int i = -1; for(;nside;nside >>= 1, i++); return i; }
int nside_round(double ns) { return (int) pow(2,floor(log(ns)/log(2.0)+0.5)); }
int nside_floor(double ns) { return (int) pow(2,floor(log(ns)/log(2.0))); }
int nside_ceil(double ns) { return (int) pow(2,ceil(log(ns)/log(2.0))); }
int npix2nside(int npix) { return (int) sqrt(npix/12.0); }

 /***************************************************
 * Input/output for maps, power spectra etc.        *
 ***************************************************/

// Read in a power spectrum from a text file. There are two possible
// orderings: Permutative ordering or diagonals first.
// Permutative:
//    tt te tb ee eb bb
//    tt te tb ee bb
//    tt te ee bb
//    tt te ee
//    tt ee
//    tt
// Diagonals first:
//    tt ee bb te tb eb
//    tt ee bb te tb
//    tt ee bb te
//    tt ee te
//    tt ee
//    tt
// When the full set is not given, the lest important parts are removed
// first.
enum PowSpecComponent { psTT, psTE, psTB, psEE, psEB, psBB, ps00 };
enum PowSpecOrdering { PERMUTATIVE, DIAGONAL_FIRST };
PowSpecComponent index2psc(int i, int n, int mode)
{
	const PowSpecComponent tab[2][6][6] = { {
		{ psTT, ps00, ps00, ps00, ps00, ps00 },
		{ psTT, psEE, ps00, ps00, ps00, ps00 },
		{ psTT, psTE, psEE, ps00, ps00, ps00 },
		{ psTT, psTE, psEE, psBB, ps00, ps00 },
		{ psTT, psTE, psTB, psEE, psBB, ps00 },
		{ psTT, psTE, psTB, psEE, psEB, psBB } },
	{
		{ psTT, ps00, ps00, ps00, ps00, ps00 },
		{ psTT, psEE, ps00, ps00, ps00, ps00 },
		{ psTT, psEE, psTE, ps00, ps00, ps00 },
		{ psTT, psEE, psBB, psTE, ps00, ps00 },
		{ psTT, psEE, psBB, psTE, psTB, ps00 },
		{ psTT, psEE, psBB, psTE, psTB, psEB } } };
	return tab[mode][n-1][i];
}
double & pslookup(PowSpec & p, int l, int n, int i, int mode)
{
	switch(index2psc(i,n,mode))
	{
	case psTT: return p.tt(l); break;
	case psTE: return p.tg(l); break;
	case psTB: return p.tc(l); break;
	case psEE: return p.gg(l); break;
	case psEB: return p.gc(l); break;
	case psBB: return p.cc(l); break;
	default: serror("Unrecognized powspec component!"); break;
	}
	return p.tt(l);
}
double pslookup(const PowSpec & p, int l, int n, int i, int mode)
{ return pslookup((PowSpec&)p,l,n,i,mode); }
void zero_powspec(PowSpec & ps) {
	if(ps.Num_specs() >= 1)
		for(int l = 0; l <= ps.Lmax(); l++) ps.tt(l) = 0;
	if(ps.Num_specs() >= 4)
		for(int l = 0; l <= ps.Lmax(); l++) ps.tg(l) = ps.gg(l) = ps.cc(l) = 0;
	if(ps.Num_specs() >= 6)
		for(int l = 0; l <= ps.Lmax(); l++) ps.tc(l) = ps.gc(l) = 0;
}

void resize_powspec(PowSpec & p, int lmax)
{
	if(p.Lmax() == lmax) return;
	arr<double> tmp(lmax+1);
	tmp.fill(0);
	for(int l = 0; l <= p.Lmax() && l <= lmax; l++)
		tmp[l] = p.tt(l);
	p.Set(tmp);
}
void resize_alm(Almc & alm, int lmax)
{
	if(alm.Lmax() == lmax) return;
	Almc tmp(lmax, lmax);
	tmp.SetToZero();
	for(int l = 0; l <= alm.Lmax() && l <= lmax; l++)
		for(int m = 0; m <= alm.Mmax() && m <= l; m++)
			tmp(l,m) = alm(l,m);
	alm = tmp;
}

/* Some handy operations on Data variants */
/* TODO: Templatize these */
std::vector<Map> fits2map(fitshandle & f, int hdu = 2)
{
	f.goto_hdu(hdu);
	std::vector<Map> ms;
	for(int i = 1; i <= f.ncols(); i++)
	{
		Map tmp;
		read_Healpix_map_from_fits(f,tmp,i);
		ms.push_back(tmp);
	}
	return ms;
}
std::vector<Map> fits2map(char * fn, int hdu = 2)
{ fitshandle f; f.open(fn); return fits2map(f, hdu); }

template<class T>
std::vector<Sparse_Healpix_Map<T> > fits2smap(fitshandle & f, int hdu = 2)
{
	typedef Sparse_Healpix_Map<T> SMap;
	f.goto_hdu(hdu);
	std::string ordering; f.get_key("ORDERING", ordering);
	Healpix_Ordering_Scheme scheme = ordering=="RING" ? RING : NEST;
	int nside;       f.get_key("NSIDE", nside);
	int npix = 12*nside*nside;
	std::string indexschm;f.get_key("INDXSCHM", indexschm);
	std::vector<int> pixels;
	int col = 1;
	if(indexschm == "IMPLICIT")
		for(int i = 0; i < npix; i++) pixels.push_back(i);
	else
	{
		pixels.resize(f.nelems(col));
		f.read_column_raw(col, &pixels[0], f.nelems(col));
		col++;
	}
	std::vector<SMap> ms;
	std::vector<T> data(pixels.size());
	for(; col <= f.ncols(); col++)
	{
		f.read_column_raw(col, &data[0], pixels.size());
		ms.push_back(SMap(nside, pixels, scheme, nside_dummy(), data));
	}
	return ms;
}
template <class T>
std::vector<Sparse_Healpix_Map<T> > fits2smap(const char * fn, int hdu = 2)
{ fitshandle f(fn); return fits2smap<T>(f, hdu); }


/* This function needs ms to actually contain compatible maps */
template <class T>
void map2fits(fitshandle & f, const std::vector<Healpix_Map<T> > & ms,
	std::vector<int> cols = std::vector<int>())
{
	int i;
	if(!cols.size()) for(int i = 1; i <= ms.size(); i++) cols.push_back(i);
	arr<std::string> colname(cols.size());
	for(i = 0; i < cols.size() && COLNAMES[i]; i++)
		colname[i] = COLNAMES[i];
	for(; i < cols.size(); i++)
		colname[i] = "unknown";
	prepare_Healpix_fitsmap (f, ms[0], fitstype<T>(), colname);
	for(i = 0; i < cols.size(); i++)
		f.write_column(cols[i],ms[i].Map());
}
template <class T>
void map2fits(const char * fn, const std::vector<Healpix_Map<T> > & ms,
	std::vector<int> cols = std::vector<int>()) {
	fitshandle f; f.create(std::string("!") + fn);
	map2fits(f, ms, cols);
}

template <class T>
void smap2fits(fitshandle & f, const std::vector<Sparse_Healpix_Map<T> > & ms)
{
	typedef Sparse_Healpix_Map<T> Smap;
	int ncol; for(ncol = 0; COLNAMES[ncol]; ncol++);
	std::vector<fitscolumn> cs;
	cs.push_back(fitscolumn("INDEX","Pixel index",1,PLANCK_INT32));
	for(int i = 0; i < ms.size(); i++)
		cs.push_back(fitscolumn(i >= ncol ? "unknown" : COLNAMES[i], "unknown", 1, fitstype<T>()));
	f.insert_bintab(cs);

	f.set_key("PIXTYPE",std::string("HEALPIX"),"HEALPIX pixelisation");
	f.set_key("ORDERING",std::string(ms[0].Scheme()==RING?"RING":"NESTED"),"Pixel ordering scheme, either RING or NESTED");
	f.set_key("NSIDE", ms[0].Nside(), "Resolution parameter for HEALPIX");
	f.set_key("FIRSTPIX",0,"First pixel # (0 based)");
	f.set_key("LASTPIX",ms[0].Npix()-1,"Last pixel # (0 based)");
	f.set_key("INDXSCHM",std::string("EXPLICIT"),"Indexint: IMPLICIT or EXPLICIT");

	int n = ms[0].size();
	f.write_column_raw(1, &ms[0].get_sparsity()[0], n);
	for(int i = 0; i < ms.size(); i++)
		f.write_column_raw(2+i,&ms[i](0),n);
}
template <class T>
void smap2fits(const char * fn, const std::vector<Sparse_Healpix_Map<T> > & m) {
	fitshandle f; f.create(std::string("!") + fn);

	smap2fits(f, m);
}

// Translates between two views of power spectra: the array view
// represents a multi-component power spectrum as an array of
// one-component spectra, which is easy to deal with and manipulate.
// The other view, needed by certain functions uses all the components
// of PowSpec, but is harder to deal with, and limited.
std::vector<PowSpec> powspec2array(const PowSpec & ps, int mode, int n)
{
	std::vector<PowSpec> res;
	PowSpec tmp(1, ps.Lmax());
	for(int i = 0; i < n; i++)
	{
		for(int l = 0; l <= ps.Lmax(); l++)
			tmp.tt(l) = pslookup(ps, l, n, i, mode);
		res.push_back(tmp);
	}
	return res;
}
PowSpec array2powspec(const std::vector<PowSpec> & array, int mode)
{
	int n = array.size();
	int m = n > 1 ? (n > 4 ? 6 : 4) : 1;
	PowSpec res(m, array[0].Lmax());
	zero_powspec(res);
	for(int i = 0; i < n; i++)
	{
		for(int l = 0; l <= res.Lmax(); l++)
			pslookup(res, l, n, i, mode) = array[i].tt(l);
	}
	return res;
}

// This one could be extended to support arbitrary many components,
// but I'll not do that now.
std::vector<Almc> powspec2alm(const std::vector<PowSpec> & ps, int mode, planck_rng & rng)
{
	int lmax = ps[0].Lmax();
	std::vector<Almc> alms;
	if(ps.size() == 1) {
		Almc a(lmax, lmax);
		create_alm(ps[0], a, rng);
		alms.push_back(a);
	}
	else {
		PowSpec p = array2powspec(ps, mode);
		alms = std::vector<Almc>(3, Almc(lmax,lmax));
		create_alm_pol(p, alms[0], alms[1], alms[2], rng);
	}
	return alms;
}
std::vector<Almc> powspec2alm(const std::vector<PowSpec> & ps, int mode = PERMUTATIVE) { return powspec2alm(ps,mode,rng); }

double rand_cauchy(planck_rng & rng) { return tan(M_PI*(rng.rand_uni()-0.5)); }

std::vector<Almc> powspec2alm_cauchy(const std::vector<PowSpec> & ps, int mode, planck_rng & rng)
{
	int lmax = ps[0].Lmax();
	std::vector<Almc> alms;
	if(ps.size() == 1) {
		Almc a(lmax, lmax);
		for(int l = 0; l <= lmax; l++)
		{
			double s = sqrt(ps[0].tt(l));
			double s2 = s / sqrt(2.0);
			a(l,0).Set(rand_cauchy(rng)*s,0);
			for(int m = 1; m <= l; m++)
			{
				double r1 = rand_cauchy(rng), r2 = rand_cauchy(rng);
				a(l,m).Set(r1*s2,r2*s2);
			}
		}
		alms.push_back(a);
	}
	else {
		serror("Cauchy realization not implemented for polarization!");
	}
	return alms;
}

// Non-warning generating version of the Healpix one
PowSpec extract_crosspow (const Almc & a, const Almc & b)
{
	int lmax = a.Lmax(), mmax = a.Mmax();
	PowSpec res(1,lmax);
	for (int l = 0; l <= lmax; l++)
	{
		res.tt(l) = a(l,0).re*b(l,0).re;
		int lim = std::min(l,mmax);
		for (int m = 1; m <= lim; m++)
			res.tt(l) += 2 * (a(l,m)*conj(b(l,m))).re;
		res.tt(l) /= (2*l+1);
	}
	return res;
}

std::vector<PowSpec> alm2powspec(const std::vector<Almc> & alms, int mode = PERMUTATIVE, int nspec = 0)
{
	if(alms.size() != 1 && alms.size() != 3)
		serror("powspec2alm: Only 1 or 3 components supported, but got %d!", alms.size());
	if(nspec == 0) nspec = alms.size() > 1 ? 6 : 1;
	int ncalc = nspec <= 1 ? 1 : nspec <= 4 ? 4 : 6;
	int lmax = alms[0].Lmax();
	std::vector<PowSpec> ps;
	if(alms.size() == 1) {
		PowSpec p(1, lmax);
		extract_powspec(alms[0], p);
		ps.push_back(p);
	}
	else {
		for(int i = 0; i < alms.size(); i++)
			for(int j = i; j < alms.size(); j++)
				ps.push_back(extract_crosspow(alms[i],alms[j]));
		if(mode != PERMUTATIVE)
		{
			PowSpec moo = array2powspec(ps,PERMUTATIVE);
			ps = powspec2array(moo,mode,nspec);
			ps = powspec2array(array2powspec(ps,PERMUTATIVE),mode,nspec);
		}
	}
	return ps;
}

Map alm2map(const Almc & alms, int nside, Healpix_Ordering_Scheme scheme = RING)
{
	Map res(nside, RING, nside_dummy());
	alm2map(alms, res);
	if(scheme != RING) res.swap_scheme();
	return res;
}
std::vector<Map> alm2map(const std::vector<Almc> & alms, int nside, Healpix_Ordering_Scheme scheme = RING)
{
	std::vector<Map> res(alms.size(), Map(nside, RING, nside_dummy()));
	switch(alms.size())
	{
	case 1: alm2map(alms[0], res[0]); break;
	case 3: alm2map_pol(alms[0], alms[1], alms[2], res[0],
		res[1], res[2]); break;
	default: serror("Only 1 or 3 components supported!");
	}
	if(scheme != RING)
		for(int i = 0; i < res.size(); i++)
			res[i].swap_scheme();
	return res;
}

Almc map2alm(const Map & map, int lmax = -1, int mmax=-1)
{
	if(lmax<0) lmax = 3*map.Nside();
	if(mmax<0) mmax=lmax;
	Almc res(lmax,mmax);
	if(map.Scheme() == RING) map2alm(map, res, get_weights(map.Nside(),lmax), false);
	else { Map tmp = map; tmp.swap_scheme(); map2alm(tmp, res, get_weights(map.Nside(),lmax), false); }
	return res;
}
std::vector<Almc> map2alm(const std::vector<Map> & maps, int lmax, int mmax=-1)
{
	if(mmax<0) mmax = lmax;
	std::vector<Almc> res(maps.size(), Almc(lmax,mmax));
	std::vector<Map> tmap;
	if(maps[0].Scheme() != RING)
	{
		tmap = maps;
		for(int i = 0; i < maps.size(); i++)
			tmap[i].swap_scheme();
	}
	const std::vector<Map> & ms = maps[0].Scheme() == RING ? maps : tmap;
	switch(maps.size())
	{
	case 1: map2alm(ms[0], res[0], get_weights(maps[0].Nside(),lmax),false); break;
	case 3: map2alm_pol(ms[0], ms[1], ms[2], res[0],
		res[1], res[2], get_weights(maps[0].Nside(),lmax)); break;
	default: serror("Only 1 or 3 components supported!");
	}
	return res;
}

void apply_beam(std::vector<Almc> & alms, double beam)
{
	switch(alms.size())
	{
	case 1: smoothWithGauss(alms[0], beam); break;
	case 3: smoothWithGauss(alms[0], alms[1], alms[2], beam); break;
	default: serror("Only 1 or 3 components supported!");
	}
}

// Convert ascii file to powspec array, with a choice of whether to
// interpret the first column as indices or not (skipl), an lmax cap,
// and wheter to assume that the data is scaled or not.
std::vector<PowSpec> ascii2powspec(const char * filename, int lmax_max = 0, int skipl = -1, bool scale = true)
{
	if(skipl < 0) skipl = 0;
	
	std::vector<int> ls; std::vector<std::vector<double> > rows;
	std::ifstream file(filename);
	std::istream & f = filename == "-" ? std::cin : file;
	std::string line;
	int l = 0;
	while(getline(f, line))
	{
		if(line[0] == '#') continue;
		std::stringstream ss; ss << line;
		if(skipl) ls.push_back(l++);
		else { ss >> l; ls.push_back(l); }
		rows.push_back(std::vector<double>());
		double d;
		while(ss >> d) rows.back().push_back(d);
	}
	// File has been read in. Transfer to powspec.
	if(ls.empty()) serror("Error reading powspec file %s!", filename);
	if(rows.front().empty()) serror("No data columns in file %s!", filename);
	int n = rows.front().size();
	int lmax = (lmax_max && lmax_max < ls.back()) ? lmax_max : ls.back();
	std::vector<PowSpec> res;
	for(int s = 0; s < n; s++)
	{
		PowSpec tmp(1, lmax); zero_powspec(tmp);
		for(int i = 0; i < ls.size() && ls[i] <= lmax; i++)
		{
			int l = ls[i];
			double fact = scale ? l ? 2*M_PI/l/(l+1) : 0 : 1;
			tmp.tt(l) = rows[i][s]*fact;
		}
		res.push_back(tmp);
	}
	return res;
}

void powspec2ascii(const std::vector<PowSpec> & ds, const std::string & filename, int skipl = -1, bool scale = true)
{
	if(skipl < 0) skipl = 0;
	
	FILE * file = filename == "-" ? stdout : fopen(filename.c_str(),"w");
	int lmax = ds[0].Lmax();
	for(int l = 0; l <= lmax; l++)
	{
		double fact = scale ? 2*M_PI/l/(l+1) : 1;
		if(!skipl) fprintf(file, "%5d", l);
		for(int i = 0; i < ds.size(); i++)
			fprintf(file, "%17.5le", ds[i].tt(l)/fact);
		fprintf(file,"\n");
	}
	fclose(file);
}

std::vector<PowSpec> fits2powspec(const std::string & filename, int lmax_max = 0, int skipl = -1)
{
	fitshandle f; f.open(filename);
	f.goto_hdu(2);

	if(skipl < 0) skipl = !f.key_present("PDMTYPE");

	int nspec = f.ncols()-1+skipl, nrow = f.nrows();

	arr<int> idx(nrow);
	if(skipl) for(int l = 0; l < nrow; l++) idx[l] = l;
	else f.read_column(1, idx);
	// Find the true lmax by finding the max of this array
	int lmax = 0;
	for(int i = 0; i < nrow; i++) if(idx[i] > lmax) lmax = idx[i];
	lmax = (lmax_max && lmax > lmax_max) ? lmax_max : lmax;

	std::vector<PowSpec> res;
	PowSpec tmp(1, lmax);
	arr<double> a(nrow);
	int first_data_col = skipl ? 1 : 2;
	for(int i = 0; i < nspec; i++)
	{
		zero_powspec(tmp);
		f.read_column(i+first_data_col, a);
		for(int j = 0; j < nrow; j++)
			if(idx[j] <= lmax)
				tmp.tt(idx[j]) = a[j];
		res.push_back(tmp);
	}
	return res;
}

void powspec2fits(const std::vector<PowSpec> & ds, const std::string & filename, int skipl = -1)
{
	if(skipl < 0) skipl = 0;
	
	fitshandle f;
	f.create(std::string("!") + filename);

	std::vector<fitscolumn> coldescs;
	if(!skipl)
		coldescs.push_back(fitscolumn("l", "multipole order", 1, PLANCK_INT32));
	char buffer[0x100];
	for(int i = 0; i < ds.size(); i++)
	{
		sprintf(buffer, "COMPONENT%02d", i+1);
		coldescs.push_back(fitscolumn(buffer, "unknown", 1, PLANCK_FLOAT64));
	}
	f.insert_bintab(coldescs);
	int lmax = ds[0].Lmax();
	int j = 1;
	if(!skipl) {
		arr<int> ls(lmax+1);
		for(int l = 0; l < ls.size(); l++) ls[l] = l;
		f.write_column(j++, ls);
	}
	for(int i = 0; i < ds.size(); i++, j++)
		f.write_column(j, ds[i].tt());
	if(!skipl)
		f.set_key("PDMTYPE", std::string("POWERSPEC"), "Planck data model type");
}

void powspec2file(const std::vector<PowSpec> & ds, const std::string & filename, int skipl = -1, bool scale = true, int type = 0)
{
	if(type == 0 || filename.substr(filename.size()-5) != ".fits")
		powspec2ascii(ds, filename, skipl, scale);
	else
		powspec2fits (ds, filename, skipl);
}

 /****************************************************
 *  Support for working with vectors of maps, etc.   *
 ****************************************************/

typedef double (*dbinfun)(double,double);
inline double add(double a, double b) { return a+b; }
inline double sub(double a, double b) { return a-b; }
inline double mul(double a, double b) { return a*b; }
inline double div(double a, double b) { return a/b; }
Map & apply(Map & a, const Map & b, dbinfun f)
{
	for(int j = 0; j < a.Npix(); j++) a[j] = f(a[j], b[j]);
	return a;
}
Map & apply(Map & a, double b, dbinfun f)
{
	for(int j = 0; j < a.Npix(); j++) a[j] = f(a[j], b);
	return a;
}
Map & fill(Map & m, double val) { for(int j = 0; j < m.Npix(); j++) m[j] = val; return m; }
// Similar for power spectrum
PowSpec & apply(PowSpec & a, const PowSpec & b, dbinfun f)
{
	for(int j = 0; j <= a.Lmax(); j++) a.tt(j) = f(a.tt(j), b.tt(j));
	return a;
}
PowSpec & apply(PowSpec & a, double b, dbinfun f)
{
	for(int j = 0; j <= a.Lmax(); j++) a.tt(j) = f(a.tt(j), b);
	return a;
}
PowSpec & fill(PowSpec & m, double val) {
	for(int j = 0; j <= m.Lmax(); j++) m.tt(j) = val;
	return m;
}

// This did not take well to templatizing, so I'll do it manually
template <class T>
Map & operator+=(Map & a, const T & b) { return apply(a,b,add); }
template <class T>
Map & operator-=(Map & a, const T & b) { return apply(a,b,sub); }
template <class T>
Map & operator*=(Map & a, const T & b) { return apply(a,b,mul); }
template <class T>
Map & operator/=(Map & a, const T & b) { return apply(a,b,div); }
Map operator+(const Map & a, const Map & b) { Map res = a; return res += b; }
Map operator-(const Map & a, const Map & b) { Map res = a; return res -= b; }
Map operator*(const Map & a, const Map & b) { Map res = a; return res *= b; }
Map operator/(const Map & a, const Map & b) { Map res = a; return res /= b; }
Map pow(const Map & a, const Map & b) { Map res = a; apply(res,b,pow); return res; }
Map operator/(double a, const Map & b) { Map res(b); fill(res,a); return res/=b; }
Map operator/(const Map & a, double b) { Map tmp(a); fill(tmp,b); return a/tmp; }
Map operator*(double a, const Map & b) { Map res(b); fill(res,a); return res*=b; }
Map operator*(const Map & a, double b) { return b*a; }
Map operator+(double a, const Map & b) { Map res(b); fill(res,a); return res+=b; }
Map operator+(const Map & a, double b) { return b+a; }
Map operator-(double a, const Map & b) { Map res(b); fill(res,a); return res-=b; }
Map operator-(const Map & a, double b) { return b-a; }
Map pow(double a, const Map & b) { Map res(b); fill(res,a); apply(res,b,pow); return res; }
Map pow(const Map & a, double b) { Map res(a); apply(res,b,pow); return res; }

template <class T>
PowSpec & operator+=(PowSpec & a, const T & b) { return apply(a,b,add); }
template <class T>
PowSpec & operator-=(PowSpec & a, const T & b) { return apply(a,b,sub); }
template <class T>
PowSpec & operator*=(PowSpec & a, const T & b) { return apply(a,b,mul); }
template <class T>
PowSpec & operator/=(PowSpec & a, const T & b) { return apply(a,b,div); }
PowSpec operator+(const PowSpec & a, const PowSpec & b)
{ PowSpec res = a; return res += b; }
PowSpec operator-(const PowSpec & a, const PowSpec & b)
{ PowSpec res = a; return res -= b; }
PowSpec operator*(const PowSpec & a, const PowSpec & b)
{ PowSpec res = a; return res *= b; }
PowSpec operator/(const PowSpec & a, const PowSpec & b)
{ PowSpec res = a; return res /= b; }
PowSpec pow(const PowSpec & a, const PowSpec & b)
{ PowSpec res = a; apply(res,b,pow); return res; }
PowSpec operator/(double a, const PowSpec & b)
{ PowSpec res(b); fill(res,a); return res/=b; }
PowSpec operator/(const PowSpec & a, double b)
{ PowSpec tmp(a); fill(tmp,b); return a/tmp; }
PowSpec operator*(double a, const PowSpec & b)
{ PowSpec res(b); fill(res,a); return res*=b; }
PowSpec operator*(const PowSpec & a, double b) { return b*a; }
PowSpec operator+(double a, const PowSpec & b)
{ PowSpec res(b); fill(res,a); return res+=b; }
PowSpec operator+(const PowSpec & a, double b) { return b+a; }
PowSpec operator-(double a, const PowSpec & b)
{ PowSpec res(b); fill(res,a); return res-=b; }
PowSpec operator-(const PowSpec & a, double b) { return b-a; }
PowSpec pow(double a, const PowSpec & b)
{ PowSpec res(b); fill(res,a); apply(res,b,pow); return res; }
PowSpec pow(const PowSpec & a, double b)
{ PowSpec res(a); apply(res,b,pow); return res; }

// These are nice to have for general vectors. Note: The left hand side wins!
template <class T, class S> std::vector<T> & operator+=(std::vector<T> & a, const std::vector<S> & b)
{ for(int i = 0; i < a.size() && i < b.size(); i++) a[i] += b[i]; return a;}
template <class T, class S> std::vector<T> & operator-=(std::vector<T> & a, const std::vector<S> & b)
{ for(int i = 0; i < a.size() && i < b.size(); i++) a[i] -= b[i]; return a;}
template <class T, class S> std::vector<T> & operator*=(std::vector<T> & a, const std::vector<S> & b)
{ for(int i = 0; i < a.size() && i < b.size(); i++) a[i] *= b[i]; return a;}
template <class T, class S> std::vector<T> & operator/=(std::vector<T> & a, const std::vector<S> & b)
{ for(int i = 0; i < a.size() && i < b.size(); i++) a[i] /= b[i]; return a;}

template <class T, class S> std::vector<T> & operator+=(std::vector<T> & a, const S & b)
{ for(int i = 0; i < a.size(); i++) a[i] += b; return a;}
template <class T, class S> std::vector<T> & operator-=(std::vector<T> & a, const S & b)
{ for(int i = 0; i < a.size(); i++) a[i] -= b; return a;}
template <class T, class S> std::vector<T> & operator*=(std::vector<T> & a, const S & b)
{ for(int i = 0; i < a.size(); i++) a[i] *= b; return a;}
template <class T, class S> std::vector<T> & operator/=(std::vector<T> & a, const S & b)
{ for(int i = 0; i < a.size(); i++) a[i] /= b; return a;}

template <class T, class S> std::vector<T> operator+(const std::vector<T> & a, const std::vector<S>&b)
{ std::vector<T> res = a; return res += b; }
template <class T, class S> std::vector<T> operator-(const std::vector<T> & a, const std::vector<S>&b)
{ std::vector<T> res = a; return res -= b; }
template <class T, class S> std::vector<T> operator*(const std::vector<T> & a, const std::vector<S>&b)
{ std::vector<T> res = a; return res *= b; }
template <class T, class S> std::vector<T> operator/(const std::vector<T> & a, const std::vector<S>&b)
{ std::vector<T> res = a; return res /= b; }
template <class T, class S> std::vector<T> pow(const std::vector<T> & a, const std::vector<S>&b)
{ std::vector<T> res(a.size()); for(int i = 0; i < a.size() && i < b.size(); i++) res[i] = pow(a[i],b[i]);}

template <class T, class S> std::vector<T> operator+(const std::vector<T> & a, const S &b)
{ std::vector<T> res = a; return res += b; }
template <class T, class S> std::vector<T> operator-(const std::vector<T> & a, const S &b)
{ std::vector<T> res = a; return res -= b; }
template <class T, class S> std::vector<T> operator*(const std::vector<T> & a, const S &b)
{ std::vector<T> res = a; return res *= b; }
template <class T, class S> std::vector<T> operator/(const std::vector<T> & a, const S &b)
{ std::vector<T> res = a; return res /= b; }
template <class T, class S> std::vector<T> pow(const std::vector<T> & a, const S &b)
{ std::vector<T> res(a.size()); for(int i = 0; i < a.size(); i++) res[i] = pow(a[i],b);}

template <class T, class S> std::vector<T> operator+(const T & b, const std::vector<S> & a)
{ std::vector<T> res = a; return res += b; }
template <class T, class S> std::vector<T> operator-(const T & b, const std::vector<S> & a)
{ std::vector<T> res(a.size(),b); return res -= a; }
template <class T, class S> std::vector<T> operator*(const T & b, const std::vector<S> & a)
{ std::vector<T> res = a; return res *= b; }
template <class T, class S> std::vector<T> operator/(const T & b, const std::vector<S> & a)
{ std::vector<T> res(a.size(),b); return res /= a; }
template <class T, class S> std::vector<T> pow(const T & b, const std::vector<S> & a)
{ std::vector<T> res(a.size()); for(int i = 0; i < a.size(); i++) res[i] = pow(b,a[i]);}

template <class T> std::vector<T> operator-(const std::vector<T> & a)
{ std::vector<T> res(a.size()); for(int i = 0; i < a.size(); i++) res[i] = -a[i]; return res; }

// Sparsity operations
Pixset getpix(const std::vector<SMap> & a)
{
	if(a.empty()) return Pixset();
	Pixset r = a[0].get_sparsity();
	for(int i = 0; i < a.size(); i++) r &= a[i].get_sparsity();
	return r;
}

void setpix(std::vector<SMap> & a, const Pixset & pix)
{ for(int i = 0; i < a.size(); i++) a[i].set_sparsity(pix); }


#endif
