#include <seph.h>
#include <iostream>
#include <serror.h>
#include <cstdio>
#include <cstring>

using namespace std;
using namespace skn;

typedef Vector<double,3> Vec;

string to_minsec(double r, double tot, bool positive)
{
	if(positive && r < 0) r += 2*M_PI;
	r *= tot/(2*M_PI);
	int sign = r < 0 ? -1 : 1;
	r *= sign;
	int h = int(r);
	double tmp = (r-h)*60;
	int m = int(tmp);
	double s = (tmp-m)*60;
	h *= sign;
	char buf[0x100];
	sprintf(buf,"%4d:%02d:%05.2lf", h, m, s);
	return buf;
}

void help()
{
	fprintf(stderr,
"eph, a simple command line ephemeris utility, built on ephcom.\n"
"Syntax: eph [options] object mjd [mjd2 [npoint]].\n"
"Options:\n"
"  -rec: Rectangular coordinates. (Default is polar (degrees)).\n"
"  -r:   Polar coordinates (radians).\n"
"  -h:   Polar coordinates (degrees, minutes, seconds).\n"
"  -t:   As -h, but the first angle is measures in hours.\n"
"  -cel: Celestial coordinates. (Default is galactic).\n"
"  -gal: Galactic coordinates.\n"
"  -s:   Sampling rate. Interval between points when in table mode.\n"
"  -npoint: Number of points to use in table mode. Conflicts with -s.\n"
"  -db:  Location of the 405 ephemeris file from JPL.\n"
"Objects: sun, venus, earth, moon, mars, jupiter, saturn, uranus,\n"
"         neptune and pluto. Names can be shortened, but must be\n"
"         lower case.\n"
"Examples:\n"
"  * Get the galactic coordinates of Venus in degrees at mjd 54322.32\n"
"     eph venus 54322.32\n"
"    output:   54322.3200000  52.525556  42.911547   0.29771\n"
"  * Get a table of the rectangular celestial coordinates of jupiter\n"
"    between mjd 54000 and 56000 with intervals of 10 days:\n"
"     eph jup -s 10 -cel -rec 54000 56000\n"
"    output:\n"
"     54000.0000000  -4.1196936  -4.0513872  -1.6605850\n"
"     54010.0000000  -4.0448767  -4.2440856  -1.7454376\n"
"       .....\n"
"     55990.0000000   4.3725822   3.1046392   1.2461877\n"
"     56000.0000000   4.3529667   3.3098037   1.3362565\n");
}

enum Outformat { rectangular, radians, degrees, degminsec, traditional };

void do_mjd(double mjd, Ephfile & eph, int obj, bool galactic, Outformat ofmt)
{
	double jd = mjd+mdiff;
	Vec rec = get_pos_rel(obj, jd, eph);
	if(galactic) rec = equ2gal * rec;
	Vec pos = rec2pol(rec);
	switch(ofmt)
	{
	case rectangular:
		printf("%15.7lf %11.7lf %11.7lf %11.7lf\n", mjd, rec[0], rec[1], rec[2]);
		break;
	case radians:
		printf("%15.7lf %11.7lf %11.7lf %9.5lf\n", mjd, pos[0], pos[1], pos[2]);
		break;
	case degrees:
		printf("%15.7lf %10.6lf %10.6lf %9.5lf\n", mjd, pos[0]*180/M_PI, pos[1]*180/M_PI, pos[2]);
		break;
	case degminsec:
		printf("%15.7lf %s %s %9.5lf\n", mjd, to_minsec(pos[0],360,true).c_str(), to_minsec(pos[1],360,false).c_str(), pos[2]);
		break;
	case traditional:
		printf("%15.7lf %s %s %9.5lf\n", mjd, to_minsec(pos[0],24,true).c_str(), to_minsec(pos[1],360,false).c_str(), pos[2]);
		break;
	}
}

int main(int argn, char ** argv)
{
	try {
		vector<char*> args;
		char * database = "/projects/quiet/external_data/unix.405";
		int npoint = 0;
		double srate = 0;
		bool galactic = true;
		Outformat ofmt = degrees;
		for(char ** i = argv+1; *i; i++)
			if(!strcmp("-db",*i)) database = *++i;
			else if(!strcmp("-c",*i) || match(*i,"-rectangular",3))
				ofmt = rectangular;
			else if(!strcmp("-r",*i)) ofmt = radians;
			else if(!strcmp("-d",*i)) ofmt = degrees;
			else if(!strcmp("-h",*i)) ofmt = degminsec;
			else if(!strcmp("-t",*i)) ofmt = traditional;
			else if(!strcmp("-n",*i)) npoint = atoi(*++i);
			else if(!strcmp("-s",*i)) srate = atof(*++i);
			else if(!strcmp("-gal",*i)) galactic = true;
			else if(!strcmp("-cel",*i)) galactic = false;
			else if(!strcmp("-h",*i)) { help(); return 0; }
			else args.push_back(*i);

		if(args.size() < 1) { help(); return 1; }
		string object = args[0];
		int obj = lookup_object(object);
		Ephfile eph(database);

		double mjd;
		if(args.size() == 1)
			while(cin >> mjd) do_mjd(mjd, eph, obj, galactic, ofmt);
		else
		{
			double mjd1 = atof(args[1]); double mjd2 = mjd1;
			if(args.size() >= 3) mjd2 = atof(args[2]);
			if(args.size() >= 4) npoint = atoi(args[3]);
			if(npoint == 0 && srate != 0) {
				npoint = int((mjd2-mjd1)/srate)+1;
				mjd2 = mjd1 + srate*(npoint-1);
			}
			else if(npoint == 0) npoint = 1;

			eph.assert_in_range(mjd1+mdiff);
			eph.assert_in_range(mjd2+mdiff);

			for(int i = 0; i < npoint; i++)
			{
				double mjd = npoint == 1 ? mjd1 :
					mjd1+(mjd2-mjd1)*i/(npoint-1);
				do_mjd(mjd, eph, obj, galactic, ofmt);
			}
		}
	} catch(Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); }
}
