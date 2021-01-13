#include <seph.h>
#include <serror.h>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>
using namespace std;
using namespace skn;

int main(int argn, char ** argv)
{
	try {
		vector<char*> args;
		char * database = "/projects/quiet/external_data/unix.405",
			* patches = "/projects/quiet/auxilliary/patches_galactic.txt";
		bool ces_detect_format = false;
		for(char ** i = argv+1; *i; i++)
			if(!strcmp("-db",*i)) database = *++i;
			else if(!strcmp("-c",*i)) ces_detect_format ^= true;
			else if(!strcmp("-patches",*i)) patches = *++i;
			else args.push_back(*i);

		Object test = { "jupiter", Vector<double,2>(0,0), 0.1, false, 1 };
		vector<Object> objdb = read_objs(patches);
		Ephfile eph(database);

		if(!ces_detect_format)
		{
			double mjd, lon, lat;
			while(cin >> mjd >> lon >> lat)
				cout << identify(eph, objdb, mjd, lon*M_PI/180, lat*M_PI/180) << endl;
		}
		else
		{
			string line;
			while(getline(cin,line))
			{
				double mjd1, mjd2, az, el, dk, lon, lat, dev;
				stringstream ss(line);
				double d;
				if(ss >> mjd1 >> mjd2 >> az >> el >> dk >> lon >> lat >> dev)
				{
					string p = identify(eph, objdb, (mjd1+mjd2)/2, lon*M_PI/180, lat*M_PI/180,&d);
					printf("%s %-10s %8.4lf\n", line.c_str(), p.c_str(), d*180/M_PI);
				}
			}
		}
	} catch(Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); }
}
