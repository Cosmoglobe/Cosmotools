#include "tabcommon.h"
#include "serror.h"
using namespace std;
using namespace skn;

int main(int argn, char ** argv)
{
	try {
		vector<int> tcol;
		vector<vector<int> > cols;
		vector<char *> filenames;
		for(char ** i = argv+1; *i; i++)
			if(!(tcol = parse_sel(*i)).empty()) {
				if(cols.empty()) serror("Syntax: file [cols] file [cols] ...");
				cols.back() = tcol;
			}
			else {
				filenames.push_back(*i);
				cols.push_back(tcol);
			}
		// Split into input/output parts
		vector<char*> infiles;
		char * outfile;
		vector<vector<int> > incols;
		vector<int> outcols;
		for(int i = 0; i < filenames.size(); i++)
		{
			infiles.push_back(filenames[i]);
			incols.push_back(cols[i]);
		}
		outfile = "-";

		vector<istream*> bar(infiles.size());
		for(int i = 0; i < infiles.size(); i++)
			if(infiles[i] && strcmp(infiles[i],"-"))
				bar[i] = new ifstream(infiles[i]);
			else bar[i] = &cin;
		FILE * out = outfile && strcmp(outfile,"-") ? fopen(outfile,"w") : stdout;

		while(true)
		{
			vector<string> tot_row;
			for(int i = 0; i < infiles.size(); i++)
			{
				vector<string> row = read_colrowT<string>(*bar[i], cols[i]);
				if(row.empty()) goto done;
				for(int j = 0; j < row.size(); j++) tot_row.push_back(row[j]);
			}
			if(outcols.empty()) for(int i = 0; i < tot_row.size(); i++)
				fprintf(out, "%s ", tot_row[i].c_str());
			else for(int i = 0; i < outcols.size(); i++)
				fprintf(out, "%s ", tot_row[outcols[i]-1].c_str());
			fprintf(out, "\n");
		}
		done: for(int i = 0; i < bar.size(); i++)
			if(bar[i] != &cin) delete bar[i];
	} catch(Error & e) { fprintf(stderr, e.msg.c_str()); }
}
