#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

using namespace std;

void eatline(FILE * f)
{
	int ret;
	while((ret = fgetc(f)) != '\n' && ret != EOF);
}

bool read_header(FILE * file, int & n)
{
	const int bufsize = 0x1000;
	char line[bufsize];
	fscanf(file, " ");
	if(!fgets(line, bufsize, file)) return false;
	if(sscanf(line, " #Curve %*d of %*d, %d points\n", &n) == 1) return 1;
	else if(sscanf(line, " #Curve %*d, %d points\n", &n) == 1) return 1;
	printf("%s is impossible to understand!\n", line);
	return false;
}

int main(int argn, char ** argv)
{
	char * infile = 0, * outfile = 0, * fmt = "%25.17le";
	vector<char*> args;
	for(char ** i = argv+1; *i; i++)
		if(!strcmp(*i,"-fmt")) fmt = *++i;
		else args.push_back(*i);

	if(args.size() > 0) infile = args[0];
	if(args.size() > 1) outfile = args[1];

	FILE * in = infile && strcmp(infile,"-") != 0 ? fopen(infile,"r") : stdin;

	vector<vector<double> > table;

	// Read in table
	int num, count = 0;
	while(read_header(in, num))
	{
		if(count == 0) table.push_back(vector<double>(num));
		table.push_back(vector<double>(num));
		if(table[0].size() != table.back().size()) {
			fprintf(stderr, "Column %d's length (%d) and %d (%d) do not match!",
				table.size()-1, table.back().size(), 0, table[0].size());
			return 1;
		}
		fscanf(in, "#x y type\n");
		double x, y; char mode;
		for(int i = 0; i < num; i++)
		{
			fscanf(in, "%lf %lf %c\n", &x, &y, &mode);
			if(count == 0) table[0][i] = x;
			else if(table[0][i] != x) {
				fprintf(stderr, "Inconsistent x values!");
				return 1;
			}
			table.back()[i] = y;
		}
		count++;
	}
	// And output it
	if(table.empty()) return 1;
	FILE * out = outfile && strcmp(outfile,"-") != 0 ? fopen(outfile,"w") : stdout;
	for(int i = 0; i < table[0].size(); i++)
	{
		for(int j = 0; j < table.size(); j++)
			fprintf(out, fmt, table[j][i]);
		fprintf(out, "\n");
	}
}
