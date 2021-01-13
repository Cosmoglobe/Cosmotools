#include <cstdlib>
#include <cstdio>
#include <glob.h>
#include <vector>
#include <signal.h>

using namespace std;

int mod(int a, int b) { int c = a%b; return c < 0 ? c+b : c; }

int main(int argc, char ** argv)
{
	vector<char*> args;
	for(char **i = argv+1; *i; i++)
		if(false);
		else args.push_back(*i);

	if(args.empty()) {
		fprintf(stderr,"Syntax: mpifor [options] command [args]\n");
		return 1;
	}

	int nprocs = atoi(getenv("SLURM_NPROCS")),
		procid = atoi(getenv("SLURM_PROCID"));

	// And parse and run
	char* command = args[0];
	for(int i = 1; i < args.size(); i++)
	{
		int k = i-1;
		glob_t p;
		p.gl_offs = 1;
		glob(args[i], GLOB_BRACE | GLOB_TILDE | GLOB_NOMAGIC, 0, &p);
		for(int j = mod(procid-k,nprocs); j < p.gl_pathc; j += nprocs)
		{
			setenv("i", p.gl_pathv[j], 1);
			int status = system(command);
			if (WIFSIGNALED(status) &&
				(WTERMSIG(status) == SIGINT || WTERMSIG(status) == SIGQUIT))
				return 1;
		}
		globfree(&p);
	}
	return 0;
}
