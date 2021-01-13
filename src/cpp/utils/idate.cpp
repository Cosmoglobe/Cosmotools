#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>

using namespace std;

int main(int argn, char **argv)
{
	if(argn < 4) {
		fprintf(stderr, "Usage: idate date input_format output_format\n");
		return 1;
	}
	char * input = argv[1];
	char * ifmt = argv[2];
	char * ofmt = argv[3];
	char buffer[0x1000];
	struct tm tim;
	for(int i = 0; i < sizeof(tim); i++) ((char*)&tim)[i] = 0;
	strptime(input, ifmt, &tim);
	strftime(buffer, 0x1000, ofmt, &tim);
	printf("%s\n", buffer);
	return 0;
}
