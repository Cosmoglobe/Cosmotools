#ifndef TABCOMMONINCL
#define TABCOMMONINCL

#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <serror.h>

// Parse list of the form 0:1:2::3 (where empty places become -1)
// or 2:, but not plain numbers, to avoid problems with files with
// plain numbers as names
std::vector<int> parse_sel(char * sel)
{
	std::vector<int> res; char *i,*j;
	for(i = sel, j = sel; *i; i++)
	{
		if(*i == ':') {
			if(i == j) res.push_back(-1);
			else {
				char * k;
				int l = strtol(j, &k, 10);
				if(k == i) res.push_back(l);
				else return std::vector<int>();
			}
			j = i+1;
		}
	}
	if(*j && j != sel) {
		char *k; int l = strtol(j, &k, 10);
		if(*k) return std::vector<int>(); // Trailing characters
		res.push_back(l);
	}
	return res;
}

std::vector<std::string> split(const std::string & s, const std::string & sep)
{
	std::vector<std::string> res;
	int i = 0, j = 0;
	for(int i = 0, j = 0; i < s.size() && j != std::string::npos; i = j+1)
	{
		j = s.find_first_of(sep, i);
		res.push_back(s.substr(i, j == std::string::npos ?
			std::string::npos : j-i));
	}
	return res;
}

std::vector<std::vector<int> > parse_sel_cross(char * sel)
{
	std::vector<std::string> a = split(sel,":");
	std::vector<std::vector<int> > res;
	if(a.size() <= 1) return res;
	for(int i = 0; i < a.size(); i++)
	{
		if(a[i].empty()) continue;
		std::vector<std::string> tmp = split(a[i],"x-*");
		std::vector<int> work;
		for(int j = 0; j < tmp.size(); j++)
		{
			char *b;
			long l = strtol(tmp[j].c_str(), &b, 10);
			if(*b) return std::vector<std::vector<int> >();
			work.push_back(atoi(tmp[j].c_str()));
		}
		if(!work.empty()) res.push_back(work);
	}
	return res;
}

int count0(const std::string & s)
{
	std::stringstream ss; ss << s;
	int c = 0; double d;
	while(ss >> d) c++;
	return c;
}

int count(const std::string & s)
{
	const char * b = s.c_str();
	int c = 0, r = 1, n; double d;
	for(int i = 0; r && i < s.size(); c += r != 0, i+= n) 
		r = sscanf(b+i,"%lg%n", &d, &n);
	return c;
}

// Read one row from selected columns of file. Return empty if done
template <class T>
std::vector<T> read_colrowT(std::istream & ifile, std::vector<int> & cols)
{
	std::vector<T> table(cols.size());
	std::string line;
	if(std::getline(ifile,line))
	{
		std::stringstream ss; ss << line;
		if(cols.empty()) {
			// process everything
			int n = count(line);
			for(int i = 0; i < n; i++) cols.push_back(i+1);
			table.resize(cols.size());
		}
		std::vector<T> nums;
		T num;
		while(ss >> num) nums.push_back(num);
		for(int i = 0; i < cols.size(); i++)
			if(cols[i] > 0) table[i] = nums[cols[i]-1];
			else table[i] = T();
	}
	else return std::vector<T>();
	return table;
}
std::vector<double> read_colrow(std::istream & ifile, std::vector<int> & cols)
{ return read_colrowT<double>(ifile, cols); }

std::vector<std::vector<double> > read_cols(std::istream & ifile, std::vector<int> & cols)
{
	std::vector<std::vector<double> > table;
	for(std::vector<double> row = read_colrow(ifile, cols); !row.empty(); row = read_colrow(ifile,cols))
	{
		if(row.size() != table.size()) table.resize(row.size()); // can only happen at the beginning
		for(int i = 0; i < row.size(); i++) table[i].push_back(row[i]);
	}
	return table;
}

std::vector<std::vector<double> > read_cols(char * filename, std::vector<int> & cols)
{
	if(filename && strcmp(filename,"-")) { std::ifstream in(filename); return read_cols(in, cols); }
	else return read_cols(std::cin, cols);
}

#endif
