#ifndef COMN_H
#define COMN_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

const double PI = 3.14159265358979;

void mkdir(const char *folder);

template <class T>
void str_to_num(const char *str, T &num)
{
	std::stringstream ss;
	ss << str;
	ss >> num;
}

template <class T>
void num_to_str(const T &num, char * str)
{
	std::stringstream ss;
	ss << num;
	ss >> str;
}

#endif
