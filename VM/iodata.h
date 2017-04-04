#ifndef IODATA_H
#define IODATA_H
#include "node.h"

void read_snap(	std::string infile,
				int idx_frame,
				int &nbird,
				double &Lx,
				double &Ly,
				std::vector<float> &x,
				std::vector<float> &y,
				std::vector<float> &theta);

void ini_output(double eta,
				double epsilon,
				unsigned long long seed,
				int nStep, 
				int nCell,
				const std::string &snap_mode,
				int snap_interval);

void output(const Node *bird, int step);
#endif
