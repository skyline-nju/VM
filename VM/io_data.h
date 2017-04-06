#ifndef IO_DATA_H
#define IO_DATA_H
#include <functional>
#include "node.h"
#include "cmdline.h"
#include <ctime>

class OrderPara
{
public:
	OrderPara(
		double eta,
		double epsilon,
		double rho0,
		double Lx,
		double Ly,
		unsigned long long seed, 
		int phi_dt);

	void out(const Node *bird, int nBird, int step);

private:
	std::ofstream fout;
	int dt;
};

class iSnapshot
{
public:
	iSnapshot(const std::string file);
	~iSnapshot();
	void read();
	void read_block(int idx_frame);
	void from_file(
		int idx_frame, 
		double &_Lx, 
		double &_Ly, 
		std::vector<float> &x,
		std::vector<float> &y, 
		std::vector<float> &theta);
private:
	std::ifstream fin;
	bool is_read_block;
	double Lx;
	double Ly;
	int nBird;
	float *buff;
};

class oSnapshot
{
public:
	oSnapshot(
		int nBird,
		double _eta, 
		double _eps, 
		double _rho0, 
		double _Lx, 
		double _Ly, 
		unsigned long long int _seed,
		int out_snap_dt,
		const std::string out_mode);

	void write(const Node *bird, int nBird);
	void to_file(const Node *bird, int nBird, int step);
private:
	std::ofstream fout;
	bool is_one_file;
	int out_dt;
	double eta;
	double eps;
	double rho0;
	double Lx;
	double Ly;
	unsigned long long int seed;
};

class Output
{
public:
	Output(double rho0, double Ly, int nBird, const cmdline::parser &cmd);
	~Output();
	void out(const Node *bird, int nBird, int step);
private:
	int log_interval;
	std::ofstream fout;
	std::time_t beg_time;
	std::vector<std::function<void(const Node*, int, int)>> fout_vec;
	OrderPara *phi;
	oSnapshot *snap;
};
#endif
