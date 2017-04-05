#ifndef IO_DATA_H
#define IO_DATA_H
#include <functional>
#include "node.h"

class Log
{
public:
	Log(double eta, 
		double epsilon, 
		double rho0, 
		double Lx, 
		double Ly, 
		unsigned long long seed,
		int log_dt);
	void out(const Node *bird, int nBird, int step);
	template <typename T>
	Log &operator << (const T &item);

private:
	std::ofstream fout;
	int dt;
};

template<typename T>
inline Log & Log::operator << (const T & item)
{
	cout << item;
	fout << item;
	Log &tmp = *this;
	return tmp;
}

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
	void from_file(int idx_frame, double &_Lx, double &_Ly, 
		std::vector<float> &x, std::vector<float> &y, std::vector<float> &theta);
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
	Output(
		double eta,
		double eps,
		double rho0,
		double Lx,
		double Ly,
		int nBird,
		unsigned long long int seed,
		int nStep,
		int log_dt,
		int phi_dt,
		int snap_dt,
		const std::string snap_mode);
	~Output();
	void write(const Node *bird, int nBird, int step);
private:
	std::vector<std::function<void(const Node*, int, int)>> out_vec;
	Log *log;
	OrderPara *phi;
	oSnapshot *snap;
};
#endif
