#include "io_data.h"

using namespace std;

#ifdef _MSC_VER
string delimiter("\\");
#else
string delimiter("/");
#endif

OrderPara::OrderPara(
	double eta,
	double eps, 
	double rho0, 
	double Lx, 
	double Ly, 
	unsigned long long seed,
	int phi_dt)
{
	mkdir("phi");
	char buff[100];
	snprintf(buff, 100, 
		"phi%sp_%g_%g_%g_%g_%g_%llu.dat", 
		delimiter.c_str(), eta, eps, rho0, Lx, Ly, seed);
	fout.open(buff);
	dt = phi_dt;
}

void OrderPara::out(const Node * bird, int nBird, int step)
{
	if (step % dt == 0)
	{
		double svx = 0;
		double svy = 0;
		for (int i = 0; i < nBird; i++)
		{
			svx += bird[i].vx;
			svy += bird[i].vy;
		}
		double phi = sqrt(svx * svx + svy * svy) / Node::N;
		double theta = atan2(svy, svx);
		fout << fixed << std::setw(16) << setprecision(10)
			<< step << "\t" << phi << "\t" << theta << endl;
	}

}

oSnapshot::oSnapshot(
	int nBird,
	double _eta,
	double _eps,
	double _rho0,
	double _Lx,
	double _Ly,
	unsigned long long int _seed,
	int out_snap_dt,
	const std::string out_mod)
{
	eta = _eta;
	eps = _eps;
	rho0 = _rho0;
	Lx = _Lx;
	Ly = _Ly;
	seed = _seed;
	out_dt = out_snap_dt;
	if (out_mod == "one")
	{
		mkdir("snap_one");
		char filename[100];
		snprintf(filename, 100, "snap_one%sso_%g_%g_%g_%g_%d_%d_%llu.bin",
			delimiter.c_str(), eta, eps, Lx, Ly, nBird, out_dt, seed);
		is_one_file = true;
		fout.open(filename, ios::binary);
	}
	else if (out_mod == "mult")
	{
		mkdir("snap");
		is_one_file = false;
	}
}

void oSnapshot::write(const Node * bird, int nBird)
{
	float *buff = new float[3 * nBird];
	for (int j = 0; j < nBird; j++)
	{
		buff[3 * j] = bird[j].x;
		buff[3 * j + 1] = bird[j].y;
		buff[3 * j + 2] = atan2(bird[j].vy, bird[j].vx);
	}
	fout.write((char *)buff, sizeof(float) * nBird * 3);
	delete[] buff;
	buff = nullptr;
}

void oSnapshot::to_file(const Node *bird, int nBird, int step)
{
	if (step % out_dt == 0)
	{
		if (is_one_file)
			write(bird, nBird);
		else
		{
			char filename[100];
			snprintf(filename, 100, "snap%ss_%g_%g_%g_%g_%g_%llu_%08d.bin",
				delimiter.c_str(), eta, eps, rho0, Lx, Ly, seed, step);
			fout.open(filename, ios::binary);
			write(bird, nBird);
			fout.close();
		}
	}
}

iSnapshot::iSnapshot(const string infile)
{
	vector<string> str_list = split(infile, "_");
	string folder;
	if (str_list[0] == "s")
	{
		folder = "snap";
		is_read_block = false;
		str_to_num(str_list[4], Lx);
		str_to_num(str_list[5], Ly);
		buff = nullptr;
	}
	else
	{
		folder = "snap_one";
		is_read_block = true;
		str_to_num(str_list[5], nBird);
		str_to_num(str_list[3], Lx);
		str_to_num(str_list[4], Ly);
		buff = new float[nBird * 3];
	}
	string path(folder + delimiter + infile);
	fin.open(path, ios::binary | ios::ate);
	if (!fin.is_open())
	{
		cout << "Error, failed to open " << infile << endl;
		exit(1);
	}
	else
		cout << "Open snapshots: " << infile << endl;
}

iSnapshot::~iSnapshot()
{
	fin.close();
	delete[]buff;
	buff = nullptr;
}

void iSnapshot::read()
{
	streamsize count = fin.tellg();
	nBird = count / sizeof(float) / 3;
	buff = new float[nBird * 3];
	fin.seekg(0, ios::beg);
	fin.read((char *)buff, count);
}

void iSnapshot::read_block(int idx_frame)
{
	streamsize count = fin.tellg();
	streamsize frame_size = sizeof(float) * nBird * 3;
	if (idx_frame < 0)
		fin.seekg(frame_size * idx_frame, ios::end);
	else
		fin.seekg(frame_size * idx_frame, ios::beg);
	fin.read((char *)buff, frame_size);
}

void iSnapshot::from_file(
	int idx_frame,
	double & _Lx, 
	double & _Ly, 
	vector<float> &x, 
	vector<float> &y, 
	vector<float> &theta)
{
	if (is_read_block)
		read_block(idx_frame);
	else
		read();
	_Lx = Lx;
	_Ly = Ly;

	x.reserve(nBird);
	y.reserve(nBird);
	theta.reserve(nBird);
	for (int i = 0; i < nBird; i++)
	{
		x.push_back(buff[3 * i]);
		y.push_back(buff[3 * i + 1]);
		theta.push_back(buff[3 * i + 2]);
	}
}

Output::Output(double rho0, double Ly, int nBird, const cmdline::parser &cmd)
{
	using std::placeholders::_1;
	using std::placeholders::_2;
	using std::placeholders::_3;
	
	double eta = cmd.get<double>("eta");
	double eps = cmd.get<double>("eps");
	double Lx = cmd.get<double>("Lx");
	unsigned long long seed = cmd.get<unsigned long long>("seed");
	
	time(&beg_time);
	mkdir("log");
	char logfile[100];
	snprintf(logfile, 100,
		"log%sl_%g_%g_%g_%g_%g_%llu.dat",
		delimiter.c_str(), eta, eps, rho0, Lx, Ly, seed);
	fout.open(logfile);
	log_interval = cmd.get<int>("log_dt");

	phi = new OrderPara(eta, eps, rho0, Lx, Ly, seed, cmd.get<int>("phi_dt"));
	fout_vec.push_back(bind(&OrderPara::out, phi, _1, _2, _3));

	string snap_mode = cmd.get<string>("snap_mode");
	if (snap_mode != "none")
	{
		int snap_dt = cmd.get<int>("snap_dt");
		snap = new oSnapshot(
			nBird, eta, eps, rho0, Lx, Ly, seed, snap_dt, snap_mode);
		fout_vec.push_back(bind(&oSnapshot::to_file, snap, _1, _2, _3));
	}

	fout << "Started at " << asctime(localtime(&beg_time));
	fout << "Parameters:\n";
	fout << "Particle number = " << nBird << endl;
	fout << "density = " << rho0 << endl;
	fout << "eta = " << eta << endl;
	fout << "epsilon = " << eps << endl;
	fout << "Lx = " << Lx << endl;
	fout << "Ly = " << Ly << endl;
	fout << "seed = " << seed << endl;
	fout << "total time steps = " << cmd.get<int>("nstep") << endl;
	if (cmd.exist("file"))
		fout << "ini mode: " << cmd.get<string>("file") << endl;
	else
		fout << "ini mode: " << cmd.get<string>("ini_mode") << endl;
	fout << "log dt = " << log_interval << endl;
	fout << "phi dt = " << cmd.get<int>("phi_dt") << endl;
	fout << "snap mode = " << snap_mode << endl;
	fout << "snap dt = " << cmd.get<int>("snap_dt") << endl;
	fout << endl;
	fout << "-------- Run --------\n";
}
Output::~Output()
{
	delete phi;
	delete snap;
	time_t end_time = time(nullptr);
	fout << "Finished at " << asctime(localtime(&end_time));
	fout.close();
}

void Output::out(const Node * bird, int nBird, int step)
{
	for (auto f: fout_vec)
	{
		f(bird, nBird, step);
	}
	if (step % log_interval == 0)
	{
		double dt = difftime(time(nullptr), beg_time);
		int hour = int(dt / 3600);
		int min = int((dt - hour * 3600) / 60);
		double sec = dt - hour * 3600 - min * 60;
		fout << step << "\t" << hour << ":" << min << ":" << sec << endl;
	}
}
