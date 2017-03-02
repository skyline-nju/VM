#include "output.h"
using namespace std;

ofstream fout_phi;
ofstream fout_log;
ofstream fout_snap;
char para[100];
clock_t t_begin;

void output_ini(double eta, double epsilon, int seed, int Nstep)
{
	t_begin = clock();
	time_t nowtime = time(NULL);
	tm t;
#ifdef _MSC_VER
	localtime_s(&t, &nowtime);
#else
	localtime_r(&nowtime, &t);
#endif

	snprintf(para, 100, "%g_%g_%g_%d_%d_%d", eta, epsilon, Node::rho_0, int(Node::Lx), int(Node::Ly), seed);
	char buff_phi[100];
	char buff_log[100];
#ifdef _MSC_VER
	mkdir("..\\phi");
	mkdir("..\\log");
	mkdir("..\\snap");
	snprintf(buff_phi, 100, "..\\phi\\p_%s.dat", para);
	snprintf(buff_log, 100, "..\\log\\l_%s.dat", para);
#else
	mkdir("../phi");
	mkdir("../log");
	mkdir("../snap");
	snprintf(buff_phi, 100, "../phi/p_%s.dat", para);
	snprintf(buff_log, 100, "../log/l_%s.dat", para);
#endif

	fout_phi.open(buff_phi);
	fout_log.open(buff_log);

	fout_log << t.tm_year + 1900 << "/" << t.tm_mon + 1 << "/" << t.tm_mday << "/" << t.tm_hour << ":" << t.tm_min << endl;
	fout_log << "eta = " << eta << endl;
	fout_log << "epsilon = " << epsilon << endl;
	fout_log << "rho_0 = " << Node::rho_0 << endl;
	fout_log << "Lx = " << Node::Lx << endl;
	fout_log << "Ly = " << Node::Ly << endl;
	fout_log << "N = " << Node::N << endl;
	fout_log << "nCell = " << Grid::mm << endl;
	fout_log << "seed = " << seed << endl;
	fout_log << "total step: " << Nstep << endl;
	fout_log << endl;
	fout_log << "step\ttime" << endl;
}

void output_phi(const Node * bird, int step)
{
	double svx = 0;
	double svy = 0;
	for (int i = 0; i < Node::N; i++)
	{
		svx += bird[i].vx;
		svy += bird[i].vy;
	}
	double phi = sqrt(svx*svx + svy*svy) / Node::N;
	double theta = atan2(svy, svx);
	fout_phi << fixed << std::setw(16) << setprecision(10) << step << "\t" << phi << "\t" << theta << endl;
}

void output_log(int step)
{
	clock_t t_now = clock();
	double dt = double(t_now - t_begin) / CLOCKS_PER_SEC / 3600;
	fout_log << step << "\t" << dt << endl;
#ifdef _MSC_VER
	cout << step << "\t" << dt << endl;
#endif
}

void output_snap(const Node * bird, int step)
{
	char file[100];
#ifdef _MSC_VER
	snprintf(file, 100, "..\\snap\\s_%s_%08d.bin", para, step);
#else
	snprintf(file, 100, "../snap/s_%s_%08d.bin", para, step);
#endif
	fout_snap.open(file, ios::binary);
	float *buff = new float[3 * Node::N];
	for (int j = 0; j < Node::N; j++)
	{
		buff[3 * j] = bird[j].x;
		buff[3 * j + 1] = bird[j].y;
		buff[3 * j + 2] = atan2(bird[j].vy, bird[j].vx);
	}
	fout_snap.write((char*)&buff[0], sizeof(float) * Node::N * 3);
	fout_snap.close();
	delete[] buff;
	buff = NULL;
}

void output(const Node * bird, int step)
{
	if (step % 100 == 0)
	{
		output_phi(bird, step);
		if (step % 100000 == 0)
		{
			output_log(step);
			output_snap(bird, step);
		}
	}
}

void mkdir(char * folder)
{
#ifdef _MSC_VER
	if (_access(folder, 0) != 0)
#else
	if (access(folder, 0) != 0)
#endif
	{
		char command[100];
		snprintf(command, 100, "mkdir %s", folder);
		if (system(command))
			cout << "create folder: " << folder << " successfully" << endl;
	}
	else
		cout << "folder: " << folder << " already exists" << endl;
}
