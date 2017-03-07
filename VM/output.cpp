#include "output.h"
#define COUT_PHI
using namespace std;

ofstream fout_phi;
ofstream fout_log;
ofstream fout_snap;
char para[100];
clock_t t_begin;
int tot_step;
double sum_phi = 0;
int count_phi = 0;

void output_ini(double eta, double epsilon, unsigned long long seed, int nStep)
{
    t_begin = clock();
    time_t nowtime = time(NULL);
    tm t;
    tot_step = nStep;
#ifdef _MSC_VER
    localtime_s(&t, &nowtime);
#else
    localtime_r(&nowtime, &t);
#endif

    char seed_str[30];
    num_to_str(seed, seed_str);
    snprintf(para, 100, "%g_%g_%g_%d_%d_%s", eta, epsilon, Node::rho_0, int(Node::Lx), int(Node::Ly), seed_str);
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
    fout_log << "total step: " << nStep << endl;
    fout_log << endl;
    fout_log << "step\th:m:s" << endl;
}

void output_phi(const Node *bird, int step)
{
    double svx = 0;
    double svy = 0;
    for (int i = 0; i < Node::N; i++)
    {
	svx += bird[i].vx;
	svy += bird[i].vy;
    }
    double phi = sqrt(svx * svx + svy * svy) / Node::N;
    double theta = atan2(svy, svx);
    fout_phi << fixed << std::setw(16) << setprecision(10) << step << "\t" << phi << "\t" << theta << endl;

    if (step > 50000)
    {
	sum_phi += phi;
	count_phi += 100;
    }
}

void output_log(int step)
{
    clock_t t_now = clock();
    int dt = (t_now - t_begin) / CLOCKS_PER_SEC;
    int hour = dt / 3600;
    int min = (dt - hour * 3600) / 60;
    int sec = dt - hour * 3600 - min * 60;
    fout_log << step << "\t" << hour << ":" << min << ":" << sec << endl;
    if (step == tot_step)
    {
	fout_log << "phi = " << sum_phi / tot_step << endl;
	double speed = double(step) / dt * 3600;
	fout_log << "simulation speed: " << int(speed) << " step per hour" << endl;
    }
}

void output_snap(const Node *bird, int step)
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
    fout_snap.write((char *)&buff[0], sizeof(float) * Node::N * 3);
    fout_snap.close();
    delete[] buff;
    buff = NULL;
}

void output(const Node *bird, int step)
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

void mkdir(const char *folder)
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
