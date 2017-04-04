#include "iodata.h"

using namespace std;

#ifdef _MSC_VER
string dlm("\\");
#else
string dlm("/");
#endif

ofstream fout_phi;
ofstream fout_log;
ofstream fout_snap;

char para[100];
int dt_phi = 100;
int dt_log = 100000;
int dt_snap = 2000;
bool flag_one_snap_file = true;
int tot_step;

/* read a snapshot */
void read_snap(	string infile,
				int idx_frame,
				int &nbird, 
				double &Lx,
				double &Ly,
				vector<float> &x, 
				vector<float> &y, 
				vector<float> &theta)
{
	vector<string> str_list = split(infile, "_");
	bool is_single_frame = str_list[0] == "s" ? true : false;
	string folder = is_single_frame ? "snap" : "snap_one";
	string path(folder + dlm + infile);
	ifstream fin(path, ios::binary | ios::ate);
	if (!fin.is_open())
	{
		cout << "Error, failed to open " << infile << endl;
		exit(1);
	}
	streamsize count = fin.tellg();
	int frame_size;
	if (is_single_frame)
	{
		nbird = count / sizeof(float) / 3;
		str_to_num(str_list[4], Lx);
		str_to_num(str_list[5], Ly);
		frame_size = count;
	}
	else
	{
		str_to_num(str_list[5], nbird);
		str_to_num(str_list[3], Lx);
		str_to_num(str_list[4], Ly);
		cout << "Lx0 = " << Lx << endl;
		cout << "Ly0 = " << Ly << endl;
		frame_size = sizeof(float) * nbird * 3;
		if (count % nbird != 0)
		{
			cout << "error when read " << infile << endl;
			exit(1);
		}
	}
	cout << "size of file " << count << endl;
	cout << "num of birds = " << nbird << endl;
	cout << "frame size = " << frame_size << endl;
	cout << "idx_frame = " << idx_frame << endl;
	if (idx_frame < 0)
		fin.seekg(-frame_size * idx_frame, ios::end);
	else
		fin.seekg(frame_size * idx_frame, ios::beg);
	cout << (count - fin.tellg()) / 3 / sizeof(float) << endl;
	float *buff = new float[nbird * 3];
	fin.read((char*)&buff, frame_size);
	fin.close();
	
	x.reserve(nbird);
	y.reserve(nbird);
	theta.reserve(nbird);
	for (int i = 0; i < nbird; i++)
	{
		x.push_back(buff[3 * i]);
		y.push_back(buff[3 * i + 1]);
		theta.push_back(buff[3 * i + 2]);
	}

	delete[] buff;
	buff = nullptr;
}

/* initialize output */
void ini_output(double eta, 
				double epsilon, 
				unsigned long long seed, 
				int nStep, 
				int nCell, 
				const string &snap_mode, 
				int snap_interval)
{
	tot_step = nStep;
	dt_snap = snap_interval;
	snprintf(para, 100, "%g_%g_%g_%d_%d_%llu",
		eta, epsilon, Node::rho_0, int(Node::Lx), int(Node::Ly), seed);
	mkdir("phi");
	mkdir("log");
	char buff[100];
	snprintf(buff, 100, "phi%sp_%s.dat", dlm.c_str(), para);
	fout_phi.open(buff);
	snprintf(buff, 100, "log%sl_%s.dat", dlm.c_str(), para);
	fout_log.open(buff);

	if (snap_mode == "one")
	{
		mkdir("snap_one");
		flag_one_snap_file = true;
		snprintf(buff, 100, "snap_one%sso_%g_%g_%d_%d_%d_%d_%llu.bin",
			dlm.c_str(), eta, epsilon, int(Node::Lx), int(Node::Ly), Node::N, dt_snap, seed);
		fout_snap.open(buff, ios::binary);
	}
	else if (snap_mode == "mult")
	{
		mkdir("snap");
		flag_one_snap_file = false;
	}
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
	fout_phi << fixed << std::setw(16) << setprecision(10)
		<< step << "\t" << phi << "\t" << theta << endl;
}

void output_log(int step)
{
	fout_log << step << endl;
}

/* Output snapshots in separated files */
void output_snap(const Node *bird, int step)
{
	char file[100];
	snprintf(file, 100, "snap%ss_%s_%08d.bin", dlm.c_str(), para, step);
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
	buff = nullptr;
}

/* Output snapshots in one single file */
void output_snap(const Node *bird)
{
	float *buff = new float[3 * Node::N];
	for (int j = 0; j < Node::N; j++)
	{
		buff[3 * j] = bird[j].x;
		buff[3 * j + 1] = bird[j].y;
		buff[3 * j + 2] = atan2(bird[j].vy, bird[j].vx);
	}
	fout_snap.write((char *)&buff[0], sizeof(float) * Node::N * 3);
	delete[] buff;
	buff = nullptr;
}

void output(const Node *bird, int step)
{
	if (step % 100 == 0)
	{
		output_phi(bird, step);
		if (step % dt_snap == 0)
		{
			if(flag_one_snap_file)
				output_snap(bird);
			else
				output_snap(bird, step);
			if (step % dt_log == 0)
				output_log(step);
		}
	}
}