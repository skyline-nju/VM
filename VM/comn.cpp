#include "comn.h"

using namespace std;

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
