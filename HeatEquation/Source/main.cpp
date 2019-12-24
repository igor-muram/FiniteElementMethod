#include "Core.h"

int main()
{
	FEM fem(
		"points.txt",
		"triangles.txt",
		"boundary1.txt",
		"boundary2.txt");
	system("pause");
	return 0;
}