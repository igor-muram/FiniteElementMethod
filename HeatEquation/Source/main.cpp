#include "Core.h"

int main()
{
	FEM fem(
		"C:/data/points.txt",
		"C:/data/triangles.txt",
		"C:/data/boundary1.txt",
		"C:/data/boundary2.txt");
	system("pause");
	return 0;
}