#pragma once
#include "Matrix.h"
#include <vector>

namespace Solvers
{
	int LOS(Matrix& A, vector<double>& x, vector<double>& b);
	int BCG(Matrix& A, vector<double>& x, vector<double>& b);
}