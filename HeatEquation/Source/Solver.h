#pragma once
#include "Matrix.h"
#include <vector>

namespace Solvers
{
	int BCG(Matrix& A, std::vector<double>& x, std::vector<double>& b);
}