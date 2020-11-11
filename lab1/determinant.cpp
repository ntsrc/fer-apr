#include "matrix.h"
#include <iostream>

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "Missing arguments\n";
		return 1;
	}

	apr::matrix A(argv[1]);
	std::cout << "A =\n" << A;

	auto det_A = apr::determinant(A);
	std::cout << "\ndet(A) =\n" << det_A << '\n';
}
