#include "matrix.h"
#include <iostream>

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		std::cerr << "Missing arguments\n";
		return 1;
	}

	apr::matrix A(argv[1]);
	std::cout << "A =\n" << A;

	try
	{
		auto inverse = !A;
		std::cout << "\nA^-1 =\n" << inverse;
	} catch (...) {
		std::cout << "\nSingular matrix!\n";
	}
}
