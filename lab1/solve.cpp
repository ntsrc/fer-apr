#include "matrix.h"
#include <cstring>
#include <iostream>
#include <utility>

int main(int argc, char **argv)
{
	if (argc < 4)
	{
		std::cerr << "Missing arguments\n";
		return 1;
	}

	if (strcmp(argv[1], "-LU") && strcmp(argv[1], "-LUP"))
		std::cerr << "Invalid arguments\n";

	apr::matrix A(argv[2]);
	apr::matrix b(argv[3]);

	std::cout << "A =\n" << A;
	std::cout << "\nb =\n" << b;

	if (!strcmp(argv[1], "-LUP"))
	{
		try 
		{
		auto [LUP, P] = apr::LUP_decomposition(A);
		std::cout << "\nLU =\n" << LUP;
		std::cout << "\nP =\n" << P;

		auto y = apr::forward_supstitution(LUP, P * b);
		std::cout << "\ny =\n" << y;

		auto x = apr::backward_supstitution(LUP, y);
		std::cout << "\nx =\n" << x;
		}
		catch (...) {
			std::cout << "\nCan't solve with LUP decomposition!\n";
		}
	} else {
		try
		{
		auto LU = apr::LU_decomposition(A);
		std::cout << "\nLU =\n" << LU;

		auto y = apr::forward_supstitution(LU, b);
		std::cout << "\ny =\n" << y;

		auto x = apr::backward_supstitution(LU, y);
		std::cout << "\nx =\n" << x;
		}
		catch (...) {
			std::cout << "\nCan't solve with LU decomposition!\n";
		}
	}
}
