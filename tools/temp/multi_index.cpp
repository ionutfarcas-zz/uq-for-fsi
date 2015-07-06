#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

typedef std::vector<std::vector<int>> vec2d;

int l1_norm(const std::vector<int> mindex)
{
	int norm = 0;

	for (int elem : mindex)
	{ 
		norm += abs(elem); 
	}

	return norm;
}

vec2d mindex(const int& dimension, const int& degree)
{
	int j = 0;
	int norm = 0;

	std::vector<int> temp(dimension, 0);
	vec2d mindex_degree;

	while(true)
	{
		norm = l1_norm(temp);

		if(norm <= degree)
		{
			mindex_degree.push_back(temp);
		}

		for(j = dimension - 1 ; j >= 0 ; --j)
		{
			if(++temp[j] <= degree)
				break;
			else
				temp[j] = 0;
		}

		if( j < 0)
			break;
	}

	return mindex_degree;
}

int main()
{
	int dimension = 8;
	int degree = 6;
	int size = 0;

	vec2d mindex_degree = mindex(dimension, degree);
	size = mindex_degree.size();

	std::cout << "size is " << size << std::endl;

	for(int i  = 0 ; i < size ; ++i)
	{
		for(int j = 0 ; j < dimension ; ++j)
		{
			std::cout << mindex_degree[i][j] << " ";
		}
		std::cout << std::endl;
	}

	return 0;
}