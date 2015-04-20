#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

//user-defined constants
const double TOP_TEMP = 100.;
const double LEFT_TEMP = 75.;
const double RIGHT_TEMP = 50.;
const double INTERIOR_NODES = 49.;
const double LAMBDA = 1.5;

//calculated
const double DIM = sqrt(INTERIOR_NODES) + 2;

//error function
double error(double old, double next)
{
	return abs((next - old) / next);
}

//performs one iteration of Gauss-Seidel on the grid and returns the maximum error
double liebmann(std::vector<std::vector<double>>& grid)
{
	double old;
	double next;
	double thisError;
	double maxError = 0.;
	for (int i = DIM-1; i > 0; i--)
	{ 
		for (int j = 1; j < DIM - 1; j++)
		{
			//recalculated interior nodes
			if (i != DIM-1)
			{
				old = grid[i][j];
				next = (grid[i][j+1] + grid[i][j-1] + grid[i+1][j] + grid[i-1][j]) / 4;
				next = LAMBDA*next + (1 - LAMBDA)*(old);
				thisError = error(old, next);
				grid[i][j] = next;
				if (thisError > maxError)
				{
					maxError = thisError;
				}
			}
			else
			{
				//recalculates insulated bottom boundary temperatures
				old = grid[i][j];
				next = (grid[i][j+1] + grid[i][j-1] + 2 * grid[i-1][j]) / 4;
				next = LAMBDA*next + (1 - LAMBDA)*(old);
				thisError = error(old, next);
				grid[i][j] = next;
				if (thisError > maxError)
				{
					maxError = thisError;
				}
			}
		}
	}
	return maxError;
}

int main()
{
	using namespace std;

	//bookkeeping variables
	double err;
	vector<vector<double>> rows(DIM+1);

	//populates grid
	for (int i = 0; i < DIM; i++)
	{
		if (i == 0)
		{
			for (int j = 0; j < DIM; j++)
			{
				rows[i].push_back(TOP_TEMP);
			}
		}
		else if (i > 0 && i < DIM-1)
		{
			rows[i].push_back(LEFT_TEMP);
			for (int j = 1; j < DIM - 1; j++)
			{
				rows[i].push_back(0);
			}
			rows[i].push_back(RIGHT_TEMP);
		}
		else
		{
			rows[i].push_back(LEFT_TEMP);
			for (int j = 1; j < DIM-1; j++)
			{
				rows[i].push_back(0);
			}
			rows[i].push_back(RIGHT_TEMP);
		}
	}
	
	//performs Gauss-Seidel until maximum error drops below 1%
	int count = 0;
	do{
		err = liebmann(rows);
		count++;
	} while (err > 0.01);

	//prints grid
	cout << "Liebmann's Method (insulated bottom boundary):\n\n";
	for (vector<double> v : rows)
	{
		for (double d : v)
		{
			cout << setw(7) << d << " ";
		}
		cout << endl;
	}
	cout << "Method accomplished in " << count << " iterations.\n";


	return 0;
}