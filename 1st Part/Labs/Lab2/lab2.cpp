#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

void Triangle_Vec(vector<vector<double>> &v, int n)
{   
    /* First line */
    vector<double> Temp(n, 0);
    v.push_back(Temp);
    n--;
    /* Other lines */
    while (n >= 1)
    {
        vector<double> Temp(n, 0);
        v.push_back(Temp);
        v.push_back(Temp);
        n--;
    }
}

/* Auxiliary function */
void Print_Matrix(vector<vector<double>> &M)
{
    /* Printing matrix */
    for (int i = 0; i < M.size(); i++)
    {
        for (int j = 0; j < M[i].size(); j++)
        {
            cout << M[i][j] << " ";
        }   
    cout << endl;
    }
}

/* Might be useful later 
int Check_Pos(int i, int j, vector<vector<double>> &testVector)
{
   try
    {
       testVector.at(i).at(j);
       return(testVector[i][j]);
    }
   catch (const out_of_range& oor)
    {
        return(0);
    }
} */

bool Check_Bc(int i, int j, vector<vector<double>> &testVector)
{
   try
    {
        if(i%2 == 1)
        {
        testVector.at(i+1).at(j);
        testVector.at(i-1).at(j);
        testVector.at(i-1).at(j+1);
        return(false);
        }
        else
        {
        testVector.at(i+1).at(j);
        testVector.at(i+1).at(j-1);
        testVector.at(i-1).at(j);
        return(false);
        }
    }
   catch (const out_of_range& oor)
    {
        return(true);
    }
}

int main(){
    vector<vector<double>> r, u;
    int    n     = 100;
    double h     = 1/double(n);
    double tau   = 1e-5;
    double S     = sqrt(3)/4 * (h*h);
    double l     = 2 * sqrt(3)/6 * h;
    int    t_max = 100;
    Triangle_Vec(u, n);

    /* Initial condition */
    u[n/4 + 1][n/2 - 1] = 0.1;
    Print_Matrix(u);

    /* Looping over time */
    for (int t = 0; t < t_max; t++)
    {
        /* Update equation */
        for (int i = 0; i < u.size(); i++)
        {
            for (int j = 0; j < u[i].size(); j++)
            {
                if(Check_Bc(i, j, u))
                {
                    u[i][j] = 0;
                }
                else
                {
                    /* General case: Up/Down triangles must be treated differently! */
                    if(i%2 == 1)
                        u[i][j] = u[i][j] + (h * tau)/(l * S) * (u[i+1][j] + u[i-1][j] + u[i-1][j+1] - 3*u[i][j]);
                    else
                        u[i][j] = u[i][j] + (h * tau)/(l * S) * (u[i+1][j] + u[i+1][j-1] + u[i-1][j] - 3*u[i][j]);
                }
             }
        }
    }
    
    Print_Matrix(u);

    // Writing to text file
    fstream myfile;
    myfile.open("drum_sol.txt", fstream::out);
    for (int i = 0; i < u.size(); i++)
    {     
        for (int j = 0; j < u[i].size(); j++)
        {
        if(j < (u[i].size() - 1))
            myfile << u[i][j] << ",";
        else
            myfile << u[i][j] << endl;
        }
    }

    return(0);
}