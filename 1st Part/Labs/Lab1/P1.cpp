#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main()
{
// Definitions
float h = 0.1, tau = 0.00005, sigma = tau/(h*h); // Random values so far.
int   N_x = int(1/h);
int   N_t = int(10/tau);
vector<double> x, t;
int OMEGA_0 = 1, OMEGA_1 = 1;
float a = 0.1;

// From Von Neumann stability analysis
cout << "Parabolic Courant number (Should be ≤ 1/2) = " << sigma << " [Check: " << bool(sigma < 0.5) << "]" << endl;

// Filling up space grid vector
for(int i = 0; i < N_x; i++){
    x.push_back(i * h);
}
// Filling up time grid vector
for(int i = 0; i < N_t; i++){
    t.push_back(i * tau);
}

/*
// Printing: Debugging
// Space grid
cout << "Space grid:\n";
for(auto a : x)
    cout  << a << " ";
cout << endl;

// Time grid
cout << "Time grid:\n";
for(auto b : t)
    cout << b << " ";
cout << endl;
*/

// Initial condition and boundary conditions
vector<double> u0, u_p0, u_p1;
for(int i = 0; i < N_x; i++){
    u0.push_back(pow(sin(M_PI * x[i]), 2));
}
for(int i = 0; i < N_t; i++){
    u_p0.push_back(sin(OMEGA_0 * t[i]));
    u_p1.push_back(sin(OMEGA_1 * t[i]));
}

// Solution vector
vector<vector<double>> u;
for(int i = 0; i < N_t - 1; i++)
{
    if(i == 0)
    {
        vector<double> v = u0;
        u.push_back(v);
    }
    
    vector<double> v_Next;

    v_Next.push_back(u[i][1] - h * u_p0[i]);
    for(int j = 1; j < N_x - 1; j++)
        v_Next.push_back(u[i][j] + a * tau/(h*h) * (u[i][j - 1] - 2*u[i][j] + u[i][j + 1]));
    v_Next.push_back(u[i][x.size() - 2] + h * u_p1[i]);

    u.push_back(v_Next);
}

// Writing to text file
fstream myfile;
myfile.open("Solution.txt", fstream::out);
for (int i = 0; i < N_t; i++)
{     
    for (int j = 0; j < N_x; j++)
    {
    if(j < (N_x - 1))
        myfile << u[i][j] << ",";
    else
        myfile << u[i][j] << endl;
    }
}

// Message
cout << "Successfully printed to \"Solutions.txt\"." << endl;

return(0);
}