#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main()
{
// Definitions: grid dimensions, etc.
float h     = 0.1, tau = 0.00005;
float x_min = -1 , x_max = 1;
float t_min = 0  , t_max = 2.5;
int   N_x_E = floor((x_max - x_min)/(2*h));
int   N_x_B = floor((x_max - (x_min + h))/(2*h));
int   N_t_E = floor((t_max - t_min)/tau);
int   N_t_B = floor((t_max - (t_min + tau))/tau);
int   N_t   = max(N_t_E, N_t_B);

vector<float> x_E, x_B, t_E, t_B;
// Filling up space grid vector (E_y)
for(int i = 0; i < N_x_E; i++){
    x_E.push_back(0 * h + i * h);
}
// Filling up space grid vector (B_z)
for(int i = 0; i < N_x_B; i++){
    x_B.push_back(1 * h + i * h);
}
// Filling up time grid vector (E_y)
for(int i = 0; i < N_t_E; i++){
    t_E.push_back(0 * tau + i * tau);
}
// Filling up time grid vector (B_z)
for(int i = 0; i < N_t_B; i++){
    t_B.push_back(1 * tau + i * tau);
}

// Initial condition
vector<double> E_0;
for(int i = 0; i < N_x_E; i++){
    if(abs(x_E[i]) < 0.1){
        E_0.push_back(sin(20 * M_PI * x_E[i]));
    }
    else{
    E_0.push_back(0);
    }
}

// Remember that we want to implement periodic boundary conditions!
// 1D FDTD Implementation:
vector<vector<double>> E, B;
E.push_back(E_0);
for(int i = 0; i < N_t; i++)
{
    // First time instant is different. All (E and B) fields are 0 before the source term (in E) comes to exist.
    if(i == 0){
        // Evolve B
        vector<double> B_Next;
        B_Next.push_back(0 - tau/h*(E[int(i/2)][1] - E[int(i/2)][0]));
        continue;
    }

    if(i%2 == 0){
        // Evolve B
        vector<double> B_Next;
        for(int j = 0; j < N_x_B - 1; j++)
            B_Next.push_back(B[int(i/2) - 1][j] - tau/h*(E[int(i/2)][j + 1] - E[int(i/2)][j]));
        if(N_x_B == N_x_E)
            B_Next.push_back(B[int(i/2) - 1][N_x_B - 2] - tau/h*(E[int(i/2)][0] - E[int(i/2)][N_x_E - 1]));  // j = N_x_B - 1
        B.push_back(B_Next);
    }
    else{
        // Evolve E
        vector<double> E_Next;
        for(int j = 1; j < N_x_E - 1; j++)
            E_Next.push_back(E[floor(i/2)][j] - tau/h*(B[floor(i/2) + 1][j] - B[floor(i/2) + 1][j - 1]));
        E_Next.insert(E_Next.begin(), E[floor(i/2)][0] - tau/h*(B[floor(i/2) + 1][0] - B[floor(i/2) + 1][N_x_E - 1])); // j = 0
        if(N_x_B == N_x_E)
            E_Next.push_back(E[floor(i/2)][N_x_E - 2] - tau/h*(B[floor(i/2) + 1][0] - B[floor(i/2) + 1][N_x_B - 1])); // j = N_x_E - 1
        E.push_back(E_Next);

    }
}

// Writing to text file
fstream myfile_E;
myfile_E.open("Solution_E.txt", fstream::out);
for (int i = 0; i < N_t; i++)
{     
    for (int j = 0; j < N_x_E; j++)
    {
    if(j < (N_x_E - 1))
        myfile_E << E[i][j] << ",";
    else
        myfile_E << E[i][j] << endl;
    }
}

fstream myfile_B;
myfile_B.open("Solution_B.txt", fstream::out);
for (int i = 0; i < N_t; i++)
{     
    for (int j = 0; j < N_x_B; j++)
    {
    if(j < (N_x_B - 1))
        myfile_B << B[i][j] << ",";
    else
        myfile_B << B[i][j] << endl;
    }
}

// Message
cout << "Successfully printed to \"Solution_E.txt\" and \"Solution_B.txt\"." << endl;

return(0);
}