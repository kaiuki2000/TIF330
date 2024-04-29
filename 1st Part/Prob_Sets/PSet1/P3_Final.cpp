#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

#define SIZE 2048

using namespace std;

/* Auxiliary functions */
void Write_to_file(vector<vector<double>> &Res, string filename)
{
    cout << "Entering \"Write_to_file\"... " << endl;
    cout << "Matrix dimensions: " << size(Res) << " " << size(Res[0]) << endl;
    fstream w_File;
    w_File.open(filename, fstream::out);
    for (int i = 0; i < size(Res); i++)
    {     
        for (int j = 0; j < size(Res[0]); j++)
        {
        if(j < (size(Res[0]) - 1)) w_File << Res[i][j] << ",";
        else                       w_File << Res[i][j] << endl;
        }
    }
    w_File.close();
}

void Write_Vec_to_file(vector<double> &Vec, string filename)
{
    cout << "Entering \"Write_Vec_to_file\"... " << endl;
    fstream w_File;
    w_File.open(filename, fstream::out);
    for (int i = 0; i < size(Vec); i++)
    {     
        if(i < (size(Vec) - 1)) w_File << Vec[i] << ",";
        else                    w_File << Vec[i] << endl;
    }
}

void Add_to_Res(double* E_y, double* B_z, vector<vector<double>> &Res_E, vector<vector<double>> &Res_B)
{
    vector<double> E_y_, B_z_;
    E_y_.assign(E_y, E_y + SIZE);
    B_z_.assign(B_z, B_z + (SIZE - 1));
    Res_E.push_back(E_y_);
    Res_B.push_back(B_z_);
    E_y_.clear();
    B_z_.clear();
}

/* Main program */
int main(int argc, const char** argv)
{
/* Boundary condition type (different cases):
1 - Periodic;
2 - Perfect Magnetic Counductor;
3 - Perfect (Infinity). */
for (int BC_type = 1; BC_type < 4; BC_type++)
{
    /* Variables */
    vector<double> t_vec;
    vector<vector<double>> Res_E, Res_B;
    int            qTime, mm;
    double         E_y[SIZE]       = {0.}, B_z[SIZE - 1] = {0.};
    double         h               = 2/(2 * double(SIZE) - 2);
    double         tau             = h;
    long int       maxTime         = 2.5/tau;
    double         x[2 * SIZE - 1] = {0.};
    bool           delay           = false;
    vector<string> names           = {"Periodic_BC", "PMC_BC", "Perfect_BC"};

    /* These ended up not being needed */
    int            numberOfImages  = 500; /* We only save 500 images */
    int            stride          = int(double(maxTime)/double(numberOfImages));

    /* Printing for debugging/reference */
    printf("Parameters:                \n");
    printf("h                      = %f\n", h);
    printf("gridSize               = %d\n", SIZE);
    printf("tau                    = %f\n", tau);
    printf("2h/tau                 = %f > c = 1 (Stable)\n", (2*h)/tau);
    printf("numberOfIterations     = %ld\n", maxTime);
    printf("numberOfImages         = %d\n", numberOfImages);
    printf("stride                 = %d\n", stride);

    /* Filling up space grid vector */
    for(mm = 0; mm < (2 * SIZE - 1); mm++) x[mm] = -1 + (mm * h);

    if(BC_type < 3)
    {
        /* Initial condition */
        for(mm = 0; mm < SIZE; mm++)
        {
            if(abs(x[2 * mm]) < 0.1) E_y[mm] = sin(20 * M_PI * x[2 * mm]);
            else                     E_y[mm] = 0;
        }
    }

    /* 'Pushing back' intial condition */
    t_vec.push_back(0.);
    Add_to_Res(E_y, B_z, Res_E, Res_B);

    /* Do time stepping */
    for (qTime = 0; qTime < maxTime; qTime++)
    {

        /* Update magnetic field. Was missing factors of (tau/h) here (and in the electric field equation(s))! */
        for (mm = 0; mm < (SIZE - 1); mm++) B_z[mm] = B_z[mm] - (tau/(2*h))*(E_y[mm + 1] - E_y[mm]);

        /* Update electric field */
        switch (BC_type)
        {
        case 1:
            /* Boundary conditions - Periodic! First and last points, respectively */
            E_y[0]        = E_y[0]        - (tau/h)*(B_z[0] - B_z[SIZE - 2]);
            E_y[SIZE - 1] = E_y[SIZE - 1] - (tau/h)*(B_z[0] - B_z[SIZE - 2]);
            break;
        case 2:
            /* Perfect magnetic conductor at edges */
            E_y[0]        = E_y[0]        - (tau/(2*h))*(B_z[0] - 0);
            E_y[SIZE - 1] = E_y[SIZE - 1] - (tau/(2*h))*(0 - B_z[SIZE - 2]);
            break;
        case 3:
            /* Perfect boundary conditions */
            if(qTime > 0)
            {
                E_y[0]        = E_y[0]        - (tau/(2*h))*(B_z[0] - Res_B[(qTime + 1) - 2][0]);
                E_y[SIZE - 1] = E_y[SIZE - 1] - (tau/(2*h))*(Res_B[(qTime + 1) - 2][SIZE - 2] - B_z[SIZE - 2]);
            }
            break;
        }

        /* All the other 'middle' points */
        for (mm = 1; mm < SIZE - 1; mm++) E_y[mm] = E_y[mm] - (tau/(2*h))*(B_z[mm] - B_z[mm - 1]);

        /* Delay counter: Antenna beeps and then waits for 500 iterations before beeping again */
        if(BC_type == 3)
        {
            if ((qTime+1) % 500 == 0) delay = !delay;
            if (!delay)               E_y[SIZE/2] += 0.75 * abs(sin(M_PI * (double(qTime)/500.)));
        }

        /* Results */
        if (qTime % 1 == 0) // 'stride' was changed to 1 for now.
        {
            cout << "Progress: " << 100*(qTime/double(maxTime)) << " %" << endl;
            t_vec.push_back((qTime + 1)*tau);
            Add_to_Res(E_y, B_z, Res_E, Res_B);
        }

    } /* End of time-stepping */

    /* Writing to text file */
    switch (BC_type)
    {
    case 1:
        Write_to_file(Res_E, "E_Periodic_BC.txt");
        Write_to_file(Res_B, "B_Periodic_BC.txt");
        break;
    case 2:
        Write_to_file(Res_E, "E_PMC_BC.txt");
        Write_to_file(Res_B, "B_PMC_BC.txt");
        break;
    case 3:
        Write_to_file(Res_E, "E_Perfect_BC.txt");
        Write_to_file(Res_B, "B_Perfect_BC.txt");
        break;
    }

    /* Additionally saving the time vector, for easier plotting later */
    // ofstream outFile("timeVec.txt"); 
    // for (const auto &e : t_vec) outFile << e << "\n"; // Find out why this doesn't work!
    Write_Vec_to_file(t_vec, "timeVec.txt");
}

return(0);
}