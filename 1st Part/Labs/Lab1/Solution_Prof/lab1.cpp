#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cmath>

using namespace std;

struct FTCS1D_solver // the structure that contains data and parameters of the simulation
{
    vector<double> u; // vector with componetns being the function values at grid nodes 
    double timeStep, spaceStep;
    double omega0, omega1;
    string fileName; // file name for the output
    double t; // current time
    
    // in the following block we descrie so-called contructor, a function that is called to creat a new object of our class 
    FTCS1D_solver(int gridSize, double timeStep, double omega0 = 2, double omega1 = 3): 
    u(gridSize, 0), timeStep(timeStep), omega0(omega0), omega1(omega1){
        t = 0;
        spaceStep = 1.0/double(u.size() - 1);
    }

    // in the following block we descibe the routine for updating the sate of u by a single time step
    void advance()
    {
        double u_ = u[0], u_old = 0; // we declare two auxiliary numbers
        for(int m = 1; m < u.size() - 1; m++)
        {
            u_old = u[m];
            u[m] += 0.1*timeStep*(u_ - 2*u[m] + u[m+1])/(spaceStep*spaceStep);
            u_ = u_old;
        }
        t += timeStep;
        u[0] = u[1] - spaceStep*sin(omega0*t); // left boundary condition
        u[u.size() - 1] = u[u.size() - 2] + spaceStep*sin(omega1*t); // right boundary condition
    }

    void initOutput(string experimentName) // function that creates or clears the output files
    { 
        fileName = experimentName;
        ofstream file;
        file.open((fileName + ".dat").c_str(), ios::out);
        file.close();
        file.open((fileName + ".txt").c_str(), ios::out);
        file << "0" << endl << u.size();
        file.close();
    }

    void appendOutput() // function that appends the output file with the data for the currect state of u
    {
        ifstream myfile;
        myfile.open((fileName + ".txt").c_str());
        int image;
        myfile >> image;
        image++;
        ofstream file;
        file.open((fileName + ".txt").c_str(), ios::out);
        file << image << endl << u.size();
        file.close();
        FILE * pFile;
        pFile = fopen ((fileName + ".dat").c_str(), "ab");        
        fwrite (&u[0], sizeof(double), u.size(), pFile);
        fclose (pFile);
    }
};

double sqr(double x){ // auxiliary function
    return x*x;
};

int main()
{    
    int gridSize = 2048;
    FTCS1D_solver sim(gridSize, 4/sqr(gridSize), sqrt(2), sqrt(3));
    long int numberOfIterations = 10.0/sim.timeStep;
    cout << "gridSize = " << gridSize << endl;
    cout << "spaceStep = " << sim.spaceStep << endl;
    cout << "timeStep = " << 4/sqr(gridSize) << endl;
    cout << "numberOfIterations = " << numberOfIterations << endl;
    for(int m = 0; m < sim.u.size(); m++) sim.u[m] = sqr(sin(M_PI*m/double(sim.u.size())));
    sim.initOutput("e1");
    int numberOfImages = 100;
    int stride = int(numberOfIterations/double(numberOfImages));
    for(long int i = 0; i < numberOfIterations; i++){
        sim.advance();
        if(i%stride == 0){
            cout << 100*(i/double(numberOfIterations)) << "%" << endl;
            sim.appendOutput();
        }
    }
    int i_ = int(0.3/sim.spaceStep);    
    double f = (0.3/sim.spaceStep - i_);
    cout << sim.u[i_]*(1 - f) + sim.u[i_+1]*f << endl;
}