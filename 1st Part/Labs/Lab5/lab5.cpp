#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <stdio.h>

using namespace std;

void pyPlot2D(double** X, int Nx, int Ny, string fileName = "", string kwarg = "", bool noFrame = false, bool positive = false) // X[iy][ix]
{
    ofstream file; file.open("plot2D_tmp.dat", ios::out); file.close();
    FILE * dataFile = fopen ("plot2D_tmp.dat", "ab");
    fwrite (&Ny, sizeof(int), 1, dataFile);
    fwrite (&Nx, sizeof(int), 1, dataFile);
    for(int y = 0; y < Ny; y++) for(int x = 0; x < Nx; x++) fwrite (&(X[y][x]), sizeof(double), 1, dataFile);
    fclose (dataFile);
    ofstream pyFile; pyFile.open("plot2D_tmp.py", ios::out); 
    pyFile << "import matplotlib.pyplot as plt" << endl;
    pyFile << "import numpy as np" << endl; 
    pyFile << "import struct" << endl;
    pyFile << "f = open('plot2D_tmp.dat', 'rb')" << endl;
    pyFile << "Ny = struct.unpack('I', f.read(4))[0]" << endl;
    pyFile << "Nx = struct.unpack('I', f.read(4))[0]" << endl;
    pyFile << "v = np.empty(shape=(Ny, Nx))" << endl;
    pyFile << "for y in range(Ny):" << endl;
    pyFile << "    for x in range(Nx):" << endl;
    pyFile << "        v[y][x] = struct.unpack('d', f.read(8))[0]" << endl;
    pyFile << "fig, ax = plt.subplots()" << endl;
    
    pyFile << "vAbsMax = np.maximum(abs(np.amin(v)), abs(np.amax(v)))" << endl;
    if(!positive) kwarg += ", vmin = -vAbsMax, vmax = vAbsMax";
        else kwarg += ", vmin = 0, vmax = vAbsMax";

    pyFile << "plot = ax.imshow(v, interpolation='none', extent=[0,1,0,1], origin='lower'" << kwarg << ")" << endl;
    pyFile << "ax.set_aspect('equal')" << endl;
    pyFile << "fig.colorbar(plot, ax=ax, location='right')" << endl;
    //pyFile << "ax.contour(v, np.array([0]), colors='k', extent=[0,1,0,1])" << endl;
    //if(fileName == "") pyFile << "plt.show()" << endl;
    //    else pyFile << "plt.savefig('" << fileName << ".png')" << endl;
    if(noFrame){
        pyFile << "plt.gca().set_axis_off()" << endl;
        pyFile << "plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)" << endl;
        pyFile << "plt.margins(0,0)" << endl;
        pyFile << "plt.gca().xaxis.set_major_locator(plt.NullLocator())" << endl;
        pyFile << "plt.gca().yaxis.set_major_locator(plt.NullLocator())" << endl;
        pyFile << "plt.savefig('" << fileName << ".png', bbox_inches = 'tight', pad_inches = 0)" << endl;
    } else 
        pyFile << "plt.savefig('" << fileName << ".png')" << endl;
    
    pyFile.close();
    system("python plot2D_tmp.py");
};

void test_plot2D()
{
    cout << "test_plot2D" << endl;
    int Nx = 128, Ny = 64;
    double **X = new double*[Ny]; for(int i = 0; i < Ny; i++) X[i] = new double[Nx];
    for(int iy = 0; iy < Ny; iy++)
    for(int ix = 0; ix < Nx; ix++){
        double y = iy/double(Ny), x = ix/double(Nx);
        X[iy][ix] = sin(2*M_PI*y) + 0.0*sin(2*M_PI*(x + 2*y));
    }
    pyPlot2D(X, Nx, Ny, "test_plot2d", ", cmap='RdBu'");
};

void plot1D(vector<double> X, string fileName = "", double yRange = 0) // function for plotting a vector<double>
{
    ofstream file; file.open("plot1D_tmp.dat", ios::out); file.close();
    FILE * dataFile = fopen ("plot1D_tmp.dat", "ab");
    int N = X.size(); fwrite (&N, sizeof(int), 1, dataFile);
    for(int i = 0; i < X.size(); i++) fwrite (&(X[i]), sizeof(double), 1, dataFile);
    fclose (dataFile);
    ofstream pyFile; pyFile.open("plot1D_tmp.py", ios::out); 
    pyFile << "import matplotlib.pyplot as plt" << endl;
    pyFile << "import struct" << endl;
    pyFile << "f = open('plot1D_tmp.dat', 'rb')" << endl;
    pyFile << "N = struct.unpack('I', f.read(4))[0]" << endl;
    pyFile << "y = [0]*N" << endl;
    pyFile << "for i in range(N):" << endl;
    pyFile << "    y[i] = struct.unpack('d', f.read(8))[0]" << endl;
    pyFile << "fig, ax = plt.subplots()" << endl;
    pyFile << "ax.plot(y)" << endl;
    if(yRange != 0) pyFile << "ax.set_ylim(" << -yRange << ", "<< yRange << ")" << endl;
    if(fileName == "") pyFile << "plt.show()" << endl;
        else pyFile << "plt.savefig('" << fileName << ".png')" << endl;
    pyFile.close();
    system("python plot1D_tmp.py");
};
double randNormal(double sigma2) // generation of a pseudorandom number following normal distribution.
{ 
    return cos(2*M_PI*(rand()/double(RAND_MAX)))*sqrt(-2*sigma2*log(1 - (1 + rand())/double((long)RAND_MAX + 2)));
}
//CGS units are used
const double electronCharge = -4.80320427e-10;
const double electronMass = 9.10938215e-28;
const double lightVelocity = 29979245800.0;
double sqr(double x){return x*x;};

struct electrostaticPIC1D{
    vector<double> Ex, Jx;
    vector<double> px, x;
    double L, dt, weight;
    int nCells, nParticles;
    electrostaticPIC1D(double L, double dt, int nParticles, int nCells): L(L), dt(dt), nParticles(nParticles), nCells(nCells)
    {
        // definitions/constants
        const double a       = 1e-1;
        const double Ne      = 1e18;
        const double T       = 1e-6*electronMass*sqr(lightVelocity); 
        const double LambdaD = sqrt(T/(4*M_PI*Ne*sqr(electronCharge)));
        weight = L*Ne/(nParticles);

        Ex.resize(nCells);
        Jx.resize(nCells);
        px.resize(nParticles);
        x.resize(nParticles);

        // initialization
        for(int i = 0; i < px.size(); i++)
        {
            x[i]  = L*(rand() + 1)/double((long)RAND_MAX + 2);
            px[i] = randNormal(2 * electronMass * T); // (Sigma^2 = 2mT) here we need to generate a random value to get the normal distribution with a given temperature; function randNormal() can be helpful 
        }
        for(int i = 0; i < Ex.size(); i++)
            Ex[i] = a * 4*M_PI*L*electronCharge*Ne*sin(2*M_PI*(i*LambdaD)/L);
    }
    void advance()// this function is to advance the state of both particles and fields
    {
        for(int ix = 0; ix < Jx.size(); ix++)
            Jx[ix] = 0;

        //compute Jx
        for(int ip = 0; ip < px.size(); ip++){
            // here we need to deposite the contribution of ip-th particle to the grid value of Jx at its location
            Jx[floor(x[ip]/L * nCells)] += px[ip]/electronMass*electronCharge*(weight)/(L/nCells); // was missing weighing the macro-particles
        }   

        //advance Ex
        for(int ix = 0; ix < Ex.size(); ix++){
            // here we need to do one step for Ex[ix] to advance its state in time
            Ex[ix] += -4*M_PI*Jx[ix]*dt;
        }

        //advance particles
        for(int ip = 0; ip < px.size(); ip++){
            //here we need to advance the state of particles
            x[ip]  += px[ip]/electronMass*dt;
            while (x[ip] >= L){x[ip] -= L;};
            while (x[ip] < 0 ){x[ip] += L;};            
            px[ip] += electronCharge*Ex[floor(x[ip]/L * nCells)]*dt;
        }
    }
    void plot_xpx(double temperature, int imageNumber)
    {
        int nx = 100;
        int npx = 100;
        double pxMax = 30*sqrt(temperature*electronMass);
        double **d = new double*[npx];
        for(int i = 0; i < npx; i++){
            d[i] = new double[nx];
            for(int ix = 0; ix < nx; ix++)d[i][ix] = 0;
        }

        for(int ip = 0; ip < px.size(); ip++){
            int ix = int(nx*x[ip]/L + 0.5) % nx;
            int ipx = int(npx*(0.5 + px[ip]/pxMax));
            if((ipx > 0)&&(ipx < npx)){
                d[ipx][ix] += 1;
            }
        }
        pyPlot2D(d, nx, npx, "xpx" + to_string(imageNumber), ", cmap='plasma'", false, true);
        for(int i = 0; i < npx; i++) delete []d[i];
        delete []d;
    }
    double GetField_energy()
    {
        int Sum = 0;
        for(int ix = 0; ix < Ex.size(); ix++)
            Sum += sqr(Ex[ix]);
        return(Sum);
    }
};

int main()
{
    // here we need to create the instance of sruct "electrostaticPIC1D", set up initial conditions and call advance() to advance its state.
    cout << "lab5" << endl;
    // definitions/constants
    const double Ne      = 1e18;
    const double T       = 1e-6*electronMass*sqr(lightVelocity); 
    const double LambdaD = sqrt(T/(4*M_PI*Ne*sqr(electronCharge)));
    
    // relevant values
    const double L       = 100 * LambdaD;
    const double dt      = 1e-3 * sqrt(M_PI*electronMass/(Ne*sqr(electronCharge)));
    const int nCells     = floor(L/LambdaD);    
    const int nParticles = 1e3*nCells;

    double R1, R2;

    electrostaticPIC1D Obj(L, dt, nParticles, nCells);

    R1 = Obj.GetField_energy();
    for (int i = 0; i < 10000; i++)
    {
        cout << "Iteration " << i << "..." << endl;
        if(i % 100 == 0){
            char* filename;
            sprintf(filename, "ExFig%d", i);
            plot1D(Obj.Ex, filename);
            Obj.plot_xpx(T, i);}
        Obj.advance();
    }
    R2 = Obj.GetField_energy();
    cout << "Ratio = " << R2/R1 << endl;
    return(0);
}