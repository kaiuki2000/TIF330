#include "ensemble.h"

struct emCell
{
    vector3d E;
    vector3d B;
};

class emField
{
    // A class for storing and advancing electromagnetic field sampled with a 3D grid.
    // Physical limits: x \in [0, sizeX], y \in [0, sizeY], z \in [0, sizeZ];
    // Field values units: dimensionless (the speed of light is 1);
    // Boundaries: periodic (torus topology);
    // Grid size: sizeX, sizeY, sizeZ;
    // Allocation: AoS;
    int Nx, Ny, Nz; // grid sizes
    double sizeX, sizeY, sizeZ;
    double stepX, stepY, stepZ;
    double timeStep;
    double time;
    vector<emCell> field; //straight (collocated) grid
    int numberOfCells;
    string fileName;

private:
    int index(int x, int y, int z)
    {
        return x + Nx*(y + Ny*z); 
    }
    void nodeLocation(int x, int y, int z, vector3d &R)
    {
        R.x = (x + 0.25)*stepX;
        R.y = (y + 0.25)*stepY;
        R.z = (z + 0.25)*stepZ;
    }
public:
    emField(int Nx, int Ny, int Nz, double sizeX, double sizeY, double sizeZ, double timeStep): 
    Nx(Nx), Ny(Ny), Nz(Nz), sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), timeStep(timeStep)
    {
        stepX = sizeX/double(Nx);
        stepY = sizeY/double(Ny);
        stepZ = sizeZ/double(Nz);
        time = 0;
        numberOfCells = Nx*Ny*Nz;
        field.resize(numberOfCells);
        cout << "emField: the grid is allocated." << endl;
        calculateField(time);
    }
    void calculateField(double t) // analytical calculation of the field
    {
        vector3d R;
        double Size = 0.2;
        double amplitude;
        for(int z = 0; z < Nz; z++)
            for(int y = 0; y < Ny; y++)
                for(int x = 0; x < Nx; x++)
                {
                    field[index(x, y, z)].E = vector3d(0, 0, 0);
                    field[index(x, y, z)].B = vector3d(0, 0, 0);
                    nodeLocation(x, y, z, R);
                    amplitude = exp(-(sqr(R.x - t) + sqr(R.y - 0.5) + sqr(R.z - 0.5))/sqr(0.3));
                    field[index(x, y, z)].E.y = amplitude*sin(5*2*pi*(R.x - t + 1));
                    field[index(x, y, z)].B.z = field[index(x, y, z)].E.y;
                }
    }
    void advance()
    {
        time += timeStep;
        calculateField(time);
    }
    void appendFile()
    {        
        FILE * pFile;
        pFile = fopen (fileName.c_str(), "ab");  
        for(int y = 0; y < Ny; y++)
            for(int x = 0; x < Nx; x++)      
                fwrite(&field[index(x, y, Nz/2)].E.y, sizeof(double), 1, pFile);
                
        fclose (pFile);
    }
    void initFile(string experimentName)
    { 
        fileName = experimentName + ".dat";
        ofstream file;
        file.open(fileName.c_str(), ios::out);
        file.close();
        file.open((experimentName + ".txt").c_str(), ios::out);
        file << Nx << endl << Ny << endl;
        file.close();
    }
    void getField_NGP(vector3d r, vector3d &E_, vector3d &B_) // returns values of E and B at position r (NGP) 
    {
        E_ = field[index(int(r.x*Nx), int(r.y*Ny), int(r.z*Nz))].E; 
        B_ = field[index(int(r.x*Nx), int(r.y*Ny), int(r.z*Nz))].B;
    }
};