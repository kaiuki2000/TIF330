#include "primitives.h"

class ensemble
{
public:
    //virtual void allocate(int Size) = 0;
    virtual int addParticle(vector3d r, vector3d p) = 0;
    // interface for going through all particles in the ensemble
    virtual void begin(vector3d &r, vector3d &p) = 0;
    virtual bool next(vector3d &r, vector3d &p) = 0;
    virtual void setCurrent(vector3d r, vector3d p) = 0;
};

class ensembleSoA : public ensemble
{
public:
    vector<double> x; //coordinates
    vector<double> y;
    vector<double> z;
    vector<double> px; //momentum components
    vector<double> py;
    vector<double> pz;

    int number; // number of particles
    int size;

    int index; // for the loop
    ensembleSoA(int Size) : size(Size) // the statement after ":" is another way to assign values of variables during the creation of an object
    {
        x.resize(size);
        y.resize(size);
        z.resize(size);
        px.resize(size);
        py.resize(size);
        pz.resize(size);
        number = 0; // the allocated errays are empty and are supposed to be filled with addParticle()
        index = 0;
    }
    int addParticle(vector3d r, vector3d p)
    {
        if(number == x.size())
        {
            cout << "ensembleSoA: allocated vector is full." << endl;
            return 1; // error code that indicates running out of allocated number of particles
        }
        
        x[number] = r.x;
        y[number] = r.y;
        z[number] = r.z;
        
        px[number] = p.x;
        py[number] = p.y;
        pz[number] = p.z;
        
        number++;
        
        return 0;
    }
    void begin(vector3d &r, vector3d &p)
    {
        index = 0;
        r = vector3d(x[index], y[index], z[index]);
        p = vector3d(px[index], py[index], pz[index]);
    }
    bool next(vector3d &r, vector3d &p)
    {
        index++;
        if(index == number)
        {
            index = 0;
            return false;
        }
        r = vector3d(x[index], y[index], z[index]);
        p = vector3d(px[index], py[index], pz[index]);
        return true;
    }
    void setCurrent(vector3d r, vector3d p)
    {
        x[index] = r.x;
        y[index] = r.y;
        z[index] = r.z;
        px[index] = p.x;
        py[index] = p.y;
        pz[index] = p.z;
    }
};

class density
{
public:
    // x \in [0, 1], y \in [0, 1] 
    vector<double> data;
    int Nx;
    int Ny;
    double stepX;
    double stepY;
    density(int Nx, int Ny): Nx(Nx), Ny(Ny)
    {
        stepX = 1/double(Nx); 
        stepY = 1/double(Ny);
        data.resize(Nx*Ny);
    }
    void Null()
    {
        memset(&data[0], 0, data.size() * sizeof(double));    
    }
    void weighParticle(vector3d r, double weight) // NGP
    {
        if((abs(r.z - 0.5) < 0.1)&&(abs(r.x - 0.5) < 0.5)&&(abs(r.y - 0.5) < 0.5))
            data[int(r.x*Nx) + Nx*int(r.y*Ny)] += weight;
    }
    void appendFile(string fileName)
    {  
        FILE * pFile;
        pFile = fopen (fileName.c_str(), "ab");  
        fwrite(&data[0], sizeof(double), Nx*Ny, pFile);        
        fclose (pFile);
    }  
};

