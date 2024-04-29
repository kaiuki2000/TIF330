#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <stdio.h>
#include <direct.h>

using namespace std;

const double pi = 3.14159265358979323846;

double sqr(double x)
{
    return x*x;
}


struct vector3d 
{
    double x, y, z;

    vector3d():
        x(0), y(0), z(0)
    {}

    vector3d(double x, double y, double z):
        x(x), y(y), z(z)
    {}

    vector3d operator + (const vector3d & v) const
    {
        return vector3d(x + v.x, y + v.y, z + v.z);
    }

    vector3d& operator += (const vector3d & v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    vector3d operator - (const vector3d & v) const
    {
        return vector3d(x - v.x, y - v.y, z - v.z);
    }

    double norm2() const
    {
        return x * x + y * y + z * z;
    }

};

inline vector3d vectorProduct(vector3d a, vector3d b)
{
    return vector3d(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
};

inline vector3d operator*(double a, vector3d b)
{
    return vector3d(a*b.x, a*b.y, a*b.z);
};
