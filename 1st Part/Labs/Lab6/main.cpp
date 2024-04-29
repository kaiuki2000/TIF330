#include "emfield.h"

int main()
{
    cout << "lab-session-5" << endl;    
    
    double timeStep = 0.05;
    emField field(64, 64, 64, 1, 1, 1, timeStep);
    string Name = "out1";
    field.initFile(Name);

    int particleNumber = 10000;
    ensembleSoA Ensemble(particleNumber);
    for(int i = 0; i < Ensemble.size; i++) 
        Ensemble.addParticle(vector3d(0.5 + 0.2*rand()/double(RAND_MAX), 0.5 + 0.2*rand()/double(RAND_MAX), 0.5 + 0.2*rand()/double(RAND_MAX)), vector3d(0, 0, 0));

    density Density(64, 64);

    vector3d r, p;
    for(int i = 0; i < 10; i++)
    {
        //output for the field
        field.appendFile();

        Density.Null();
        for(Ensemble.begin(r, p); Ensemble.next(r, p);)
        {
            Density.weighParticle(r, 1000/double(particleNumber));
        }
        Density.appendFile(Name + ".dat");

        // advance particles
        for(Ensemble.begin(r, p); Ensemble.next(r, p);)
        {
            vector3d E, B, v;
            field.getField_NGP(r, E, B);
            v = (1.0/sqrt(1 + p.norm2()))*p;
            r += timeStep*v;
            p += timeStep*100*E;
            p += (-timeStep*100)*vectorProduct(v, B);
            Ensemble.setCurrent(r, p);
        }
        
        field.advance();
        cout << i << endl;
    }
}