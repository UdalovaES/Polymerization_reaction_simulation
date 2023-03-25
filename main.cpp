#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <vector>
#include <random>


#define eps  3
#define sigma 1
#define m 1.0
#define L  10
#define r_c  4.0
#define p  0.001
#define lambda  10
#define kb  0.00831
#define tau 0.001
#define T  300
#define NSTEPS 7000

struct Atom
{
    float x, y, z;
    float vx, vy, vz;
    float fx, fy, fz;
    char name;
    bool real;
//    float q;
//    float sigma;
//    float epsilon;
//    char name;
};

float transferPBC(float x)
{
    if (x < 0)
    {
        return x + L;
    }
    else if (x > L)
    {
        return x - L;
    }
    return x;
}

void saveFrame(const char* filename, const char* modifier, std::vector<Atom> atoms)
{
    FILE* out = fopen(filename, modifier);
    fprintf(out, "%ld\nA+B->A+C\n", atoms.size());
    for (int i = 0; i < atoms.size(); i++)
    {
        fprintf(out, "%c\t%f\t%f\t%f\n",
                atoms[i].name,
                atoms[i].x*10.0,
                atoms[i].y*10.0,
                atoms[i].z*10.0);
    }
    fclose(out);
}

// float distance(Atom atom1, Atom atom2) {
//     float dist = pow((atom1.x - atom2.x), 2) +
//                  pow((atom1.y - atom2.y), 2) +
//                  pow((atom1.z - atom2.z), 2);
//     return sqrtf(dist);
// }

float distance(Atom atom1, Atom atom2) {
    float dx = atom1.x - atom2.x;
    float dy = atom1.y - atom2.y;
    float dz = atom1.z - atom2.z;

    // Check these
    dx -= rint(dx/L)*L;
    dy -= rint(dy/L)*L;
    dz -= rint(dz/L)*L;

    float dist = pow(dx, 2) +
                 pow(dy, 2) +
                 pow(dz, 2);

    return sqrtf(dist);
}

int main(int argc, char* argv[])
{
    int N_b = 100, N_a = 10;
    int N_c = N_b;
    int N_start = N_b + N_a;
    int N = N_b + N_c + N_a;
    // float currentTemperature = T;
    std::vector<Atom> atoms(N);

    std::random_device randomDevice;
    std::mt19937 randomGenerator(randomDevice());
    std::uniform_real_distribution<> distributionX(0, L);
    std::normal_distribution<> distributionV(0.0, sqrtf(kb*T));
    std::normal_distribution<> distributionR(0.0, 1.0);
    std::uniform_real_distribution<> distributionKSI(0, 1);

    Atom null_atom;
    null_atom.name = 'B';
    null_atom.x = 0;
    null_atom.y = 0;
    null_atom.z = 0;    
    null_atom.vx = 0;
    null_atom.vy = 0;
    null_atom.vz = 0;
    null_atom.fx = 0.0;
    null_atom.fy = 0.0;
    null_atom.fz = 0.0;
    null_atom.real = false;

    for (int i = 0; i < N; i++)
    {

        if (i < N_a) {
            atoms[i].name = 'A';
            atoms[i].real = true;
        } else if (i >= N_a and i < N_a + N_b) {
            atoms[i].name = 'B';
            atoms[i].real = true;
        } else {
            atoms[i].name = 'C';
            atoms[i].real = false;
        }
        
        if (i < N_start) {
            atoms[i].x = distributionX(randomGenerator);
            atoms[i].y = distributionX(randomGenerator);
            atoms[i].z = distributionX(randomGenerator);
            
            atoms[i].vx = distributionV(randomGenerator);
            atoms[i].vy = distributionV(randomGenerator);
            atoms[i].vz = distributionV(randomGenerator);

            atoms[i].fx = 0.0;
            atoms[i].fy = 0.0;
            atoms[i].fz = 0.0;

        } else {
            atoms[i].x = 0;
            atoms[i].y = 0;
            atoms[i].z = 0;

            atoms[i].vx = 0;
            atoms[i].vy = 0;
            atoms[i].vz = 0;

            atoms[i].fx = 0.0;
            atoms[i].fy = 0.0;
            atoms[i].fz = 0.0;
        }
    }
    saveFrame("ABC.xyz", "w", atoms);

    // double temperature = 0.0;
    // int nTemperature = 0;

    for (int n = 0; n < NSTEPS; n++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (atoms[i].real and atoms[j].real) 
                {
                    float rij = distance(atoms[i], atoms[j]);
                    if ((atoms[i].name == 'A' and atoms[j].name == 'B') or (atoms[i].name == 'B' and atoms[j].name == 'A')) {
                        if (rij < r_c) {
                            float ksi = distributionKSI(randomGenerator);
                            if (ksi < p) {
                                if (atoms[i].name == 'B') { 
                                    atoms[i + N_b] = atoms[i];
                                    atoms[i + N_b].name = 'C';
                                    atoms[i + N_b].real = true;
                                    atoms[i] = null_atom;
                                    }
                                if (atoms[j].name == 'B') {
                                    atoms[j + N_b] = atoms[j];
                                    atoms[j + N_b].name = 'C';
                                    atoms[j + N_b].real = true;
                                    atoms[j] = null_atom;
                                    }
                            }
                        }
                    }
                } 
            } 
        }
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < i; j++)
            {
                float rij = distance(atoms[i], atoms[j]);
                
                if (atoms[i].real and atoms[j].real) 
                {
                    float c1, c2;

                    if (atoms[i].name =='C' && atoms[j].name =='C'){
                        c1 = 1.0;
                        c2 = -3.0;
                    } else {
                        c1 = 0.0;
                        c2 = 2.0;
                    }

                    float fi = eps * 2 * pow(sigma, 2) * pow(rij, -4) * (2 * c1 * pow(sigma, 2) * pow(rij, -2) + c2);

                    float dx = atoms[i].x - atoms[j].x;
                    float dy = atoms[i].y - atoms[j].y;
                    float dz = atoms[i].z - atoms[j].z;

                    // Check these
                    dx -= rint(dx/L)*L;
                    dy -= rint(dy/L)*L;
                    dz -= rint(dz/L)*L;

                    atoms[i].fx += fi*dx;
                    atoms[i].fy += fi*dy;
                    atoms[i].fz += fi*dz;

                    atoms[j].fx -= fi*dx;
                    atoms[j].fy -= fi*dy;
                    atoms[j].fz -= fi*dz;

                // float dr2 = dx*dx + dy*dy + dz*dz;
                // float sigma = 0.5*(atoms[i].sigma + atoms[j].sigma);
                // float epsilon = sqrtf(atoms[i].epsilon*atoms[j].epsilon);
                // float sigma2 = sigma*sigma;
                // float sor2 = sigma2/dr2;
                // float sor6 = sor2*sor2*sor2;
                // float df = 12.0*epsilon*sor6*(sor6 - 1.0)/dr2;

                // float dr = sqrtf(dr2);
                // df += KC*atoms[i].q*atoms[j].q/(dr2*dr);

                // atoms[i].fx += df*dx;
                // atoms[i].fy += df*dy;
                // atoms[i].fz += df*dz;

                // atoms[j].fx -= df*dx;
                // atoms[j].fy -= df*dy;
                // atoms[j].fz -= df*dz;
                }
            }
        }

        for (int i = 0; i < N; i++)
        {
        if (atoms[i].real) 
        {
            // float scale = sqrtf(1.0 - ((currentTemperature-T)*tau)/(T*relax));
            // atoms[i].vx *= scale;
            // atoms[i].vy *= scale;
            // atoms[i].vz *= scale;

            // float mult = tau/atoms[i].m;

            // atoms[i].vx = atoms[i].vx + mult*atoms[i].fx;
            // atoms[i].vy = atoms[i].vy + mult*atoms[i].fy;
            // atoms[i].vz = atoms[i].vz + mult*atoms[i].fz;
// правильные
            atoms[i].vx = 2 * tau / (2 + lambda * tau) * (atoms[i].fx + atoms[i].vx*(m/tau - lambda*m/2) + 
            sqrtf(2*kb*T*lambda/tau) * distributionR(randomGenerator));
            atoms[i].vy = 2 * tau / (2 + lambda * tau) * (atoms[i].fy + atoms[i].vy*(m/tau - lambda*m/2) + 
            sqrtf(2*kb*T*lambda/tau) * distributionR(randomGenerator));
            atoms[i].vz = 2 * tau / (2 + lambda * tau) * (atoms[i].fz + atoms[i].vz*(m/tau - lambda*m/2) + 
            sqrtf(2*kb*T*lambda/tau) * distributionR(randomGenerator));
            // atoms[i].vx = tau / m * (atoms[i].fx + atoms[i].vx*(m/tau - lambda) + 
            // sqrtf(2*kb*T*lambda/tau) * distributionR(randomGenerator));
            // atoms[i].vy = tau / m * (atoms[i].fy + atoms[i].vy*(m/tau - lambda) + 
            // sqrtf(2*kb*T*lambda/tau) * distributionR(randomGenerator));
            // atoms[i].vz = tau / m * (atoms[i].fz + atoms[i].vz*(m/tau - lambda) + 
            // sqrtf(2*kb*T*lambda/tau) * distributionR(randomGenerator));

            atoms[i].x = atoms[i].x + tau*atoms[i].vx;
            atoms[i].y = atoms[i].y + tau*atoms[i].vy;
            atoms[i].z = atoms[i].z + tau*atoms[i].vz;

            atoms[i].x = transferPBC(atoms[i].x);
            atoms[i].y = transferPBC(atoms[i].y);
            atoms[i].z = transferPBC(atoms[i].z);

            atoms[i].fx = 0.0;
            atoms[i].fy = 0.0;
            atoms[i].fz = 0.0;
        }
        }

        if (n % 100 == 0)
        {
            // currentTemperature = (temperature/(3.0*KB))/nTemperature;
            // printf("%d\t%f\n", n);
            // temperature = 0.0;
            // nTemperature = 0;
            saveFrame("ABC.xyz", "a", atoms);
        }
    }
}