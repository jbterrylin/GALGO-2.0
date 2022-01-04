#pragma once

#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

// Ring Crossover
template <typename T>
double sphere(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 0; i < particle.size(); i++) {
        sum += pow(particle[i], 2.);
    }
    return sum;
}

template <typename T>
class SphereObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = sphere<double>(xd);
        return { obj };
    }
};

template <typename T>
double axisParallelHyperEllipsoid(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 1; i <= particle.size(); i++) {
        sum += (i * pow(particle[i], 2.));
    }
    return sum;
}

template <typename T>
class AxisParallelHyperEllipsoidObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = axisParallelHyperEllipsoid<double>(xd);
        return { obj };
    }
};

template <typename T>
double rotatedHyperEllipsoid(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 1; i <= particle.size(); i++) {
        double inner (0.);
        for (int j = 1; j <= i; j++)
            inner += pow(particle[j],2);
        sum += inner;
    }

    return sum;
}

template <typename T>
class RotatedHyperEllipsoidObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = rotatedHyperEllipsoid<double>(xd);
        return { obj };
    }
};

template <typename T>
double normalizedSchwefel(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 1; i <= particle.size(); i++) {
        sum += ( -(particle[i]) * sin(sqrt(abs(particle[i]))) );
    }

    return sum;
}

template <typename T>
class NormalizedSchwefelObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = normalizedSchwefel<double>(xd);
        return { obj };
    }
};

template <typename T>
double generalizedRastrigin(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 1; i <= particle.size(); i++) {
        sum += ( pow(particle[i],2) - 10 * cos(2 * PI * particle[i]) );
    }

    return sum;
}

template <typename T>
class GeneralizedRastriginObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = generalizedRastrigin<double>(xd);
        return { obj };
    }
};

template <typename T>
double rosenbrocksValley(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 1; i <= particle.size()-1; i++) {
        sum += ( 100 * pow(particle[i+1] - pow(particle[i],2),2) + pow(particle[i] - 1,2));
    }

    return sum;
}

template <typename T>
class RosenbrocksValleyObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];
    
        double obj = rosenbrocksValley<double>(xd);
        return { obj };
    }
};

double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;

template <typename T>
class ShiftedandRotatedBentCigarObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& genes)
    {
        int i,j,k,n,m,func_num=1;
        double *f,*x;
        FILE *fpt;
        char FileName[30];
        m=1;
        n=genes.size();
        x=(double *)malloc(m*n*sizeof(double));
        f=(double *)malloc(sizeof(double)  *  m);

        sprintf(FileName, "input_data/shift_data_%d.txt", func_num);
        fpt = fopen(FileName,"r");
        if (fpt==NULL)
        {
            printf("\n Error: Cannot open input file for reading \n");
        }

        if (x==NULL)
            printf("\nError: there is insufficient memory available!\n");

        for(k=0;k<n;k++)
        {
                fscanf(fpt,"%lf",&x[k]);
                // std::cout << "x[k]" << x[k] << std::endl;
        }

        fclose(fpt);

        for (j = 0; j < n; j++)
        {
            x[1*n+j]=0.0;
        }

        z=(double *)malloc(sizeof(double)  *  genes.size());
        for (size_t i = 0; i < genes.size(); i++) z[i] = (double)genes[i];
        // for (size_t i = 0; i < genes.size(); i++) std::cout << "z[i]" << z[i] << std::endl;

        for (k = 0; k < 1; k++)
        {
            f[0] = 0;
            // for (j = 0; j < 2; j++)
            // {
            //     printf(" f%d(x[%d]) = %lf,",func_num,j+1,f[j]);
            // }
            // printf("\n");
            cec17_test_func(x, f, n,m,func_num);
            // for (j = 0; j < 2; j++)
            // {
            //     printf(" f%d(x[%d]) = %lf,",func_num,j+1,f[j]);
            // }
            // printf("\n");
        }

        // for (size_t i = 0; i < genes.size(); i++) z[i] = (double)x[i];
        // double obj = bentCigar<double>(xd);
        return { f[0] };
    }
};

template <typename T>
class ShiftedandRotatedRosenbrockObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& genes)
    {
        int i,j,k,n,m,func_num=3;
        double *f,*x;
        FILE *fpt;
        char FileName[30];
        m=1;
        n=genes.size();
        x=(double *)malloc(m*n*sizeof(double));
        f=(double *)malloc(sizeof(double)  *  m);

        sprintf(FileName, "input_data/shift_data_%d.txt", func_num);
        fpt = fopen(FileName,"r");
        if (fpt==NULL)
        {
            printf("\n Error: Cannot open input file for reading \n");
        }

        if (x==NULL)
            printf("\nError: there is insufficient memory available!\n");

        for(k=0;k<n;k++)
        {
                fscanf(fpt,"%lf",&x[k]);
                // std::cout << "x[k]" << x[k] << std::endl;
        }

        fclose(fpt);

        for (j = 0; j < n; j++)
        {
            x[1*n+j]=0.0;
        }

        z=(double *)malloc(sizeof(double)  *  genes.size());
        for (size_t i = 0; i < genes.size(); i++) z[i] = (double)genes[i];
        // for (size_t i = 0; i < genes.size(); i++) std::cout << "z[i]" << z[i] << std::endl;

        for (k = 0; k < 1; k++)
        {
            f[0] = 0;
            // for (j = 0; j < 2; j++)
            // {
            //     printf(" f%d(x[%d]) = %lf,",func_num,j+1,f[j]);
            // }
            // printf("\n");
            cec17_test_func(x, f, n,m,func_num);
            // for (j = 0; j < 2; j++)
            // {
            //     printf(" f%d(x[%d]) = %lf,",func_num,j+1,f[j]);
            // }
            // printf("\n");
        }

        // for (size_t i = 0; i < genes.size(); i++) z[i] = (double)x[i];
        // double obj = bentCigar<double>(xd);
        return { f[0] };
    }
};

#endif