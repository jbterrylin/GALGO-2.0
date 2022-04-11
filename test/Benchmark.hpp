#pragma once

#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

#include "../src/Galgo.hpp"

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
        for (size_t i = 0; i < x.size(); i++) {
            xd[i] = (double)x[i];
        }

        double obj = -sphere<double>(xd);
        return { obj };
    }
};

template <typename T>
double axisParallelHyperEllipsoid(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 0; i < particle.size(); i++) {
        sum += ((i+1) * pow(particle[i], 2.));
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
        for (size_t i = 0; i < x.size(); i++) {
            xd[i] = (double)x[i];
        }

        double obj = -axisParallelHyperEllipsoid<double>(xd);
        return { obj };
    }
};

template <typename T>
double rotatedHyperEllipsoid(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 0; i < particle.size(); i++) {
        double inner (0.);
        for (int j = 0; j <= i; j++) {
            //original:
            inner += particle[j];

            // Other papper:
            // inner += pow(particle[j],2);
        }
        //original:
        sum += pow(inner,2);

        // Other papper:    
        // sum += inner;
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
        for (size_t i = 0; i < x.size(); i++) {
            xd[i] = (double)x[i];
        }

        double obj = -rotatedHyperEllipsoid<double>(xd);
        return { obj };
    }
};

template <typename T>
double rotatedHyperEllipsoidCorrected(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 0; i < particle.size(); i++) {
        double inner (0.);
        for (int j = 0; j <= i; j++) {
            //original:
            // inner += particle[j];

            // Other papper:
            inner += pow(particle[j],2);
        }
        //original:
        // sum += pow(inner,2);

        // Other papper:    
        sum += inner;
    }
    
    return sum;
}

template <typename T>
class RotatedHyperEllipsoidObjectiveCorrected
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) {
            xd[i] = (double)x[i];
        }

        double obj = -rotatedHyperEllipsoidCorrected<double>(xd);
        return { obj };
    }
};

template <typename T>
double normalizedSchwefel(std::vector< T > particle) 
{
    double sum(0.);

    // original:
    for (int i = 0; i < particle.size(); i++) {
        sum += ( -(particle[i]) * sin(sqrt(fabs(particle[i]))) );
    }

    return sum;

    // Other papper:
    // for (int i = 0; i < particle.size(); i++) {
    //     sum += ( -(particle[i]) * sin(sqrt(fabs(particle[i]))) );
    // }

    // return sum/particle.size();
}

template <typename T>
class NormalizedSchwefelObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = -normalizedSchwefel<double>(xd);
        return { obj };
    }
};

template <typename T>
double normalizedSchwefelCorrected(std::vector< T > particle) 
{
    double sum(0.);

    // original:
    // for (int i = 0; i < particle.size(); i++) {
    //     sum += ( -(particle[i]) * sin(sqrt(fabs(particle[i]))) );
    // }

    // return sum;

    // Other papper:
    for (int i = 0; i < particle.size(); i++) {
        sum += ( -(particle[i]) * sin(sqrt(fabs(particle[i]))) );
    }

    return sum/particle.size();
}

template <typename T>
class NormalizedSchwefelObjectiveCorrected
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = -normalizedSchwefelCorrected<double>(xd);
        return { obj };
    }
};

template <typename T>
double generalizedRastrigin(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 0; i < particle.size(); i++) {
        sum += ( pow(particle[i],2) - 10 * cos(2 * PI * particle[i]) );
    }

    // original:
    return 10 * particle.size() - sum;

    // corrected
    // return 10 * particle.size() + sum;
}

template <typename T>
class GeneralizedRastriginObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = -generalizedRastrigin<double>(xd);
        return { obj };
    }
};

template <typename T>
double generalizedRastriginCorrected(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 0; i < particle.size(); i++) {
        sum += ( pow(particle[i],2) - 10 * cos(2 * PI * particle[i]) );
    }

    // original:
    // return 10 * particle.size() - sum;

    // corrected
    return 10 * particle.size() + sum;
}

template <typename T>
class GeneralizedRastriginObjectiveCorrected
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];

        double obj = -generalizedRastriginCorrected<double>(xd);
        return { obj };
    }
};

template <typename T>
double rosenbrocksValley(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 0; i < particle.size()-1; i++) {
        // original:
        sum += ( 100 - pow(particle[i+1] - pow(particle[i],2),2) + pow(1 - particle[i],2));

        // corrected:
        // sum += ( 100 * pow(particle[i+1] - pow(particle[i],2),2) + pow(1 - particle[i],2));
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
    
        double obj = -rosenbrocksValley<double>(xd);
        return { obj };
    }
};

template <typename T>
double rosenbrocksValleyCorrected(std::vector< T > particle) 
{
    double sum(0.);

    for (int i = 0; i < particle.size()-1; i++) {
        // original:
        // sum += ( 100 - pow(particle[i+1] - pow(particle[i],2),2) + pow(1 - particle[i],2));

        // corrected:
        sum += ( 100 * pow(particle[i+1] - pow(particle[i],2),2) + pow(1 - particle[i],2));
    }

    return sum;
}

template <typename T>
class RosenbrocksValleyObjectiveCorrected
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];
    
        double obj = -rosenbrocksValleyCorrected<double>(xd);
        return { obj };
    }
};

//=================================================================================================

// CEC17

double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;

template <typename T>
double cec17_entrance(std::vector< T > genes, int fun_num)
{
    int i,j,k,n,m,func_num=fun_num;
	double *f,*x;
	// FILE *fpt;
	// char FileName[30];
	m=1;
	n=genes.size();
	x=(double *)malloc(m*n*sizeof(double));
	f=(double *)malloc(sizeof(double)  *  m);
	// for (i = 0; i < 30; i++)
	// {
    // func_num=i+1;
    // sprintf(FileName, "input_data/shift_data_%d.txt", func_num);
    // fpt = fopen(FileName,"r");
    // if (fpt==NULL)
    // {
    //     printf("\n Error: Cannot open input file for reading \n");
    // }
    
    if (x==NULL)
        printf("\nError: there is insufficient memory available!\n");

    for(k=0;k<n;k++)
    {
        x[k]=genes[k];
        // fscanf(fpt,"%Lf",&x[k]);
        /*printf("%Lf\n",x[k]);*/
    }

    // fclose(fpt);

    // for (j = 0; j < n; j++)
    // {
    //     x[1*n+j]=0.0;
    //     /*printf("%Lf\n",x[1*n+j]);*/
    // }
    
    
    // for (k = 0; k < 1; k++)
    // {
    cec17_test_func(x, f, n,m,func_num);
        // for (j = 0; j < 2; j++)
        // {
        //     printf(" f%d(x[%d]) = %lf,",func_num,j+1,f[j]);
        // }
        // printf("\n");
    // }
	
	// }
    double answer = f[0];
    free(x);
    free(f);
    
    return answer;
}

template <typename T>
class ShiftedandRotatedBentCigarObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& genes)
    {
        return { -cec17_entrance(genes, 1) };
    }
};

template <typename T>
class ShiftedandRotatedZakharovObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& genes)
    {
        return { -cec17_entrance(genes, 3) };
    }
};

template <typename T>
class ShiftedandRotatedRosenbrockObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& genes)
    {
        return { -cec17_entrance(genes, 4) };
    }
};

template <typename T>
class ShiftedandRotatedRastriginObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& genes)
    {
        return { -cec17_entrance(genes, 5) };
    }
};

template <typename T>
class ShiftedandRotatedLunacekBi_RastriginObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& genes)
    {
        return { -cec17_entrance(genes, 7) };
    }
};

template <typename T>
class ShiftedandRotatedLevyObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& genes)
    {
        return { -cec17_entrance(genes, 9) };
    }
};

template <typename T>
class ShiftedandRotatedSchwefelObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& genes)
    {
        return { -cec17_entrance(genes, 10) };
    }
};

//=================================================================================================

template <typename T>
double shubert(std::vector< T > particle) 
{
    double tmp1(0.),tmp2(0.),sum(0.);

    for (int i = 1; i <= 5; i++) {
        tmp1 += i * cos((i+1) * particle[0] + i);
        tmp2 += i * cos((i+1) * particle[1] + i);
    }
    sum = tmp1 * tmp2;

    return sum;
}

template <typename T>
class ShubertObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) {
            xd[i] = (double)x[i];
        }
    
        double obj = -shubert<double>(xd);
        return { obj };
    }
};

template <typename T>
double zakharov(std::vector< T > particle) 
{
    double sum1(0.),sum2(0.),sum(0.);
    // particle = {1,2,3,4,5,6,7,8,9,10};

    for (int i=0; i<particle.size(); i++)
	{
		double xi = particle[i];
		sum1 = sum1 + pow(xi,2);
		sum2 = sum2 + 0.5*(i+1)*xi;
	}

	sum = sum1 + pow(sum2,2) + pow(sum2,4);
    // std:: cout << sum << std::endl;

    return sum;
}

template <typename T>
class ZakharovObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) {
            xd[i] = (double)x[i];
        }
    
        double obj = -zakharov<double>(xd);
        return { obj };
    }
};

template <typename T>
double dixonPrice(std::vector< T > particle) 
{
    double sum = 0;
    // particle = {1,2,3};

    double x1 = particle[0];
    double term1 = pow((x1-1),2);

    for (int i=1; i<particle.size(); i++)
    {
        double xi = particle[i];
        double xold = particle[i-1];
        double newv = (i+1) * pow((2*pow(xi,2) - xold),2);
        sum = sum + newv;
    }

    sum = term1 + sum;

    // std::cout << sum << std::endl;

    return sum;
}

template <typename T>
class DixonPriceObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++){
            xd[i] = (double)x[i];
        } 
    
        double obj = -dixonPrice<double>(xd);
        return { obj };
    }
};

//=================================================================================================

extern std::vector< galgo::Parameter<galgo::_TYPE, galgo::NBIT > > myvector;

template <typename T>
double maxProblem(const std::vector<T>& x, int oneOrZero) {
    std::string str = "";
    for (size_t i = 0; i < x.size(); i++){
        uint64_t value = (uint64_t)(galgo::Randomize<galgo::NBIT>::MAXVAL * (x[i] - myvector[i].getData()[0]) / (myvector[i].getData()[1] - myvector[i].getData()[0]));
        std::string temp = galgo::GetBinary(value);
        str += temp.substr(temp.size() - galgo::NBIT, galgo::NBIT);
    }

    int total = 0;
    for(int i=0;i< str.size();i++) {
        if(oneOrZero == 1) {
            if(str[i] == '1')
                total++;
        } else if(oneOrZero == 0) {
            if(str[i] == '0')
                total++;
        }
    }

    return { (double)total };
}

template <typename T>
class OneMaxObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        return { maxProblem(x, 1) };
    }
};

template <typename T>
class ZeroMaxObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        return { maxProblem(x, 0) };
    }
};

template <typename T>
double trapProblem(const std::vector<T>& x, int trapSize) {
    std::string str = "";
    for (size_t i = 0; i < x.size(); i++){
        uint64_t value = (uint64_t)(galgo::Randomize<galgo::NBIT>::MAXVAL * (x[i] - myvector[i].getData()[0]) / (myvector[i].getData()[1] - myvector[i].getData()[0]));
        std::string temp = galgo::GetBinary(value);
        str += temp.substr(temp.size() - galgo::NBIT, galgo::NBIT);
    }

    if (str.size() % trapSize != 0) {
        throw std::invalid_argument("bit length can't fully divide by trap size");
    }

    std::string fullzero("");
    for(int i=0; i<trapSize; i++) fullzero += "0";

    int total = 0;
    for(int i=0; i<str.size(); i+=trapSize) {
        auto block = str.substr(i, trapSize);
        if(fullzero == block) {
            total += trapSize-1;
        } else {
            for(int j=0; j<block.size(); j++) {
                if(block[j] == '1') {
                    total++;
                }
            }
        }
    }

    return { (double)total };
}

template <typename T>
class TrapThreeObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        return { trapProblem(x, 3) };
    }
};

template <typename T>
class TrapFiveObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        return { trapProblem(x, 5) };
    }
};

//=================================================================================================

template <typename T>
double rastrigrin(std::vector< T > particle) 
{
    double sum(0.);
    for (int i = 0; i < 2; i++) {
        sum = 10 + (particle[i] * particle[i]) - 10 * cos(2*PI*particle[i]);
    }

    return sum;
}

template <typename T>
class RastrigrinObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];
    
        double obj = -rastrigrin<double>(xd);
        return { obj };
    }
};

template <typename T>
double schaffersF6(std::vector< T > particle) 
{
    double sum(0.), numerator(0.), denominator(0.);

    double xysqr = pow(particle[0],2) + pow(particle[1],2);
    numerator = pow( sin( sqrt(xysqr) ), 2) - 0.5;
    denominator = pow((1+0.001*pow(xysqr,2)),2);
    sum += 0.5 + (numerator/denominator);

    return sum;
}

template <typename T>
class SchaffersF6Objective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++){
            xd[i] = (double)x[i];
        }
    
        double obj = -schaffersF6<double>(xd);
        return { obj };
    }
};

template <typename T>
double griewangks(std::vector< T > particle) 
{
    double sum(0.), s(0.), p(1.);
    for (int i=0; i<particle.size(); i++)
	{
		s += particle[i]*particle[i];
		p *= cos(particle[i]/sqrt(1.0+i));
	}
	sum = 1.0 + s/4000.0 - p;

    return sum;
}

template <typename T>
class GriewangksObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++){
            xd[i] = (double)x[i];
        }

        double obj = -griewangks<double>(xd);
        return { obj };
    }
};

template <typename T>
double hansen(std::vector< T > particle) 
{
    double sum(0.), s(0.), p(1.);
    for (int i=1; i<=5; i++)
	{
		s += i * cos( (i-1) * particle[0] + i );
		p += i * cos( (i+1) * particle[1] + i );
	}
	sum = s * p;

    return sum;
}

template <typename T>
class HansenObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++) xd[i] = (double)x[i];
    
        double obj = -hansen<double>(xd);
        return { obj };
    }
};

template <typename T>
double michalewicz(std::vector< T > particle) 
{
    double sum(0.), m(10.);
    for (int i=0; i<particle.size(); i++)
	{
		sum += sin(particle[i]) * pow( sin( (i+1) * pow(particle[i],2) / PI ), 2*m);
	}   
	sum = sum * -1;

    return sum;
}

template <typename T>
class MichalewiczObjective
{
public:
    static std::vector<double> Objective(const std::vector<T>& x)
    {
        std::vector<double> xd(x.size());
        for (size_t i = 0; i < x.size(); i++){
            xd[i] = (double)x[i];
        }
    
        double obj = -michalewicz<double>(xd);
        return { obj };
    }
};

#endif