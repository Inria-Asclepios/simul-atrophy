#ifndef ADLEM3D_HXX
#define ADLEM3D_HXX
#include<string>
/* Linear Elastic Model 3D for AD deformation. This model contains parameters and
  inputs of the following system:
  div(mu grad(v)) - grad(p) = (mu + lambda)grad(a)
  div(v)                    = -a

where,
mu and lambda are the Lame coefficients,
v is the velcoity/displacement for small time step,
and a is the atrophy
*/

class AdLem3D{
public:
    enum bcType {
        DIRICHLET, NEUMANN
    };

    AdLem3D(int mx, int my, int mz, double muCsf,
            double lambdaCsf, double muRatio, double lambdaRatio);

    int getXnum() const;
    int getYnum() const;
    int getZnum() const;
    bcType getBcType() const;

    //string should be either of "mu", "lambda" or "atrophy"
    long double dataAt(std::string dType, int x, int y, int z);


protected:
    int mXnum, mYnum, mZnum;
    AdLem3D::bcType mBc;

    //image mu, lambda and a; probably itk images!
    //parameters for Gray matter, White matter and Csf:
    long double mMuGm, mMuWm, mMuCsf;
    long double mLambdaGm, mLambdaWm, mLambdaCsf;

    long double muAt(int x, int y, int z) const;
    long double lambdaAt(int x, int y, int z) const;
    long double aAt(int x, int y, int z) const;
};

#endif // ADLEM3D_HXX
