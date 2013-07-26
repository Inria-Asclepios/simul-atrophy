#ifndef ADLEM2D_HXX
#define ADLEM2D_HXX


#include <iostream>

class AdLem2D{
public:
    enum bcType{
        DIRICHLET, NEUMANN
    };

    AdLem2D(int mx, int my, double muC, double lambdaC,
            double muR, double lambdaR);

    long double dataAt(std::string dType, int x, int y);
    int getXnum() const;
    int getYnum() const;
    bcType getBcType() const;

protected:
    int mXnum, mYnum;
    AdLem2D::bcType mBc;

    //image muc, lambdac and atrophy.
    //Probably itk image type! but for now, let's play with
    //scalars:
    long double mMuGm, mMuWm, mMuCsf;
    long double mLambdaGm, mLambdaWm, mLambdaCsf;

    long double muAt(int x, int y) const;
    long double lambdaAt(int x, int y) const;
    long double aAt(int x, int y) const;


};


#endif // ADLEM2D_HXX
