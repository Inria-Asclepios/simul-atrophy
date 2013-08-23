#include "AdLem3D.hxx"
#include<iostream>


//Initialize with Dirichlet boundary condition, no other boundary condition for now.
AdLem3D::AdLem3D(int mx, int my, int mz, double muCsf, double lambdaCsf,
                 double muRatio, double lambdaRatio):
    mXnum(mx), mYnum(my), mZnum(mz), mBc(AdLem3D::DIRICHLET),
    mMuCsf(muCsf), mLambdaCsf(lambdaCsf)
{
    mMuGm = muCsf*muRatio;      mMuWm = muCsf*muRatio;
    mLambdaGm = lambdaCsf*lambdaRatio;      mLambdaWm = lambdaCsf*lambdaRatio;
}

long double AdLem3D::muAt(int x, int y, int z) const
{
    if ( (x > (mXnum/2. - 4)) && (x < (mXnum/2. + 4))
         && (y > (mYnum/2. - 4)) && (y < (mYnum/2. + 4))
         && (z > (mZnum/2. - 4)) && (z < (mZnum/2. + 4))) {
        return mMuGm;
    }
    return mMuCsf;
}

long double AdLem3D::lambdaAt(int x, int y, int z) const
{
    if ( (x > (mXnum/2. - 4)) && (x < (mXnum/2. + 4))
         && (y > (mYnum/2. - 4)) && (y < (mYnum/2. + 4))
         && (z > (mZnum/2. - 4)) && (z < (mZnum/2. + 4))) {
        return mLambdaGm;
    }
    return mLambdaCsf;
}

long double AdLem3D::aAt(int x, int y, int z) const
{
    if ( (x > (mXnum/2. - 1)) && (x < (mXnum/2. + 1))
         && (y > (mYnum/2. - 1)) && (y < (mYnum/2. + 1))
         && (z > (mZnum/2. - 1)) && (z < (mZnum/2. + 1))) {
        return -0.2;
    } else {
        if (x == 3 && y == 3 && z == 3)
            return 0.1;
        else {
            if (x==mXnum-4 && y==3 && z==3)
                return 0.1;
            else
                return 0;
        }
    }
}


long double AdLem3D::dataAt(std::string dType, int x, int y, int z)
{
    if (dType.compare("mu") == 0)
        return muAt(x,y,z);
    else if (dType.compare("lambda") == 0)
        return lambdaAt(x,y,z);
    else if (dType.compare("atrophy") == 0)
        return aAt(x,y,z);
    else
        std::cout<<"invalid option: "<<dType<<" : for funciton dataAt"<<std::endl;
    return 0;
}

int AdLem3D::getXnum() const
{
    return mXnum;
}

int AdLem3D::getYnum() const
{
    return mYnum;
}

int AdLem3D::getZnum() const
{
    return mZnum;
}

AdLem3D::bcType AdLem3D::getBcType() const { return mBc; }
