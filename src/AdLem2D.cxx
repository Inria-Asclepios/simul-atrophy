#include"AdLem2D.hxx"

AdLem2D::AdLem2D(int mx, int my, double muCsf, double lambdaCsf,
                 double muRatio, double lambdaRatio):
    mXnum(mx), mYnum(my), mBc(AdLem2D::DIRICHLET),
    mMuCsf(muCsf), mLambdaGm(lambdaCsf) {

    mMuGm = muCsf*muRatio;      mMuWm = muCsf*muRatio;
    mLambdaGm = lambdaCsf*lambdaRatio;      mLambdaWm = lambdaCsf*lambdaRatio;
}

long double AdLem2D::muAt(int x, int y) const{
    return mMuGm;
}

long double AdLem2D::lambdaAt(int x, int y) const{
    return mLambdaGm;
}

long double AdLem2D::aAt(int x, int y) const{
    if ((x > (mXnum/2. - 1)) && (x < (mXnum/2. + 1))
            && (y > (mYnum/2. -1 )) && (y < (mYnum/2. +1))) {
        return -0.2;
    } else
        return 0;
}

long double AdLem2D::dataAt(std::string dType, int x, int y)
{
    if (dType.compare("mu") == 0)
        return muAt(x,y);
    else if (dType.compare("lambda") == 0)
        return lambdaAt(x,y);
    else if (dType.compare("atrophy") == 0)
        return aAt(x,y);
    else
        std::cout<<"invalid option: "<<dType<<" : for funciton dataAt"<<std::endl;
    return 0;
}

int AdLem2D::getXnum() const { return mXnum;}

int AdLem2D::getYnum() const { return mYnum;}

AdLem2D::bcType AdLem2D::getBcType() const { return mBc; }
