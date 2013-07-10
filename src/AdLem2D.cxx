#include"AdLem2D.hxx"

long double AdLem2D::SetAtrophy(){return 2.;}

long double AdLem2D::muAt(int x, int y) const{
    return mTst;
}

long double AdLem2D::lambdaAt(int x, int y) const{
    return mTst;
}

long double AdLem2D::aAt(int x, int y) const{
    return mTst;
}

int AdLem2D::getXnum() const { return mXnum;}

int AdLem2D::getYnum() const { return mYnum;}

AdLem2D::bcType AdLem2D::getBcType() const { return mBc; }

long double AdLem2D::getNu() const  { return 2.; }
