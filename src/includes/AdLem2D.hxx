#ifndef ADLEM2D_HXX
#define ADLEM2D_HXX

class AdLem2D{
public:
    enum bcType{
        DIRICHLET, NEUMANN
    };

    AdLem2D(int mx = 5, int my = 5):
        mTst(10.), mXnum(mx), mYnum(my), mBc(AdLem2D::DIRICHLET){}

    long double SetAtrophy();
    long double muAt(int x, int y) const;
    long double lambdaAt(int x, int y) const;
    long double aAt(int x, int y) const;
    int getXnum() const;
    int getYnum() const;
    bcType getBcType() const;
    long double getNu() const;

protected:
    long double mTst; //for quickstart test purposes; should be removed later on.
    //image muc, lambdac and atrophy. Probably itk image type!
    int mXnum, mYnum;
    AdLem2D::bcType mBc;

};


#endif // ADLEM2D_HXX
