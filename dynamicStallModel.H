/****************************************************************************\
This program is based on the openFOAM, and is developed by MaZhe.
The goal of this program is to build an dynamicStallModel class .
\****************************************************************************/

#ifndef dynamicStallModel_H
#define dynamicStallModel_H

#include "airfoilInfo.H"
#include "scalar.H"
#include <vector>


namespace Foam
{
namespace fv
{

class dynamicStallModel
{

public:

	// Constructors
    dynamicStallModel
    (
        const Time& time,
        airfoilInfo& airfoil
    ):
        time_(time),
        airfoilInfo_(airfoil),
        X_(0.0),
	    XPrev_(0.0),
	    Y_(0.0),
	    YPrev_(0.0),
	    alpha_(0.0),
	    deltaAlpha_(0.0),
	    deltaAlphaPrev_(0.0),
	    a_(1e12),
	    timePrev_(0.0),
	    D_(0.0),
	    DPrev_(0.0),
	    DP_(0.0),
	    DPPrev_(0.0),
	    CNP_(0.0),
	    CNPPrev_(0.0),
	    fPrime_(1.0),
	    fPrimePrev_(1.0),
	    DF_(0.0),
	    DFPrev_(0.0),
	    CV_(0.0),
	    CVPrev_(0.0),
	    CNV_(0.0),
	    CNVPrev_(0.0),
	    eta_(0.95),
	    stalledPrev_(false),
	    tau_(0.0),
	    tauPrev_(0.0),
	    nNewTimes_(0),
	    fCrit_(0.7),
	    CM_(0.0),
	    Re_(0.0)
    {
    }

    //- Correct lift, drag, and moment coefficients
    void correct
    (
        scalar magU,
        scalar alphaDeg,
        scalar& cl,
        scalar& cd,
        scalar& cm
    );
   
    // Parallel running

        //- Reduce to set data equal on all processors
        void reduceParallel(bool inMesh);


private:

    //- Runtime reference
    const Time& time_;

    //airfoilData
    airfoilInfo& airfoilInfo_;

    //- Chord length in meters
    scalar c_;

    //- Normal force coefficient from circulatory effects
    scalar CNC_;

    //- Equivalent angle of attack in radians
    scalar alphaEquiv_;

    //- Deficiency function for equivalent angle of attack
    scalar X_;

    //- Previous value of X
    scalar XPrev_;

    //- Deficiency function for equivalent angle of attack
    scalar Y_;

    //- Previous value of Y
    scalar YPrev_;

    //- Mach number
    scalar M_;

    //- Angle of attack (radians)
    scalar alpha_;

    //- Previous angle of attack (radians)
    scalar alphaPrev_;

    //- Change in angle of attack (radians)
    scalar deltaAlpha_;

    //- Previous value of deltaAlpha
    scalar deltaAlphaPrev_;

    //- Time step in seconds
    scalar deltaT_;

    //- Nondimensional time step
    scalar deltaS_;

    //- Normal force coefficient from impulsive effects
    scalar CNI_;

    //- Time constant for impulsive loading = c/a
    scalar TI_;

    //- Speed of sound (m/s)
    scalar a_;

    //- Previous time value
    scalar timePrev_;

    //- Deficiency function for the impulsive force response
    scalar D_;

    //- Previous value of D
    scalar DPrev_;

    //- Deficiency function for delayed pressure response
    scalar DP_;

    //- Previous value of DP
    scalar DPPrev_;

    //- Total normal force coefficient under the unsteady attached flow
    scalar CNP_;

    //- Previous value of CNP
    scalar CNPPrev_;

    //- Lagged normal force coefficient
    scalar CNPrime_;

    //- Delayed separation point
    scalar fPrime_;

    //- Previous value of fPrime
    scalar fPrimePrev_;

    //- Second order delayed separation point
    scalar fDoublePrime_;

    //- Delayed angle of attack (radians)
    scalar alphaPrime_;

    //- Deficiency function for the delayed separation point
    scalar DF_;

    //- Previous value of DF
    scalar DFPrev_;

    //- Normal force coefficient for the delayed separation point
    scalar CNF_;

    //- Vortex lift
    scalar CV_;

    //- Previous value of CV
    scalar CVPrev_;

    //- Vortex lift contribution to normal force coefficient
    scalar CNV_;

    //- Previous value of CNV
    scalar CNVPrev_;

    //- Total normal force coefficient
    scalar CN_;

    //- Tangential force coefficient
    scalar CT_;

    //- Tangential force efficiency factor
    scalar eta_;

    //- Switch to indicate stall
    bool stalled_;

    //- Previous stall state
    bool stalledPrev_;

    //- Nondimensional vortex tracking time
    scalar tau_;

    //- Previous value of tau
    scalar tauPrev_;

    //- Number of new times evaluated -- does not increment for subcycles
    label nNewTimes_;

    //- Critical trailing edge separation point
    scalar fCrit_;

    //- Magnitude of the relative velocity
    scalar magU_;

    //- Total moment coefficient
    scalar CM_;

    //- Reynolds number
    scalar Re_;

    //- Calculate the equivalent angle of attack
    void calcAlphaEquiv();

    //- Calculate unsteady quantities
    void calcUnsteady();

    //- Calculate separated flow quantities
    void calcSeparated();

    //- Update previous time step values
    void update();

};

}//end namespace fv
}//end namespace Foam

inline void Foam::fv::dynamicStallModel::calcAlphaEquiv()
{
    scalar beta = 1.0 - M_*M_;
    X_ = XPrev_*exp(-airfoilInfo_.b1()*beta*deltaS_)
       + airfoilInfo_.A1()*deltaAlpha_*exp(-airfoilInfo_.b1()*beta*deltaS_/2.0);
    Y_ = YPrev_*exp(-airfoilInfo_.b2()*beta*deltaS_)
       + airfoilInfo_.A2()*deltaAlpha_*exp(-airfoilInfo_.b2()*beta*deltaS_/2.0);
    alphaEquiv_ = alpha_ - X_ - Y_;
}

void Foam::fv::dynamicStallModel::calcUnsteady()
{
    // Calculate the circulatory normal force coefficient
    CNC_ = airfoilInfo_.Cnalpha()*alphaEquiv_;

    // Calculate the impulsive normal force coefficient
    scalar pi = Foam::constant::mathematical::pi;
    scalar kAlpha = 0.75
                  / (1.0 - M_ + pi*(1.0 - M_*M_)*M_*M_*(airfoilInfo_.A1()*airfoilInfo_.b1()+ airfoilInfo_.A2()*airfoilInfo_.b2()));
    TI_ = c_/a_;
    D_ = DPrev_*exp(-deltaT_/(kAlpha*TI_))
       + ((deltaAlpha_ - deltaAlphaPrev_)/deltaT_)
       *exp(-deltaT_/(2.0*kAlpha*TI_));
    CNI_ = 4.0*kAlpha*TI_/M_*(deltaAlpha_/deltaT_ - D_);

    // Calculate total normal force coefficient
    CNP_ = CNC_ + CNI_;

    // Apply first-order lag to normal force coefficient
    DP_ = DPPrev_*exp(-deltaS_/airfoilInfo_.Tp())
        + (CNP_ - CNPPrev_)*exp(-deltaS_/(2.0*airfoilInfo_.Tp()));
    CNPrime_ = CNP_ - DP_;

    // Calculate lagged angle of attack
    alphaPrime_ = CNPrime_/airfoilInfo_.Cnalpha();

    // Set stalled switch
    stalled_ = (mag(CNPrime_) > airfoilInfo_.Cn1());

}

void Foam::fv::dynamicStallModel::calcSeparated()
{
    // Calculate trailing-edge separation point
    if (mag(alphaPrime_) < airfoilInfo_.alpha1())
    {
        fPrime_ = 1.0 - 0.3*exp((mag(alphaPrime_) - airfoilInfo_.alpha1())/airfoilInfo_.S1());
    }
    else
    {
        fPrime_ = 0.04 + 0.66*exp((airfoilInfo_.alpha1() - mag(alphaPrime_))/airfoilInfo_.S2());
    }

    // Modify Tf time constant if necessary
    scalar Tf = airfoilInfo_.Tf();
    if (tau_ > 0 and tau_ <= airfoilInfo_.Tvl())
    {
        Tf = 0.5*airfoilInfo_.Tf();
    }
    else if (tau_ > airfoilInfo_.Tvl() and tau_ <= 2.0*airfoilInfo_.Tvl())
    {
        Tf = 4.0*airfoilInfo_.Tf();
    }
    if (mag(alpha_) < mag(alphaPrev_) and mag(CNPrime_) < airfoilInfo_.Cn1())
    {
        Tf = 0.5*airfoilInfo_.Tf();
    }

    // Calculate dynamic separation point
    DF_ = DFPrev_*exp(-deltaS_/Tf)
        + (fPrime_ - fPrimePrev_)*exp(-deltaS_/(2.0*Tf));
    fDoublePrime_ = fPrime_ - DF_;
    if (fDoublePrime_ < 0)
    {
        fDoublePrime_ = 0.0;
    }
    else if (fDoublePrime_ > 1)
    {
        fDoublePrime_ = 1.0;
    }

    // Calculate normal force coefficient including dynamic separation point
    CNF_ = airfoilInfo_.Cnalpha()*alphaEquiv_*pow(((1.0 + sqrt(fDoublePrime_))/2.0), 2.0)
         + CNI_;

    // Calculate tangential force coefficient
    if (fDoublePrime_ < fCrit_)
    {
        CT_ = eta_*airfoilInfo_.Cnalpha()*alphaEquiv_*alphaEquiv_*pow(fDoublePrime_, 1.5);
    }
    else
    {
        CT_ = eta_*airfoilInfo_.Cnalpha()*alphaEquiv_*alphaEquiv_*sqrt(fDoublePrime_);
    }

    // Compute vortex shedding process if stalled
    // Evaluate vortex tracking time
    if (not stalledPrev_)
    {
        tau_ = 0.0;
    }
    else
    {
        if (tau_ == tauPrev_)
        {
            tau_ = tauPrev_ + deltaS_;
        }
    }

    // Calculate Strouhal number time constant and set tau to zero to
    // allow multiple vortex shedding
    scalar Tst = 2.0*(1.0 - fDoublePrime_)/0.19;
    if (tau_ > (airfoilInfo_.Tvl() + Tst))
    {
        tau_ = 0.0;
    }

    // Evaluate vortex lift contributions, which are only increasing if angle
    // of attack increased in magnitude beyond a threshold
    scalar Tv = airfoilInfo_.Tv();
    if(mag(alpha_) > mag(alphaPrev_))
    {
        if (tau_ < airfoilInfo_.Tvl())
        {
            // Halve Tv if dAlpha/dt changes sign
            if (sign(deltaAlpha_) != sign(deltaAlphaPrev_))
            {
                Tv = 0.5*airfoilInfo_.Tv();
            }
            scalar KN = magSqr((1.0 + sqrt(fDoublePrime_)))/4.0;
            CV_ = CNC_*(1.0 - KN);
            CNV_ = CNVPrev_*exp(-deltaS_/Tv)
                 + (CV_ - CVPrev_)*exp(-deltaS_/(2.0*Tv));
        }
        else
        {
            Tv = 0.5*airfoilInfo_.Tv();
            CV_ = 0.0;
            CNV_ = CNVPrev_*exp(-deltaS_/Tv);
        }
    }
    else
    {
        CNV_ = 0.0;
    }
    

    // Total normal force coefficient is the combination of that from
    // circulatory effects, impulsive effects, dynamic separation, and vortex
    // lift
    CN_ = CNF_ + CNV_;

    // Calculate moment coefficient
    scalar pi = Foam::constant::mathematical::pi;
    scalar cmf = (airfoilInfo_.k0() + airfoilInfo_.k1()*(1 - fDoublePrime_)
               + airfoilInfo_.k2()*sin(pi*Foam::pow(fDoublePrime_, airfoilInfo_.exponentM())))*CNC_
               + airfoilInfo_.Cm0();
    scalar cpv = 0.20*(1 - cos(pi*tau_/airfoilInfo_.Tvl()));
    scalar cmv = -cpv*CNV_;
    CM_ = cmf + cmv;
}


void Foam::fv::dynamicStallModel::update()
{
    timePrev_ = time_.value();
    alphaPrev_ = alpha_;
    XPrev_ = X_;
    YPrev_ = Y_;
    deltaAlphaPrev_ = deltaAlpha_;
    DPrev_ = D_;
    DPPrev_ = DP_;
    CNPPrev_ = CNP_;
    DFPrev_ = DF_;
    fPrimePrev_ = fPrime_;
    CVPrev_ = CV_;
    CNVPrev_ = CNV_;
    stalledPrev_ = stalled_;
    tauPrev_ = tau_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fv::dynamicStallModel::correct
(
    scalar magU,
    scalar alphaDeg,
    scalar& cl,
    scalar& cd,
    scalar& cm
)
{
    scalar pi = Foam::constant::mathematical::pi;
    scalar time = time_.value();
    deltaT_ = time_.deltaT().value();

    // Update previous values if time has changed
    if (time != timePrev_)
    {
        nNewTimes_++;
        if (nNewTimes_ > 1)
        {
            update();
        }
    }

    if (nNewTimes_ <= 1)
    {
        alpha_ = alphaDeg/180.0*pi;
        alphaPrev_ = alpha_;
    }

    magU_ = magU;
    alpha_ = alphaDeg/180.0*pi;
    M_ = magU/a_;
    deltaAlpha_ = alpha_ - alphaPrev_;
    deltaS_ = 2*magU*deltaT_/c_;

    calcAlphaEquiv();
    calcUnsteady();
    calcSeparated();

    // Modify coefficients
    cl = CN_*cos(alpha_) + CT_*sin(alpha_);
    cd = CN_*sin(alpha_) - CT_*cos(alpha_) + airfoilInfo_.Cd0();
    cm = CM_;
}

#endif