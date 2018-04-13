/*****************************************************\
|                       ALFBM                         |
|                      flagBit                        |
|                       MaZhe                         |
\*****************************************************/

#ifndef flagBit_H
#define flagBit_H

#include "dictionary.H"

namespace ALFBM
{

class flagBit
{

public:

    flagBit(){}

    const double & airDensity() const {return airDensity_;}

    const double & alphaVP() const {return alphaVP_;}

    const bool & tipCorrection() const {return tipCorrection_;}

    const bool & dynamicStall() const {return dynamicStallFlagBit_;}

    const bool & deform() const {return deformFlagBit_;}

    const bool & dynamic() const {return dynamicFlagBit_;}

    const bool & nonlinear() const {return nonlinearDeformFlagBit_;}

    const bool & spinSoften() const {return spinSoftenFlagBit_;}

    const bool & stressStiffen() const {return stressStiffenFlagBit_;}

    const bool & preciseDeform() const {return preciseDeformFlagBit_;}

    const bool & cg() const {return cgFlagBit_;}

    const bool & modalSolve() const {return modalSolveFlagBit_;}

    const bool & damp() const {return dampFlagBit_;}

    const bool & debug01() const {return debugFlagBit01_;}

    const bool & debug02() const {return debugFlagBit02_;}

    const bool & debug03() const {return debugFlagBit03_;}

    const bool & debug04() const {return debugFlagBit04_;}

    bool read(Foam::dictionary& flagBitDict);

private:
    
//air density
    double airDensity_;

    double alphaVP_;

    bool tipCorrection_;

    bool dynamicStallFlagBit_;

    bool deformFlagBit_;

    bool dynamicFlagBit_;

    bool nonlinearDeformFlagBit_;

    bool spinSoftenFlagBit_;

    bool stressStiffenFlagBit_;

    bool preciseDeformFlagBit_;

    bool cgFlagBit_;

    bool modalSolveFlagBit_;

    bool dampFlagBit_;

    bool debugFlagBit01_;

    bool debugFlagBit02_;

    bool debugFlagBit03_;

    bool debugFlagBit04_;

};

}//end namespace ALFBM

/******************************************************************************************************************************\
|                                                                                                                              |
|                                                    function definition                                                       |
|                                                                                                                              |
\******************************************************************************************************************************/

bool ALFBM::flagBit::read(Foam::dictionary& flagBitDict)
{
    dynamicStallFlagBit_= flagBitDict.lookupOrDefault<Foam::Switch>("dynamicStall", false);

    tipCorrection_=flagBitDict.lookupOrDefault<Foam::Switch>("tipCorrection", false);

    deformFlagBit_=flagBitDict.lookupOrDefault<Foam::Switch>("deform", false);

    dynamicFlagBit_=flagBitDict.lookupOrDefault<Foam::Switch>("dynamic", false);

    nonlinearDeformFlagBit_=flagBitDict.lookupOrDefault<Foam::Switch>("nonlinear", false);

    spinSoftenFlagBit_=flagBitDict.lookupOrDefault<Foam::Switch>("spinSoften",false);

    stressStiffenFlagBit_=flagBitDict.lookupOrDefault<Foam::Switch>("stressStiffen",false);

    preciseDeformFlagBit_=flagBitDict.lookupOrDefault<Foam::Switch>("preciseDeform",false);

    cgFlagBit_=flagBitDict.lookupOrDefault<Foam::Switch>("cg",false);

    modalSolveFlagBit_=flagBitDict.lookupOrDefault<Foam::Switch>("modalSolve",false);

    dampFlagBit_=flagBitDict.lookupOrDefault<Foam::Switch>("damp",false);

    debugFlagBit01_=flagBitDict.lookupOrDefault<Foam::Switch>("debug01",false);

    debugFlagBit02_=flagBitDict.lookupOrDefault<Foam::Switch>("debug02",false);

    debugFlagBit03_=flagBitDict.lookupOrDefault<Foam::Switch>("debug03",false);

    debugFlagBit04_=flagBitDict.lookupOrDefault<Foam::Switch>("debug04",false);

    airDensity_ = flagBitDict.lookupOrDefault("airDensity",1.225);

    alphaVP_ = flagBitDict.lookupOrDefault("alphaVP",1.0);

    Foam::Info<<" "<<Foam::endl;
    Foam::Info<< "*********************Flag*Bits*********************" << Foam::endl;
    Foam::Info<< "* dynamicStall active " << dynamicStallFlagBit_ <<"                           *"<< Foam::endl;
    Foam::Info<< "* deform active " << deformFlagBit_ <<"                                 *" << Foam::endl;
    Foam::Info<< "* nonlinear active " << nonlinearDeformFlagBit_ << "                              *"<< Foam::endl;
    Foam::Info<< "*********************Flag*Bits*********************" << Foam::endl;
    Foam::Info<<" "<<Foam::endl;
    return true;
}

#endif