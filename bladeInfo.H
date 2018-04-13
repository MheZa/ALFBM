/****************************************************************************\
This program is based on the openFOAM, and is developed by MaZhe.
The goal of this program is to build an bladeInfo class .
\****************************************************************************/

#ifndef bladeInfo_H
#define bladeInfo_H

#include "scalar.H"
#include "List.H"
#include "airfoilInfo.H"
#include <vector>
#include "Eigen/Dense"

/******************************************class declaration******************************************/

namespace Foam
{
namespace fv
{

class bladeInfo
{

public:

//Constructor
	bladeInfo
	(
		dictionary& bladeInfoPath
	):
		bladeName_(bladeInfoPath.dictName())
	{
		read(bladeInfoPath);
		interfaceCalulate();
        bladeACCalculate();
	}

//blade information APIs
	const word& bladeName() const {return bladeName_;}

//blade aerodynamic characteristic APIs
	const List<point> & bladeAEP() const {return bladeAEP_;}

	const List<word> & airfoilName() const {return airfoilName_;}

	const List<List<scalar>> & bladeAEI() const {return bladeAEI_;}

//blade finite element characteristic APIs
	const List<point> & bladeNP() const {return bladeNP_;}

	const List<List<int>> & bladeEI() const {return bladeEI_;}

    const List<List<scalar>> & bladeSI() const {return bladeSI_;}

    const List<List<scalar>> & bladeSP() const {return bladeSP_;}

    const List<scalar> & bladeAC() const {return bladeAC_;}

//blade aero finite element interface APIs
    const List<List<int>> & aeroInBeamNumber() const {return aeroInBeamNumber_;}

    const List<scalar> & aeroInBeamK() const {return aeroInBeamK_;}

private:
//blade name
	word bladeName_;

//aero infomations of blade
//blade aero element position
	List<point> bladeAEP_;

//blade airfoil type
	List<word> airfoilName_;

//blade aero element info
	List<List<scalar>> bladeAEI_;

//structure informations of blade
//blade node position
	List<point> bladeNP_;

//blade element info, node number of element
	List<List<int>> bladeEI_;

//blade element section info
    List<List<scalar>> bladeSI_;

//blade element section position
    List<List<scalar>> bladeSP_;

//blade element aerodynamic center
    List<scalar> bladeAC_;

//blade aero and finite element interface
	List<List<int>> aeroInBeamNumber_;

	List<scalar> aeroInBeamK_;

//private member functions
	bool read(dictionary& bladeDict);

	bool check();

	void interfaceCalulate();

    void bladeACCalculate();

};

}//end namespace fv
}//end namespace Foam


/******************************************************************************************************************************\
|                                                                                                                              |
|                                                    function definition                                                       |
|                                                                                                                              |
\******************************************************************************************************************************/

/******************************************private member functions******************************************/

inline bool Foam::fv::bladeInfo::read(dictionary& bladeDict)
{
    //must find profileData!
    if(bladeDict.found("bladeAEP"))
    {
        bladeDict.lookup("bladeAEP") >> bladeAEP_;
        bladeDict.lookup("airfoilName") >> airfoilName_;
        bladeDict.lookup("bladeAEI") >> bladeAEI_;
        bladeDict.lookup("bladeNP") >> bladeNP_;
        bladeDict.lookup("bladeEI") >> bladeEI_;
        bladeDict.lookup("bladeSI") >> bladeSI_;
        bladeDict.lookup("bladeSP") >> bladeSP_;
        if(check())
        {
            return true;
        }
        else
        {
            Info<<"Data check failed for "<<bladeDict.name()<<". Please correct the blade data."<<endl;
            return false;
        }
    }
    else
    {
        Info<<"No profileData for "<<bladeDict.name()<<endl;
        return false;
    }

}

inline bool Foam::fv::bladeInfo::check()
{
    if(airfoilName_.size()!=bladeAEP_.size())
    {
        Info<<"Error: error occured in airfoilName data of "<<bladeName_<<endl;
        return false;
    }
    else
    {
        if(bladeAEI_.size()!=bladeAEP_.size())
        {
            Info<<"Error: error occured in elementInfo data of "<<bladeName_<<endl;
            return false;
        }
        else
        {
            if(bladeSI_.size()!=bladeNP_.size() || bladeSI_.size()!=bladeSP_.size())
            {
                Info<<"Error: error occured in blade element info and blade element section info data of "<<bladeName_<<endl;
                return false;
            }
            else
            {
                return true;
            }
        }
    }
}

inline void Foam::fv::bladeInfo::interfaceCalulate()
{
    List<int> temp;
    temp.setSize(2,0);
    aeroInBeamNumber_.setSize(bladeAEP_.size(),temp);
    aeroInBeamK_.setSize(bladeAEP_.size(),0.0);
    int j=1;
    forAll(bladeAEP_,i)
    {
        while(j!=bladeNP_.size())
        {
            vector temp01=bladeAEP_[i]-bladeNP_[j-1];
            vector temp02=bladeAEP_[i]-bladeNP_[j];
            if((temp01&temp02)<=0.0)
            {
                aeroInBeamNumber_[i][0]=j-1;
                aeroInBeamNumber_[i][1]=j;
                aeroInBeamK_[i]=mag(temp01)/mag(bladeNP_[j]-bladeNP_[j-1]);
                break;
            }
            else
            {
                j+=1;
            }
        }
        if(j==bladeNP_.size())
        {
            Info<<"Error: Error occured in interfaceCalulate function for aero element number "<<i<<" of blade "<<bladeName_<<endl;
        }
    }
    Info<<"interface number of "<<bladeName_<<" is "<<aeroInBeamK_<<endl;
}

inline void Foam::fv::bladeInfo::bladeACCalculate()
{
    bladeAC_.setSize(bladeAEP_.size(),0.0);
    forAll(bladeAC_,i)
    {
        scalar temp=(1-aeroInBeamK_[i])*bladeSP_[aeroInBeamNumber_[i][0]][0]+aeroInBeamK_[i]*bladeSP_[aeroInBeamNumber_[i][1]][0]-0.25;
        bladeAC_[i]=temp*bladeAEI_[i][0];
    }
}

#endif

