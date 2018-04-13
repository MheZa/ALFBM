/****************************************************************************\
This program is based on the openFOAM, and is developed by MaZhe.
The goal of this program is to build an airfoilInfo class .
\****************************************************************************/

#ifndef airfoilInfo_H
#define airfoilInfo_H

#include "scalar.H"
#include "List.H"
#include "mathematicalConstants.H"
#include "simpleMatrix.H"

namespace Foam
{
namespace fv
{

class airfoilInfo
{

public:
	//Constructor
	airfoilInfo
	(
		dictionary& airfoilInfoPath
	):
		airfoilName_(airfoilInfoPath.dictName()),
		Tf_(airfoilInfoPath.lookupOrDefault("Tf", 3.0)),
		Tv_(airfoilInfoPath.lookupOrDefault("Tv", 6.0)),
		Tp_(airfoilInfoPath.lookupOrDefault("Tp", 1.7)),
		Tvl_(airfoilInfoPath.lookupOrDefault("Tvl", 11.0)),
		A1_(airfoilInfoPath.lookupOrDefault("A1", 0.3)),
		A2_(airfoilInfoPath.lookupOrDefault("A2", 0.7)),
		b1_(airfoilInfoPath.lookupOrDefault("b1", 0.14)),
		b2_(airfoilInfoPath.lookupOrDefault("b2", 0.53)),
		k0_(airfoilInfoPath.lookupOrDefault("k0", 0.0)),
		exponentM_(airfoilInfoPath.lookupOrDefault("exponentM", 2.0))
	{
		read(airfoilInfoPath);
	}

	//APIs
	const word& airfoilName() const {return airfoilName_;}
	const scalar& alpha0() const {return alpha0_;}
	const scalar& alpha1() const {return alpha1_;}
	const scalar& alpha2() const {return alpha2_;}
	const scalar& alphaSS() const {return alphaSS_;}
	const scalar& Cnalpha() const {return Cnalpha_;}
	const scalar& Tf() const {return Tf_;}
	const scalar& Tv() const {return Tv_;}
	const scalar& Tp() const {return Tp_;}
	const scalar& Tvl() const {return Tvl_;}
	const scalar& A1() const {return A1_;}
	const scalar& A2() const {return A2_;}
	const scalar& b1() const {return b1_;}
	const scalar& b2() const {return b2_;}
	const scalar& S1() const {return S1_;}
	const scalar& S2() const {return S2_;}
	const scalar& k0() const {return k0_;}
	const scalar& k1() const {return k1_;}
	const scalar& k2() const {return k2_;}
	const scalar& Cd0() const {return Cd0_;}
	const scalar& Cm0() const {return Cm0_;}
	const scalar& Cn1() const {return Cn1_;}
	const scalar& exponentM() const {return exponentM_;}

	const List<List<scalar>>& profileData() const {return profileData_;}

protected:

	word airfoilName_;

	//0-lift angle of attack
	scalar alpha0_;

	//angle of attack at f=0.7 AOA<alpha0
	scalar alphaSS_;

	//angle of attack at f=0.7 AOA>alpha0
	scalar alpha1_;

	//angle of attack at f=0.7 AOA<alpha0
	scalar alpha2_;

	//Cnalpha
	scalar Cnalpha_;

	//Tf default=3.0
	scalar Tf_;

	//Tv default=6.0
	scalar Tv_;

	//Tp default=1.7
	scalar Tp_;

	//Tvl default=11
	scalar Tvl_;

	//A1 default=0.3
	scalar A1_;

	//A2 default=0.7
	scalar A2_;

	//b1 default=0.14
	scalar b1_;

	//b2 default=0.53
	scalar b2_;

	//S1
	scalar S1_;

	//S2
	scalar S2_;

	//k0
	scalar k0_;

	//k1
	scalar k1_;

	//k2
	scalar k2_;

	//Cd0
	scalar Cd0_;

	//Cm0
	scalar Cm0_;

	//Cn1
	scalar Cn1_;

	//m
	scalar exponentM_;
/*
	//Cn2
	scalar Cn2_;

*/
	//angle, lift coefficient, drag coefficient, moment coefficient
	List<List<scalar>> profileData_;

	//List of Cn
	List<scalar> Cn_;

	scalar pi=Foam::constant::mathematical::pi;

	bool read(dictionary& airfoilDict)
	{
		//must find profileData!
		if(airfoilDict.found("profileData"))
		{
			airfoilDict.lookup("profileData") >> profileData_;
			calcCn();		
			
			//find or calculate the following parameters
			if(airfoilDict.found("alpha0"))
			{
				airfoilDict.lookup("alpha0") >> alpha0_;
			}
			else
			{
				findZeroLiftAngle();
			}
			//Info << "zero lift angle checked for "<<airfoilName_<<". "<< endl;

			if(airfoilDict.found("Cd0"))
			{
				airfoilDict.lookup("Cd0") >> Cd0_;
			}
			else
			{
				calcCdZero();
			}
			//Info << "zero drag coefficient checked for "<<airfoilName_<<". "<< endl;

			if(airfoilDict.found("Cm0"))
			{
				airfoilDict.lookup("Cm0") >> Cm0_;
			}
			else
			{
				calcCmZero();
			}
			//Info << "zero moment coefficient checked for "<<airfoilName_<<". "<< endl;

			if(airfoilDict.found("alphaSS"))
			{
				airfoilDict.lookup("alphaSS") >> alphaSS_;
			}
			else
			{
				calcAlphaSS();
			}
			//Info << "static stall angle checked for "<<airfoilName_<<". "<< endl;

			if(airfoilDict.found("alpha1"))
			{
				airfoilDict.lookup("alpha1") >> alpha1_;
			}
			else
			{
				alpha1_=0.87*alphaSS_;
			}
			//Info << "01 angle checked for "<<airfoilName_<<". "<< endl;

			if(airfoilDict.found("alpha2"))
			{
				airfoilDict.lookup("alpha2") >> alpha2_;
			}
			else
			{
				alpha2_=-alpha1_;
			}
			//Info << "02 angle checked for "<<airfoilName_<<". "<< endl;

			if(airfoilDict.found("Cnalpha"))
			{
				airfoilDict.lookup("Cnalpha") >> Cnalpha_;
			}
			else
			{
				calcCnalpha();
			}
			//Info << "slope of normal coefficient checked for "<<airfoilName_<<". "<< endl;

			if(airfoilDict.found("Cn1"))
			{
				airfoilDict.lookup("Cn1") >> Cn1_;
			}
			else
			{
				calcCnOne();
			}
			//Info << "01 normal coefficient checked for "<<airfoilName_<<". "<< endl;

			if(airfoilDict.found("S1"))
			{
				airfoilDict.lookup("S1") >> S1_;
			}
			else
			{
				calcSOne();
			}
			//Info << "s1 coefficient checked for "<<airfoilName_<<". "<< endl;

			if(airfoilDict.found("S2"))
			{
				airfoilDict.lookup("S2") >> S2_;
			}
			else
			{
				calcSTwo();
			}
			//Info << "s2 coefficient checked for "<<airfoilName_<<". "<< endl;

			if(airfoilDict.found("k1") && airfoilDict.found("k2"))
			{
				airfoilDict.lookup("k1") >> k1_;
				airfoilDict.lookup("k2") >> k2_;
			}
			else
			{
				calcKOneKTwo();
			}
			//Info << "k1 and k2 coefficient checked for "<<airfoilName_<<". "<< endl;

			return true;
		}
		else
		{
			Info<<"No profileData for "<<airfoilDict.name()<<endl;
			return false;
		}

	}

	void calcCn()
	{
		Cn_.setSize(profileData_.size(),0.0);
		forAll(Cn_,i)
		{
			Cn_[i]=sin(profileData_[i][0]*pi/180)*profileData_[i][2]+cos(profileData_[i][0]*pi/180)*profileData_[i][1];
		}
	}

	void findZeroLiftAngle()
	{
		alpha0_=0.0;
		forAll(profileData_,i)
		{
			if(profileData_[i][0]<20.0 && profileData_[i][0]>-10.0)
			{
				if(i==(profileData_.size()-1))
				{
					alpha0_=0.0;
					break;
				}
				else if(profileData_[i][1] * profileData_[i+1][1]<0)
				{
					alpha0_=((0.0-profileData_[i][1])*profileData_[i+1][0]+(profileData_[i+1][1]-0.0)*profileData_[i][0])
						/(profileData_[i+1][1]-profileData_[i][1]);
					break;
				}
			}
		}
		if(alpha0_==0.0)
		{
			Info << "Warning:findZeroLiftAngle function may has error! Please check the profileData of "<<airfoilName_<<"."<< endl;
		}
	}

	void calcCdZero()
	{
		Cd0_=0.0;
		forAll(profileData_,i)
		{
			if(profileData_[i][0]<20.0 && profileData_[i][0]>-10.0)
			{
				if(i==(profileData_.size()-1))
				{
					Cd0_=0.0;
					break;
				}
				else if(((profileData_[i][0]-alpha0_) * (profileData_[i+1][0]-alpha0_))<0)
				{
					Cd0_=((alpha0_-profileData_[i][0])*profileData_[i+1][2]+(profileData_[i+1][0]-alpha0_)*profileData_[i][2])
						/(profileData_[i+1][0]-profileData_[i][0]);
					break;
				}
			}
		}
		if(Cd0_==0.0)
		{
			Info << "Warning:calcCdZero function may has error! Please check the profileData of "<<airfoilName_<<"."<< endl;
		}
	}

	void calcCmZero()
	{
		Cm0_=0.0;
		forAll(profileData_,i)
		{
			if(profileData_[i][0]<20.0 && profileData_[i][0]>-10.0)
			{
				if(i==(profileData_.size()-1))
				{
					Cm0_=0.0;
					break;
				}
				else if(((profileData_[i][0]-alpha0_) * (profileData_[i+1][0]-alpha0_))<0)
				{
					Cm0_=((alpha0_-profileData_[i][0])*profileData_[i+1][3]+(profileData_[i+1][0]-alpha0_)*profileData_[i][3])
						/(profileData_[i+1][0]-profileData_[i][0]);
					break;
				}
			}
		}
		if(Cm0_==0.0)
		{
			Info << "Warning:calcCmZero function may has error! Please check the profileData of "<<airfoilName_<<"."<< endl;
		}
	}

	void calcAlphaSS()
	{
		alphaSS_=0.0;
		forAll(profileData_,i)
		{
			if(profileData_[i][0]>alpha0_)
			{
				if(i==(profileData_.size()-1))
				{
					alphaSS_=0.0;
					break;
				}
				else if((profileData_[i+1][2]-profileData_[i][2])/(profileData_[i+1][0]-profileData_[i][0])>0.03)
				{
					alphaSS_=profileData_[i][0];
					break;
				}
			}
		}
		if(alphaSS_==0.0)
		{
			Info << "Warning:calcAlphaSS function may has error! Please check the profileData of "<<airfoilName_<<"."<< endl;
		}
	}

	void calcCnalpha()
	{
		Cnalpha_=0.0;
		simpleMatrix<scalar> Temp(2);
		Temp[0][0]=0.0;
		Temp[0][1]=0.0;
		Temp[1][0]=0.0;
		Temp[1][1]=0.0;
		Temp.source()[0]=0.0;
		Temp.source()[1]=0.0;
		forAll(profileData_,i)
		{
			if(profileData_[i][0]>alpha0_ && profileData_[i][0]<0.5*alphaSS_)
			{
				Temp[0][0]+=1;
				Temp[0][1]+=profileData_[i][0];
				Temp[1][0]+=profileData_[i][0];
				Temp[1][1]+=sqr(profileData_[i][0]);
				Temp.source()[0]+=Cn_[i];
				Temp.source()[1]+=Cn_[i]*profileData_[i][0];
			}
		}
		Cnalpha_=Temp.solve()[1];
		if(Cnalpha_==0.0)
		{
			Info << "Warning:calcCnalpha function may has error! Please check the profileData of "<<airfoilName_<<"."<< endl;
		}
	}

	void calcCnOne()
	{
		scalar f=0.7;
		Cn1_ = Cnalpha_*alpha1_*pow((1.0 + sqrt(f))/2.0, 2);
	}

	scalar calcF(scalar Cn,scalar alpha)
	{
		scalar temp=Cnalpha_*(alpha-alpha0_);
		if(temp==0)
		{
			Info << "error:error occur in function calcF! Please check the profileData of "<<airfoilName_<<"."<< endl;
		}
		scalar f=pow(2*Foam::sqrt(Cn/temp)-1,2);
		if(f>=1)
			f=1-1e-6;
		if(f<=0)
			f=1e-6;
		return f;
	}

	void calcSOne()
	{
		S1_=0.0;
		scalar a=0.0;
		scalar b=0.0;
		forAll(profileData_,i)
		{
			if(profileData_[i][0]>0.5 && profileData_[i][0]<alpha1_)
			{
				scalar x=profileData_[i][0]-alpha1_;
				scalar y=fabs((calcF(Cn_[i],profileData_[i][0])-1)/(-0.3));
				a+=x*x*Foam::log(y);
				b+=x*x*y;
			}
		}
		if(a==0.0)
		{
			S1_=1e6;
		}
		else
		{
			S1_=b/a;
		}
		if(S1_==0.0)
		{
			Info << "Warning:calcSOne function may has error! Please check the profileData of "<<airfoilName_<<"."<< endl;
		}
	}

	void calcSTwo()
	{
		S2_=0.0;
		scalar a=0.0;
		scalar b=0.0;
		forAll(profileData_,i)
		{
			if(profileData_[i][0]>alpha1_ && profileData_[i][0]<25.0)
			{
				scalar x=alpha1_ - profileData_[i][0];
				scalar y=Foam::mag((calcF(Cn_[i],profileData_[i][0])-0.04)/0.66);
				a+=x*x*Foam::log(y);
				b+=x*x*y;
			}
		}
		if(a==0.0)
		{
			S2_=1e6;
		}
		else
		{
			S2_=b/a;
		}
		if(S2_==0.0)
		{
			Info << "Warning:calcSTwo function may has error! Please check the profileData of "<<airfoilName_<<"."<< endl;
		}
	}

	void calcKOneKTwo()
	{
		k1_=0.0;
		k2_=0.0;
		pi=Foam::constant::mathematical::pi;
		simpleMatrix<scalar> Temp(2);
		Temp[0][0]=0.0;
		Temp[0][1]=0.0;
		Temp[1][0]=0.0;
		Temp[1][1]=0.0;
		Temp.source()[0]=0.0;
		Temp.source()[1]=0.0;
		forAll(profileData_,i)
		{
			if(profileData_[i][0]>0.5 && profileData_[i][0]<25.0)
			{
				scalar f=calcF(Cn_[i],profileData_[i][0]);
				Temp[0][0]+=magSqr(1.0 - f);
				Temp[0][1]+=sin(pi*pow(f,exponentM_));
				Temp[1][0]+=sin(pi*pow(f,exponentM_))*(1.0 - f);
				Temp[1][1]+=sqr(sin(pi*pow(f,exponentM_)));
				Temp.source()[0]+=(1.0 - f)*(profileData_[i][3]/Cn_[i]-k0_);
				Temp.source()[1]+=sin(pi*pow(f,exponentM_))*(profileData_[i][3]/Cn_[i]-k0_);
			}
		}
		List<scalar> sol=Temp.solve();
		k1_=sol[0];
		k2_=sol[1];
		if(k1_==0.0 || k2_==0.0)
		{
			Info << "Warning:calcKOneKTwo function may has error! Please check the profileData of "<<airfoilName_<<"."<< endl;
		}
	}

};

}//end namespace fv
}//end namespace Foam

#endif

