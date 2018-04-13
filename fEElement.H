/*****************************************************\
|                       ALFBM                         |
|                 finiteElementElement                |
|                       MaZhe                         |
\*****************************************************/

#ifndef fEElement_H
#define fEElement_H

#include "fENode.H"
#include "fESection.H"

namespace ALFBM
{
	
class fEElement
{

public:

/*******************\
|    constructor    |
\*******************/

	fEElement
	(
		int i,
		fENode & N0,
		fENode & N1,
        double dampR,
		Foam::List<Foam::scalar> sI1,
        Foam::List<Foam::scalar> sI2
	):
		eN_(i),
		node0_(N0),
		node1_(N1),
        dampR_(dampR),
		fES0_(sI1),
		fES1_(sI2),
		eL_(eLCalculation())
	{
        zeroInitial();
        propertyInitial();
        valueInitial();
    }

    fEElement
    (
        int i,
        fENode & N0,
        fENode & N1,
        double dampR,
        Foam::List<Foam::scalar> sI1,
        bool flag
    ):
        eN_(i),
        node0_(N0),
        node1_(N1),
        dampR_(dampR),
        fES0_(sI1,flag?true:false),
        fES1_(sI1,flag?true:false),
        eL_(eLCalculation())
    {
        zeroInitial();
        propertyInitial();
        valueInitial();
    }

	virtual ~fEElement()=default;

/*******************\
|  public functions |
\*******************/

//element number
	const int & eN() const {return eN_;}
//node0 and node1
    const fENode & node0() const {return node0_;}
    const fENode & node1() const {return node1_;}
//coordinates transform matrix
    const Eigen::Matrix<double,12,12> & eT() const {return eT_;}
//element length
	double elementLength() const {return eL_;}
//element matrixes in local coordinates
    //element stiffness matrix 
	const Eigen::Matrix<double,12,12> & eSIL() const {return eSIL_;}
    //element mass matrix
	const Eigen::Matrix<double,12,12> & eMIL() const {return eMIL_;}
    //element damp matrix
    const Eigen::Matrix<double,12,12> & eCIL() const {return eCIL_;}
    //element stress stiffening matrix
    const Eigen::Matrix<double,12,12> & eStressStiffenIL() const {return eStressStiffenIL_;}
    //element spin softening matrix
    const Eigen::Matrix<double,12,12> & eSpinSoftenIL() const {return eSpinSoftenIL_;}
    //element rotate coriolis matrix
    const Eigen::Matrix<double,12,12> & eRCIL() const {return eRCIL_;}
//element matrixes in global coordinates
    //element stiffness matrix 
    const Eigen::Matrix<double,12,12> & eSIG() const {return eSIG_;}
    //element mass matrix
    const Eigen::Matrix<double,12,12> & eMIG() const {return eMIG_;}
    //element damp matrix
    const Eigen::Matrix<double,12,12> & eCIG() const {return eCIG_;}
    //element stress stiffening matrix
    const Eigen::Matrix<double,12,12> & eStressStiffenIG() const {return eStressStiffenIG_;}
    //element spin softening matrix
    const Eigen::Matrix<double,12,12> & eSpinSoftenIG() const {return eSpinSoftenIG_;}
    //element rotate coriolis matrix
    const Eigen::Matrix<double,12,12> & eRCIG() const {return eRCIG_;}
//load vector in local coordinates
    //element gravity
    const Eigen::Matrix<double,12,1> & eGIL() const {return eGIL_;}
    //element centrifugal center force
    const Eigen::Matrix<double,12,1> & eCCIL() const {return eCCIL_;}
//load vector in global coordinates
    //element gravity
    const Eigen::Matrix<double,12,1> & eGIG() const {return eGIG_;}
    //element centrifugal center force
    const Eigen::Matrix<double,12,1> & eCCIG() const {return eCCIG_;}
//restore load vectore
    //local coordinates
    const Eigen::Matrix<double,12,1> & eRIL() const {return eRIL_;}
    //global coordinates
    const Eigen::Matrix<double,12,1> & eRIG() const {return eRIG_;}
//access to local matrixes
    //access to local element stiffness matrix
    Eigen::Matrix<double,12,12> & eSILAccess() {return eSIL_;}
    //access to local element mass matrix
    Eigen::Matrix<double,12,12> & eMILAccess() {return eMIL_;}
    //access to local element damp matrix
    Eigen::Matrix<double,12,12> & eCILAccess() {return eCIL_;}
//basic matrix construction
    void bMC(double pitch);
//stress stiffening matrix construction
    void sSMC(double pitch, const Eigen::Matrix<double,12,1> nDis);
//rotate matrix construction
    void rMC(double pitch, const Eigen::Matrix<double,3,1> omega , const Eigen::Matrix<double,3,1> center);
//restroe load construction
    void rLC(double pitch, const Eigen::Matrix<double,12,1> nDis);
	
private:

/*******************\
| private variables |
\*******************/

//const parameters
    //element number	
	const int eN_;
    //node0 and node1
    const fENode & node0_;
    const fENode & node1_;
    //element damping ratio
    const double dampR_;
    //sectional parameters for section0
    fESection fES0_;
    //sectional parameters for section1
    fESection fES1_;
//element length   
	const double eL_;
//element transform matrix
    Eigen::Matrix<double,12,12> eT_;
//element property matrix
    //stress-strain matrix
    Eigen::Matrix<double,3,3> d0_;
    //stress-strain matrix
    Eigen::Matrix<double,1,1> d1_;
    //rho
    Eigen::Matrix<double,3,3> r0_;
    //rho
    Eigen::Matrix<double,1,1> r1_;
    //nu
    Eigen::Matrix<double,3,3> nu0_;
    //nu
    Eigen::Matrix<double,1,1> nu1_;
//element matrixes in local coordinates
    //element stiffness matrix 
    Eigen::Matrix<double,12,12> eSIL_;
    //element mass matrix
    Eigen::Matrix<double,12,12> eMIL_;
    //element damp matrix
    Eigen::Matrix<double,12,12> eCIL_;
    //element stress stiffening matrix
    Eigen::Matrix<double,12,12> eStressStiffenIL_;
    //element spin softening matrix
    Eigen::Matrix<double,12,12> eSpinSoftenIL_;
    //element rotate coriolis matrix
    Eigen::Matrix<double,12,12> eRCIL_;
//element matrixes in global coordinates
    //element stiffness matrix 
    Eigen::Matrix<double,12,12> eSIG_;
    //element mass matrix
    Eigen::Matrix<double,12,12> eMIG_;
    //element damp matrix
    Eigen::Matrix<double,12,12> eCIG_;
    //element stress stiffening matrix
    Eigen::Matrix<double,12,12> eStressStiffenIG_;
    //element spin softening matrix
    Eigen::Matrix<double,12,12> eSpinSoftenIG_;
    //element rotate coriolis matrix
    Eigen::Matrix<double,12,12> eRCIG_;
//load vector in local coordinates
    //element gravity
    Eigen::Matrix<double,12,1> eGIL_;
    //element centrifugal center force
    Eigen::Matrix<double,12,1> eCCIL_;
//load vector in global coordinates
    //element gravity
    Eigen::Matrix<double,12,1> eGIG_;
    //element centrifugal center force
    Eigen::Matrix<double,12,1> eCCIG_;
//load restored
    //local coordinate
    Eigen::Matrix<double,12,1> eRIL_;
    //global coordinate
    Eigen::Matrix<double,12,1> eRIG_;

/*******************\
| private functions |
\*******************/

//element length calculation
    double eLCalculation();
//shape function N0
    Eigen::Matrix<double,3,12> n0Cal(double k);
//shape function N1
    Eigen::Matrix<double,1,12> n1Cal(double k);
//linear strain-displacement function B0
    Eigen::Matrix<double,3,12> b0Cal(double k);
//linear strain-displacement function B1
    Eigen::Matrix<double,1,12> b1Cal(double k);
//nonlinear strain-displacement function Bnl
    Eigen::Matrix<double,3,12> bnlCal(double k);
//calculation of element stress0 in local coordinates
    Eigen::Matrix<double,3,1> stress0ILCal(double k , const Eigen::Matrix<double,12,1> nDis);
//calculation of element stress0 in global coordinates
    Eigen::Matrix<double,3,1> stress0IGCal(double k , const Eigen::Matrix<double,12,1> nDis);
//calculation of element strain0 in local coordinates
    Eigen::Matrix<double,3,1> strain0ILCal(double k , const Eigen::Matrix<double,12,1> nDis);
//calculation of element strain0 in global coordinates
    Eigen::Matrix<double,3,1> strain0IGCal(double k , const Eigen::Matrix<double,12,1> nDis);
//calculation of element stress1 in local coordinates
    Eigen::Matrix<double,1,1> stress1ILCal(double k , const Eigen::Matrix<double,12,1> nDis);
//calculation of element stress1 in global coordinates
    Eigen::Matrix<double,3,1> stress1IGCal(double k , const Eigen::Matrix<double,12,1> nDis);
//calculation of element strain1 in local coordinates
    Eigen::Matrix<double,1,1> strain1ILCal(double k , const Eigen::Matrix<double,12,1> nDis);
//calculation of element strain1 in global coordinates
    Eigen::Matrix<double,3,1> strain1IGCal(double k , const Eigen::Matrix<double,12,1> nDis);
//gauss integral funtions for matrix calculation
    //N0^T*midMatrix*N0
    void gaussIntegralN0N0(Eigen::Matrix<double,12,12> & targetMatrix, Eigen::Matrix<double,3,3> & midMatrix );
    //N1^T*midMatrix*N1
    void gaussIntegralN1N1(Eigen::Matrix<double,12,12> & targetMatrix, Eigen::Matrix<double,1,1> & midMatrix );
    //B0^T*midMatrix*B0
    void gaussIntegralB0B0(Eigen::Matrix<double,12,12> & targetMatrix, Eigen::Matrix<double,3,3> & midMatrix );
    //B1^T*midMatrix*B1
    void gaussIntegralB1B1(Eigen::Matrix<double,12,12> & targetMatrix, Eigen::Matrix<double,1,1> & midMatrix );
    //N0^T*midMatrix
    void gaussIntegralN0(Eigen::Matrix<double,12,3> & targetMatrix, Eigen::Matrix<double,3,3> & midMatrix );
    //Bnl^T*midMatrix*Bnl
    void gaussIntegralBnlBnl(Eigen::Matrix<double,12,12> & targetMatrix, Eigen::Matrix<double,3,3> & midMatrix);
//initialization of element property matrix
    void propertyInitial();
//Initialization with value 
    void valueInitial();
//Initialization with zero
    void zeroInitial();
//Calculations of all the matrixes
    //calculation of element transformation matrix
    void eTCal(double pitch);
    //calculation of element stiffness matrix in local coordinates
    void eSILCal();
    //calculation of element stiffness matrix in global coordinates
    void eSIGCal();
    //calculation of element mass matrix in local coordinates
    void eMILCal();
    //calculation of element mass matrix in global coordinates
    void eMIGCal();
    //calculation of element damp matrix in local coordinates
    void eCILCal();
    //calculation of element damp matrix in global coordinates
    void eCIGCal();
    //calculation of element stress stiffening matrix in local coordinates
    void eStressStiffenILCal(Eigen::Matrix<double,12,1> nDis);
    //calculation of element stress stiffening matrix in global coordinates
    void eStressStiffenIGCal(Eigen::Matrix<double,12,1> nDis);
    //calculation of element spin softening matrix in local coordinates
    void eSpinSoftenILCal(const Eigen::Matrix<double,3,1> omega);
    //calculation of element spin softening matrix in global coordinates
    void eSpinSoftenIGCal(const Eigen::Matrix<double,3,1> omega);
    //calculation of element rotate coriolis matrix in local coordinates
    void eRCILCal(const Eigen::Matrix<double,3,1> omega);
    //calculation of element rotate coriolis matrix in global coordinates
    void eRCIGCal(const Eigen::Matrix<double,3,1> omega);
    //calculation of element gravity in local coordinates
    void eGILCal();
    //calculation of element gravity in global coordinates
    void eGIGCal();
    //calculation of element centrifugal center force in local coordinates
    void eCCILCal(const Eigen::Matrix<double,3,1> omega , const Eigen::Matrix<double,3,1> center);
    //calculation of element centrifugal center force in global coorfinates
    void eCCIGCal(const Eigen::Matrix<double,3,1> omega , const Eigen::Matrix<double,3,1> center);
    //calculation of element restore load in local coordinates
    void eRILCal(const Eigen::Matrix<double,12,1> nDis);
    //calculation of element restore load in global coordinates
    void eRIGCal(const Eigen::Matrix<double,12,1> nDis);

};
		
}//end namespace ALFBM

inline double ALFBM::fEElement::eLCalculation()
{
    Eigen::Matrix<double,6,1> temp;
    temp=node1_.nP()-node0_.nP();
    return temp.block(0,0,3,1).norm();
}

inline Eigen::Matrix<double,3,12> ALFBM::fEElement::n0Cal(double k)
{
    Eigen::Matrix<double,3,12> temp;
    temp<<1-k,0,0,0,0,0,k,0,0,0,0,0,
        0,1-3*pow(k,2)+2*pow(k,3),0,0,0,eL_*(k-2*pow(k,2)+pow(k,3)),0,3*pow(k,2)-2*pow(k,3),0,0,0,eL_*(-pow(k,2)+pow(k,3)),
        0,0,1-3*pow(k,2)+2*pow(k,3),0,eL_*(-k+2*pow(k,2)-pow(k,3)),0,0,0,3*pow(k,2)-2*pow(k,3),0,eL_*(pow(k,2)-pow(k,3)),0;
    return temp;
}

inline Eigen::Matrix<double,1,12> ALFBM::fEElement::n1Cal(double k)
{
    Eigen::Matrix<double,1,12> temp;
    temp<<0,0,0,1-k,0,0,0,0,0,k,0,0;
    return temp;
}

inline Eigen::Matrix<double,3,12> ALFBM::fEElement::b0Cal(double k)
{
    Eigen::Matrix<double,3,12> temp;
    temp<<-1/eL_,0,0,0,0,0,1/eL_,0,0,0,0,0,
        0,(12*k-6)/pow(eL_,2),0,0,0,(6*k-4)/eL_,0,(-12*k+6)/pow(eL_,2),0,0,0,(6*k-2)/eL_,
        0,0,(12*k-6)/pow(eL_,2),0,(-6*k+4)/eL_,0,0,0,(-12*k+6)/pow(eL_,2),0,(-6*k+2)/eL_,0;
    return temp;
}

inline Eigen::Matrix<double,1,12> ALFBM::fEElement::b1Cal(double k)
{
    Eigen::Matrix<double,1,12> temp;
    temp<<0,0,0,-1/eL_,0,0,0,0,0,1/eL_,0,0;
    return temp;
}

inline Eigen::Matrix<double,3,12> ALFBM::fEElement::bnlCal(double k)
{
    Eigen::Matrix<double,3,12> temp;
    temp<<0,0,0,0,0,0,0,0,0,0,0,0,
        0,(-6*k+6*pow(k,2))/eL_,0,0,0,1-4*k+3*pow(k,2),0,(6*k-6*pow(k,2))/eL_,0,0,0,-2*k+3*pow(k,2),
        0,0,(-6*k+6*pow(k,2))/eL_,0,1-4*k+3*pow(k,2),0,0,0,(6*k-6*pow(k,2))/eL_,0,-2*k+3*pow(k,2),0;
    return temp;
}

inline Eigen::Matrix<double,3,1> ALFBM::fEElement::strain0ILCal(double k,const Eigen::Matrix<double,12,1> nDis)
{
    Eigen::Matrix<double,3,1> temp;
    temp=b0Cal(k)*eT_*nDis;
    return temp;
}

inline Eigen::Matrix<double,3,1> ALFBM::fEElement::strain0IGCal(double k,const Eigen::Matrix<double,12,1> nDis)
{
    return eT_.block(0,0,3,3).transpose()*strain0IGCal(k,nDis);
}

inline Eigen::Matrix<double,3,1> ALFBM::fEElement::stress0ILCal(double k,const Eigen::Matrix<double,12,1> nDis)
{
    return d0_*strain0ILCal(k,nDis);
}

inline Eigen::Matrix<double,3,1> ALFBM::fEElement::stress0IGCal(double k,const Eigen::Matrix<double,12,1> nDis)
{
    return eT_.block(0,0,3,3).transpose()*stress0IGCal(k,nDis);
}

inline Eigen::Matrix<double,1,1> ALFBM::fEElement::strain1ILCal(double k,const Eigen::Matrix<double,12,1> nDis)
{
    Eigen::Matrix<double,1,1> temp;
    temp=b1Cal(k)*eT_*nDis;
    return temp;
}

inline Eigen::Matrix<double,3,1> ALFBM::fEElement::strain1IGCal(double k,const Eigen::Matrix<double,12,1> nDis)
{
    Eigen::Matrix<double,3,1> temp;
    temp.setZero(3,1);
    temp.block(0,0,1,1)=strain1IGCal(k,nDis);
    return eT_.block(0,0,3,3).transpose()*temp;
}

inline Eigen::Matrix<double,1,1> ALFBM::fEElement::stress1ILCal(double k,const Eigen::Matrix<double,12,1> nDis)
{
    return d1_*strain1ILCal(k,nDis);
}

inline Eigen::Matrix<double,3,1> ALFBM::fEElement::stress1IGCal(double k,const Eigen::Matrix<double,12,1> nDis)
{
    Eigen::Matrix<double,3,1> temp;
    temp.setZero(3,1);
    temp.block(0,0,1,1)=stress1IGCal(k,nDis);
    return eT_.block(0,0,3,3).transpose()*temp;
}

inline void ALFBM::fEElement::gaussIntegralN0N0(Eigen::Matrix<double,12,12> & targetMatrix, Eigen::Matrix<double,3,3> & midMatrix )
{
    targetMatrix=eL_*0.5*(
                0.1713244924*n0Cal(0.5+0.5*(-0.9324695142)).transpose()*midMatrix*n0Cal(0.5+0.5*(-0.9324695142))+
                0.3607615730*n0Cal(0.5+0.5*(-0.6612093865)).transpose()*midMatrix*n0Cal(0.5+0.5*(-0.6612093865))+
                0.4679139346*n0Cal(0.5+0.5*(-0.2386191861)).transpose()*midMatrix*n0Cal(0.5+0.5*(-0.2386191861))+
                0.1713244924*n0Cal(0.5+0.5*(0.9324695142)).transpose()*midMatrix*n0Cal(0.5+0.5*(0.9324695142))+
                0.3607615730*n0Cal(0.5+0.5*(0.6612093865)).transpose()*midMatrix*n0Cal(0.5+0.5*(0.6612093865))+
                0.4679139346*n0Cal(0.5+0.5*(0.2386191861)).transpose()*midMatrix*n0Cal(0.5+0.5*(0.2386191861))
                );
}

inline void ALFBM::fEElement::gaussIntegralN1N1(Eigen::Matrix<double,12,12> & targetMatrix, Eigen::Matrix<double,1,1> & midMatrix )
{
    targetMatrix=eL_*0.5*(
                0.1713244924*n1Cal(0.5+0.5*(-0.9324695142)).transpose()*midMatrix*n1Cal(0.5+0.5*(-0.9324695142))+
                0.3607615730*n1Cal(0.5+0.5*(-0.6612093865)).transpose()*midMatrix*n1Cal(0.5+0.5*(-0.6612093865))+
                0.4679139346*n1Cal(0.5+0.5*(-0.2386191861)).transpose()*midMatrix*n1Cal(0.5+0.5*(-0.2386191861))+
                0.1713244924*n1Cal(0.5+0.5*(0.9324695142)).transpose()*midMatrix*n1Cal(0.5+0.5*(0.9324695142))+
                0.3607615730*n1Cal(0.5+0.5*(0.6612093865)).transpose()*midMatrix*n1Cal(0.5+0.5*(0.6612093865))+
                0.4679139346*n1Cal(0.5+0.5*(0.2386191861)).transpose()*midMatrix*n1Cal(0.5+0.5*(0.2386191861))
                );
}

inline void ALFBM::fEElement::gaussIntegralB0B0(Eigen::Matrix<double,12,12> & targetMatrix, Eigen::Matrix<double,3,3> & midMatrix )
{
    targetMatrix=eL_*0.5*(
                0.1713244924*b0Cal(0.5+0.5*(-0.9324695142)).transpose()*midMatrix*b0Cal(0.5+0.5*(-0.9324695142))+
                0.3607615730*b0Cal(0.5+0.5*(-0.6612093865)).transpose()*midMatrix*b0Cal(0.5+0.5*(-0.6612093865))+
                0.4679139346*b0Cal(0.5+0.5*(-0.2386191861)).transpose()*midMatrix*b0Cal(0.5+0.5*(-0.2386191861))+
                0.1713244924*b0Cal(0.5+0.5*(0.9324695142)).transpose()*midMatrix*b0Cal(0.5+0.5*(0.9324695142))+
                0.3607615730*b0Cal(0.5+0.5*(0.6612093865)).transpose()*midMatrix*b0Cal(0.5+0.5*(0.6612093865))+
                0.4679139346*b0Cal(0.5+0.5*(0.2386191861)).transpose()*midMatrix*b0Cal(0.5+0.5*(0.2386191861))
                );
}

inline void ALFBM::fEElement::gaussIntegralB1B1(Eigen::Matrix<double,12,12> & targetMatrix, Eigen::Matrix<double,1,1> & midMatrix )
{
    targetMatrix=eL_*0.5*(
                0.1713244924*b1Cal(0.5+0.5*(-0.9324695142)).transpose()*midMatrix*b1Cal(0.5+0.5*(-0.9324695142))+
                0.3607615730*b1Cal(0.5+0.5*(-0.6612093865)).transpose()*midMatrix*b1Cal(0.5+0.5*(-0.6612093865))+
                0.4679139346*b1Cal(0.5+0.5*(-0.2386191861)).transpose()*midMatrix*b1Cal(0.5+0.5*(-0.2386191861))+
                0.1713244924*b1Cal(0.5+0.5*(0.9324695142)).transpose()*midMatrix*b1Cal(0.5+0.5*(0.9324695142))+
                0.3607615730*b1Cal(0.5+0.5*(0.6612093865)).transpose()*midMatrix*b1Cal(0.5+0.5*(0.6612093865))+
                0.4679139346*b1Cal(0.5+0.5*(0.2386191861)).transpose()*midMatrix*b1Cal(0.5+0.5*(0.2386191861))
                );
}

inline void ALFBM::fEElement::gaussIntegralN0(Eigen::Matrix<double,12,3> & targetMatrix, Eigen::Matrix<double,3,3> & midMatrix )
{
    targetMatrix=eL_*0.5*(
                0.1713244924*n0Cal(0.5+0.5*(-0.9324695142)).transpose()*midMatrix+
                0.3607615730*n0Cal(0.5+0.5*(-0.6612093865)).transpose()*midMatrix+
                0.4679139346*n0Cal(0.5+0.5*(-0.2386191861)).transpose()*midMatrix+
                0.1713244924*n0Cal(0.5+0.5*(0.9324695142)).transpose()*midMatrix+
                0.3607615730*n0Cal(0.5+0.5*(0.6612093865)).transpose()*midMatrix+
                0.4679139346*n0Cal(0.5+0.5*(0.2386191861)).transpose()*midMatrix
                );
}

inline void ALFBM::fEElement::gaussIntegralBnlBnl(Eigen::Matrix<double,12,12> & targetMatrix, Eigen::Matrix<double,3,3> & midMatrix )
{
    targetMatrix=eL_*0.5*(
                0.1713244924*bnlCal(0.5+0.5*(-0.9324695142)).transpose()*midMatrix*bnlCal(0.5+0.5*(-0.9324695142))+
                0.3607615730*bnlCal(0.5+0.5*(-0.6612093865)).transpose()*midMatrix*bnlCal(0.5+0.5*(-0.6612093865))+
                0.4679139346*bnlCal(0.5+0.5*(-0.2386191861)).transpose()*midMatrix*bnlCal(0.5+0.5*(-0.2386191861))+
                0.1713244924*bnlCal(0.5+0.5*(0.9324695142)).transpose()*midMatrix*bnlCal(0.5+0.5*(0.9324695142))+
                0.3607615730*bnlCal(0.5+0.5*(0.6612093865)).transpose()*midMatrix*bnlCal(0.5+0.5*(0.6612093865))+
                0.4679139346*bnlCal(0.5+0.5*(0.2386191861)).transpose()*midMatrix*bnlCal(0.5+0.5*(0.2386191861))
                );
}

inline void ALFBM::fEElement::propertyInitial()
{
    d0_<<0.5*(fES0_.eS()+fES1_.eS()),0,0,
        0,0.5*(fES0_.edgS()+fES1_.edgS()),0,
        0,0,0.5*(fES0_.flpS()+fES1_.flpS());
    d1_<<0.5*(fES0_.tS()+fES1_.tS());
    r0_<<0.5*(fES0_.rho()+fES1_.rho()),0,0,
        0,0.5*(fES0_.rho()+fES1_.rho()),0,
        0,0,0.5*(fES0_.rho()+fES1_.rho());
    r1_<<0.5*(fES0_.iner()+fES1_.iner());
    nu0_<<dampR_,0,0,
        0,dampR_,0,
        0,0,dampR_;
    nu1_<<dampR_;
}

inline void ALFBM::fEElement::valueInitial()
{
    eSILCal();
    eMILCal();
    eCILCal();
}

inline void ALFBM::fEElement::zeroInitial()
{
    eT_.setZero(12,12);
    eSIG_.setZero(12,12);
    eMIG_.setZero(12,12);
    eCIG_.setZero(12,12);
    eStressStiffenIL_.setZero(12,12);
    eStressStiffenIG_.setZero(12,12);
    eSpinSoftenIL_.setZero(12,12);
    eSpinSoftenIG_.setZero(12,12);
    eRCIL_.setZero(12,12);
    eRCIG_.setZero(12,12);
    eGIL_.setZero(12,1);
    eGIG_.setZero(12,1);
    eCCIL_.setZero(12,1);
    eCCIG_.setZero(12,1);
    eRIL_.setZero(12,1);
    eRIG_.setZero(12,1);
}

inline void ALFBM::fEElement::eTCal(double pitch)
{
    eT_.setZero(12,12);

    Eigen::Matrix<double,3,3> transX_;
    Eigen::Matrix<double,3,3> transY_;
    Eigen::Matrix<double,3,3> transZ_;
    Eigen::Matrix<double,3,3> transPitch_;

    Eigen::Matrix<double,3,1> Epos0;
    Epos0 = node1_.nP().block(0,0,3,1)-node0_.nP().block(0,0,3,1);

    double theta1=atan2(Epos0(1,0),Epos0(0,0));
    double theta2=atan2(Epos0(2,0),sqrt(pow(Epos0(0,0),2)+pow(Epos0(1,0),2)));

    transY_<<cos(theta2),0,-sin(-theta2),
            0,1.0,0,
            sin(-theta2),0,cos(theta2);
    transZ_<<cos(theta1),sin(theta1),0,
            -sin(theta1),cos(theta1),0,
            0,0,1.0;

    Eigen::Matrix<double,3,3> tempET;
    tempET=transY_*transZ_;

    Eigen::Matrix<double,3,1> Epos1;
    Epos1 = tempET* 0.5* (node1_.nP().block(3,0,3,1) + node0_.nP().block(3,0,3,1));

    double theta3= Epos1(0,0);

    transX_<<1.0,0,0,
            0,cos(theta3),sin(theta3),
            0,-sin(theta3),cos(theta3);

    double pit=pitch*3.14159265857/180.0;
    transPitch_<<1.0,0,0,
                0,cos(pit),sin(pit),
                0,-sin(pit),cos(pit);

    tempET=transPitch_*transX_*transY_*transZ_;

    for(int i=0;i<4;++i)
        eT_.block(i*3,i*3,3,3)=tempET;
}

inline void ALFBM::fEElement::eSILCal()
{
    Eigen::Matrix<double,12,12> K0;
    Eigen::Matrix<double,12,12> K1;

    gaussIntegralB0B0(K0,d0_);
    gaussIntegralB1B1(K1,d1_);

    eSIL_ = K0 + K1;
}

inline void ALFBM::fEElement::eSIGCal()
{
    eSIG_=eT_.transpose()*eSIL_*eT_;
}

inline void ALFBM::fEElement::eMILCal()
{   
    Eigen::Matrix<double,12,12> M0;
    Eigen::Matrix<double,12,12> M1;

    gaussIntegralN0N0(M0,r0_);
    gaussIntegralN1N1(M1,r1_);

    eMIL_ = M0 + M1;
}

inline void ALFBM::fEElement::eMIGCal()
{
    eMIG_=eT_.transpose()*eMIL_*eT_;
}

inline void ALFBM::fEElement::eCILCal()
{
    Eigen::Matrix<double,12,12> C0;
    Eigen::Matrix<double,12,12> C1;

    gaussIntegralN0N0(C0,nu0_);
    gaussIntegralN1N1(C1,nu1_);

    eCIL_ = C0 + C1;
}

inline void ALFBM::fEElement::eCIGCal()
{
    eCIG_=eT_.transpose()*eCIL_*eT_;
}

inline void ALFBM::fEElement::eStressStiffenILCal(Eigen::Matrix<double,12,1> nDis)
{
    eStressStiffenIL_=eL_*0.5*(
                0.1713244924*bnlCal(0.5+0.5*(-0.9324695142)).transpose()*stress0ILCal(0.5+0.5*(-0.9324695142),nDis)[0]*bnlCal(0.5+0.5*(-0.9324695142))+
                0.3607615730*bnlCal(0.5+0.5*(-0.6612093865)).transpose()*stress0ILCal(0.5+0.5*(-0.6612093865),nDis)[0]*bnlCal(0.5+0.5*(-0.6612093865))+
                0.4679139346*bnlCal(0.5+0.5*(-0.2386191861)).transpose()*stress0ILCal(0.5+0.5*(-0.2386191861),nDis)[0]*bnlCal(0.5+0.5*(-0.2386191861))+
                0.1713244924*bnlCal(0.5+0.5*(0.9324695142)).transpose()*stress0ILCal(0.5+0.5*(0.9324695142),nDis)[0]*bnlCal(0.5+0.5*(0.9324695142))+
                0.3607615730*bnlCal(0.5+0.5*(0.6612093865)).transpose()*stress0ILCal(0.5+0.5*(0.6612093865),nDis)[0]*bnlCal(0.5+0.5*(0.6612093865))+
                0.4679139346*bnlCal(0.5+0.5*(0.2386191861)).transpose()*stress0ILCal(0.5+0.5*(0.2386191861),nDis)[0]*bnlCal(0.5+0.5*(0.2386191861))
                );
}

inline void ALFBM::fEElement::eStressStiffenIGCal(Eigen::Matrix<double,12,1> nDis)
{
    eStressStiffenILCal(nDis);
    eStressStiffenIG_=eT_.transpose()*eStressStiffenIL_*eT_;
}

inline void ALFBM::fEElement::eSpinSoftenILCal(const Eigen::Matrix<double,3,1> omega)
{
    Eigen::Matrix<double,3,3> omegaD;
    omegaD << 0,-omega(2,0),omega(1,0),
            omega(2,0),0,-omega(0,0),
            -omega(1,0),omega(0,0),0;

    Eigen::Matrix<double,3,3> tempMid;
    tempMid=2*r0_*eT_.block(0,0,3,3)*omegaD*omegaD*eT_.block(0,0,3,3).transpose();

    gaussIntegralN0N0(eSpinSoftenIL_,tempMid);
}

inline void ALFBM::fEElement::eSpinSoftenIGCal(const Eigen::Matrix<double,3,1> omega)
{
    eSpinSoftenILCal(omega);
    eSpinSoftenIG_=eT_.transpose()*eSpinSoftenIL_*eT_;
}

inline void ALFBM::fEElement::eRCILCal(const Eigen::Matrix<double,3,1> omega)
{
    Eigen::Matrix<double,3,3> omegaD;
    omegaD << 0,-omega(2,0),omega(1,0),
            omega(2,0),0,-omega(0,0),
            -omega(1,0),omega(0,0),0;

    Eigen::Matrix<double,3,3> tempMid;
    tempMid=2*r0_*eT_.block(0,0,3,3)*omegaD*eT_.block(0,0,3,3).transpose();

    gaussIntegralN0N0(eRCIL_,tempMid);
}

inline void ALFBM::fEElement::eRCIGCal(const Eigen::Matrix<double,3,1> omega)
{
    eRCILCal(omega);
    eRCIG_=eT_.transpose()*eRCIL_*eT_;
}

inline void ALFBM::fEElement::eGIGCal()
{
    Eigen::Matrix<double,3,1> g;
    g << 0,-9.8,0;

    Eigen::Matrix<double,12,3> nRho;
    gaussIntegralN0(nRho,r0_);

    eGIG_=eT_.transpose()*nRho*eT_.block(0,0,3,3)*g;
}

inline void ALFBM::fEElement::eGILCal()
{
    eGIGCal();
    eGIL_=eT_.transpose()*eGIG_;
}

inline void ALFBM::fEElement::eCCIGCal(const Eigen::Matrix<double,3,1> omega , const Eigen::Matrix<double,3,1> center)
{
    Eigen::Matrix<double,3,3> omegaD;
    omegaD << 0,-omega(2,0),omega(1,0),
            omega(2,0),0,-omega(0,0),
            -omega(1,0),omega(0,0),0;

    Eigen::Matrix<double,12,3> rhoNOmegaOmega;
    Eigen::Matrix<double,3,3> tempMid;
    tempMid=r0_*eT_.block(0,0,3,3)*omegaD*omegaD;

    gaussIntegralN0(rhoNOmegaOmega,tempMid);

    eCCIG_=eT_.transpose()*rhoNOmegaOmega*center;
}

inline void ALFBM::fEElement::eCCILCal(const Eigen::Matrix<double,3,1> omega , const Eigen::Matrix<double,3,1> center)
{
    eCCIGCal(omega,center);
    eCCIL_=eT_.transpose()*eCCIG_;
}

inline void ALFBM::fEElement::eRILCal(const Eigen::Matrix<double,12,1> nDis)
{
    eRIL_+=eL_*0.5*(
                0.1713244924*b0Cal(0.5+0.5*(-0.9324695142)).transpose()*stress0ILCal(0.5+0.5*(-0.9324695142),nDis)+
                0.3607615730*b0Cal(0.5+0.5*(-0.6612093865)).transpose()*stress0ILCal(0.5+0.5*(-0.6612093865),nDis)+
                0.4679139346*b0Cal(0.5+0.5*(-0.2386191861)).transpose()*stress0ILCal(0.5+0.5*(-0.2386191861),nDis)+
                0.1713244924*b0Cal(0.5+0.5*(0.9324695142)).transpose()*stress0ILCal(0.5+0.5*(0.9324695142),nDis)+
                0.3607615730*b0Cal(0.5+0.5*(0.6612093865)).transpose()*stress0ILCal(0.5+0.5*(0.6612093865),nDis)+
                0.4679139346*b0Cal(0.5+0.5*(0.2386191861)).transpose()*stress0ILCal(0.5+0.5*(0.2386191861),nDis)
                )
            +
            eL_*0.5*(
                0.1713244924*b1Cal(0.5+0.5*(-0.9324695142)).transpose()*stress1ILCal(0.5+0.5*(-0.9324695142),nDis)+
                0.3607615730*b1Cal(0.5+0.5*(-0.6612093865)).transpose()*stress1ILCal(0.5+0.5*(-0.6612093865),nDis)+
                0.4679139346*b1Cal(0.5+0.5*(-0.2386191861)).transpose()*stress1ILCal(0.5+0.5*(-0.2386191861),nDis)+
                0.1713244924*b1Cal(0.5+0.5*(0.9324695142)).transpose()*stress1ILCal(0.5+0.5*(0.9324695142),nDis)+
                0.3607615730*b1Cal(0.5+0.5*(0.6612093865)).transpose()*stress1ILCal(0.5+0.5*(0.6612093865),nDis)+
                0.4679139346*b1Cal(0.5+0.5*(0.2386191861)).transpose()*stress1ILCal(0.5+0.5*(0.2386191861),nDis)
                );
}

inline void ALFBM::fEElement::eRIGCal(const Eigen::Matrix<double,12,1> nDis)
{
    eRILCal(nDis);
    eRIG_=eT_.transpose()*eRIL_;
}

void ALFBM::fEElement::bMC(double pitch)
{
    eTCal(pitch);
    eSIGCal();
    eMIGCal();
    eCIGCal();
}

void ALFBM::fEElement::sSMC(double pitch,const Eigen::Matrix<double,12,1> nDis)
{
    eTCal(pitch);
    eStressStiffenIGCal(nDis);
}

void ALFBM::fEElement::rMC(double pitch,const Eigen::Matrix<double,3,1> omega , const Eigen::Matrix<double,3,1> center)
{
    eTCal(pitch);
    eSpinSoftenIGCal(omega);
    eRCIGCal(omega);
    eGIGCal();
    eCCIGCal(omega, center);
}

void ALFBM::fEElement::rLC(double pitch,const Eigen::Matrix<double,12,1> nDis)
{
    eTCal(pitch);
    eRIGCal(nDis);
}

#endif
