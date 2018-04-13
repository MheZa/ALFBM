#ifndef fESection_H
#define fESection_H

#include "fENode.H"
#include "List.H"

namespace ALFBM
{
  
class fESection
{

public:

/*******************\
|    constructor    |
\*******************/

    fESection(double es,double ts,double flps,double edgs,double rho,double iner):
        eS_(es),tS_(ts),flpS_(flps),edgS_(edgs),rho_(rho),iner_(iner){}

    fESection(Foam::List<Foam::scalar> sI):
        eS_(sI[0]),tS_(sI[1]),flpS_(sI[2]),edgS_(sI[3]),rho_(sI[4]),iner_(sI[5]){}

    fESection(Foam::List<Foam::scalar> sI,bool flag):
        eS_(sI[0]),tS_(flag?sI[1]:10),flpS_(sI[2]),edgS_(sI[3]),rho_(10),iner_(10){}

    virtual ~fESection()=default;

/*******************\
|  public functions |
\*******************/

//section elastic stiffness
    const double & eS() {return eS_;}
//section torsion stiffness
    const double & tS() {return tS_;}
//section flap stiffness
    const double & flpS() {return flpS_;}
//section edge stiffness
    const double & edgS() {return edgS_;}
//section density
    const double & rho() {return rho_;}
//section rotational inertia
    const double & iner() {return iner_;}
    
private:

/*******************\
| private variables |
\*******************/

//const parameters
    //elastic stiffness
    const double eS_;
    //torsion stiffness
    const double tS_;
    //flap stiffness
    const double flpS_;
    //edge stiffness
    const double edgS_;
    //density
    const double rho_;
    //rotational inertia
    const double iner_;

/*******************\
| private functions |
\*******************/


};
        
}//end namespace ALFBM

#endif