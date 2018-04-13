/*****************************************************\
|                       ALFBM                         |
|                  finiteElementNode                  |
|                       MaZhe                         |
\*****************************************************/

#ifndef fENode_H
#define fENode_H

#include "Eigen/Dense"
#include "List.H"

namespace ALFBM
{

class fENode
{
public:
	//constructor
	fENode
	(
		int i
	):
		nN_(i)
	{
		nP_=Eigen::MatrixXd::Zero(6,1);//node position
	}

	virtual ~fENode()=default;

/*******************\
|  public functions |
\*******************/

//node number
	int nN() const {return nN_;}
//node position
	const Eigen::Matrix<double,6,1> & nP() const {return nP_;}
//access to node position
	void nPSet(const Eigen::Matrix<double,6,1> nPin) {nP_ = nPin;}
	void nPSet(const Foam::point nPin) {nP_<<nPin.x(),nPin.y(),nPin.z(),0,0,0;}
	void nPAdd(const Eigen::Matrix<double,6,1> nPin) {nP_ += nPin;}
	//Eigen::Matrix<double,6,1> & nPAccess() {return nP_;}
    
private:

/*******************\
| private variables |
\*******************/
//const parameters
    //node number
	const int nN_;
//node position
	Eigen::Matrix<double,6,1> nP_;

};

}//end namespace ALFBM
#endif
