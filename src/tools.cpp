#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;


Tools::Tools() {}

Tools::~Tools() {}

/******************************************************************************
 * Calculate the RMSE .
 * @param estimations matrix
 * @param ground_truth matrix
 * @return  rmse
 ******************************************************************************/
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        std::cout << "Invalid estimation or ground_truth data" << std::endl;
        return rmse;
    }

    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;

    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}


/*******************************************************************************
 * Calculate Jacobian Matrix Hj
 * @param x_state : take state matrix and compute Jacobian Matrix Hj
 *******************************************************************************/
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

    MatrixXd Hj(3,4);
    //recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    //

    //accelartor variables
    double c1 = px*px + py*py;
    double c2 = sqrt(c1);
    double c3 = c1*c2;


    //check division by zero
    if (fabs(c1) < 0.00001){
        VectorXd x_stemp(4);
        // in case of division by zero, set state to small float value
        x_stemp << 0.001, 0.01, 0.001, 0.001;
        std::cout<<"Error - Division by Zero (px and py equal zero)"<< std::endl;
        return CalculateJacobian(x_stemp);
    }
    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
            -(py/c1), (px/c1), 0, 0,
            py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj;

}
