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

    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);
    MatrixXd Hj = MatrixXd::Zero(3,4);
    double dist = sqrt(pow(px,2) + pow(py,2));

    //check zero cases
    float tol = 0.00001;
    if (dist < tol) {
        dist = tol;
    }

    double dist_2 = dist*dist;
    double dist_3 = dist_2*dist;
    double cross_mult = vx*py - vy*px;
    double px_dist = px/dist;
    double py_dist = py/dist;

    Hj(0,0) = px_dist;
    Hj(0,1) = py_dist;
    Hj(1,0) = -py/dist_2;
    Hj(1,1) = px/dist_2;
    Hj(2,0) = py*cross_mult/dist_3;
    Hj(2,1) = -px*cross_mult/dist_3;
    Hj(2,2) = px_dist;
    Hj(2,3) = py_dist;

    return Hj;
}
