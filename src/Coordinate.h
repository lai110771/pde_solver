#ifndef COORDINATE_H
#define COORDINATE_H
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <string>

class Coordinate
{
    private:
        /*
        _m : number of nodes in x or radial direction
        _n : number of nodes in y or angular direction
        _mat : the matrix of the linear system
        _vec : the right hand value of the linear system
        */

        int _m;
        int _n;
        std::string _type;
        Eigen::MatrixXd _mat;
        Eigen::VectorXd _vec;
    
    public:
        /*
        cartesianMat() : create a matrix for cartesian coordinate
        cartesianVec() : create a vector of right hand side values for cartesian cooridnate
        polarMat() : create a matrix for polar coordinate
        polarVec() : create a vector of right hand side values for polar coordinate
        */
        Coordinate();
        Coordinate(std::string c_type, int m, int n);
        void cartesianMat();
        void cartesianVec();
        void polarMat();
        void polarVec();
        int get_m()
        {
            return _m;
        }
        int get_n()
        {
            return _n;
        }
        Eigen::MatrixXd& get_mat()
        {
            return _mat;
        }
        Eigen::VectorXd& get_vec()
        {
            return _vec;
        }
};

#endif