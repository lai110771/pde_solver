#ifndef SOLVER_H
#define SOLVER_H
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <string>
#include "Coordinate.h"

class Solver
{
    private:
        /*
        _rtol : tolerance for the error
        _maxiter : maximum number of iteration
        _T : the matrix which stores the solution
        */
        double _rtol = pow(10,-4);
        int _maxiter = pow(10,4);
        std::string _solver;
        Eigen::MatrixXd _T;
    
    public:
        /*
        GaussSeide : Call Gauss-Seidel method. it returns the solution as a matrix
        Jacobi : call Jacobi method. it returns the solution as a matrix
        */
        Solver(){}
        Solver(Coordinate& c);
        Solver(Coordinate& c, std::string method);
        Eigen::MatrixXd GaussSeidel(int Nx, int Ny, Eigen::VectorXd b);
        Eigen::MatrixXd Jacobi(int Nx, int Ny, Eigen::MatrixXd A, Eigen::VectorXd b);
        Eigen::MatrixXd get_T()
        {
            return _T;
        }
};

#endif