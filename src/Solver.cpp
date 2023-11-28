#include "Solver.h"



Solver::Solver(Coordinate& c)
{
    // User input
    std::cout << "choose method : GaussSeidel(g) or Jacobi(j)" << std::endl;
    std::cin >> _solver;
    // Choose either Gauss-Seidel method or Jacobi method
    // Then call the respect function 
    // if a wrong name is input, show error message and exit
    if(_solver == "g")
    {
        _T = GaussSeidel(c.get_m(), c.get_n(), c.get_vec());
    }
    else if(_solver == "j")
    {
        _T = Jacobi(c.get_m(), c.get_n(), c.get_mat(), c.get_vec());
    }
    else
    {
        std::cout << "unrecognized name" << std::endl;
        exit(EXIT_FAILURE);
    }
}


Solver::Solver(Coordinate& c, std::string method)
    : _solver(method)
{
    if(_solver == "g")
    {
        _T = GaussSeidel(c.get_m(), c.get_n(), c.get_vec());
    }
    else if(_solver == "j")
    {
        _T = Jacobi(c.get_m(), c.get_n(), c.get_mat(), c.get_vec());
    }
    else
    {
        std::cout << "unrecognized name" << std::endl;
        exit(EXIT_FAILURE);
    }
}


// Gauss-Seidel method
Eigen::MatrixXd Solver::GaussSeidel(int Nx, int Ny, Eigen::VectorXd b)
{
    // Initialize T matrix which stores the temperatures
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(Nx+2,Ny+2);

    /*
    invhx2 : inverse of hx = 1/(Nx+1)^2
    invhy2 : inverse of hy = 1/(Ny+1)^2
    diagCoeff : Coefficients of the diagonal elements
    xFactor : the factor for the x values
    yFactor : the factor for the y values
    currIter : current number of iteration

    we define the residual norm as : residual = ( 1/(Nx*Ny) * sum(b - sum(Ax))^2 )^(1/2)
    toSquare : the value inside of the square of residual norm
    */
    double invhx2 = pow(Nx+1,2);
    double invhy2 = pow(Ny+1,2);
    double diagCoeff = -2*(invhx2 + invhy2);
    double xFactor = invhx2/diagCoeff;
    double yFactor = invhy2/diagCoeff;
    double toSquare = 0;
    int currIter = 0;
    
    


    while(true)
    {
        /*
        Consider 5 stencil
        C(i,j)T(i,j) + C(i-1,j)T(i-1,j) + C(i+1,j)T(i+1,j) + C(i,j-1)T(i,j-1) + C(i,j+1)T(i,j+1) = b(j*Nx+i)
        update each T elements iteratively
        */
        for(auto c = 1; c < Ny+1; c++)
        {
            for(auto r = 1; r < Nx+1; r++)
            {
                T(r, c) = b(r - 1 + (c-1)*Nx)/diagCoeff - xFactor*(T(r-1, c) + T(r+1, c)) - yFactor*(T(r, c-1) + T(r, c+1));
            }
        }

        currIter = currIter + 1;
        double residual = 0;

        // Calculate the residual
        for(auto i = 1; i < Nx+1; i++)
        {
            for(auto j = 1; j < Ny+1; j++)
            {
                toSquare = b(i - 1 + (j-1)*Nx) - diagCoeff*T(i, j) - invhx2*(T(i-1, j) + T(i+1, j)) - invhy2*(T(i, j-1) + T(i, j+1));
                residual = residual + toSquare*toSquare;
            }
        }     

        // if residual is below the tolerance, end computation
        // if the iteration number exceed its maximum, then stop computation and exit
        if(residual/Nx/Ny < _rtol*_rtol)
        {
            break;
        }
        else if(currIter > _maxiter)
        {
            std::cout << residual << std::endl;
            std::cout << "Exceeded maximum iteration number" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    return T;
}


// Jacobi method
Eigen::MatrixXd Solver::Jacobi(int Nx, int Ny, Eigen::MatrixXd A, Eigen::VectorXd b)
{
    /*
    n : the length of the matrix A
    l : lower triangular matrix
    r : upper triangular matrix
    lr : represents an element of l * r
    */
    int n = Nx*Ny;
    Eigen::MatrixXd l = Eigen::MatrixXd::Zero(n,n);
    Eigen::MatrixXd r = Eigen::MatrixXd::Zero(n,n);

    double lr = 0;

    // Calculation of LU for A
    for(auto i = 0; i < n; i++)
    {
        for(auto k = 0; k < i; k++)
        {
            lr = 0;
            for(auto j = 0; j < k; j++)
            {
                lr = lr + l(i, j) * r(j, k);
            }
            l(i, k) = (A(i, k) - lr) / r(k, k);
        }
        l(i, i) = 1;        // diagonal elements in l are fixed to 1
        for(auto k = i; k < n; k++)
        {
            lr = 0;
            for(auto j = 0; j < i; j++)
            {
                lr = lr + l(i, j) * r(j, k);
            }
            r(i, k) = A(i, k) - lr;
        }
    }

    
    // Solve the linear system Ax = b with LU 
    

    // solution after forward substitution
    Eigen::VectorXd y = Eigen::VectorXd::Zero(n);
    // the sum of the known elements in each forward steps
    double er = 0;

    // Forward Substitution
    for(auto i = 0; i < n; i++)
    {
        er = 0;
        if(i == 0)
        {
            y(i) = b(i)/l(i, i);
        }
        else
        {
            for(auto j = 0; j < i; j++)
            {
                er = er + l(i, j) * y(j);
            }
            y(i) = (b(i) - er)/l(i, i);
        }
    }


    
    // final solution
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);

    // Backward Substitution
    for(auto m = 0; m < n; m++)
    {
        er = 0;
        auto i = n - m - 1;
        if(i == n-1)
        {
            x(i) = y(i)/r(i, i);
        }
        else
        {
            for(auto j = i+1; j < n; j++)
            {
                er = er + r(i, j) * x(j);
            }
            x(i) = (y(i) - er)/r(i, i);
        }
    }


    

    // eliminate the errors with iterative refinement
    // residual for Ax = b
    Eigen::VectorXd residual = b - A * x;
    // w is the solution for the residual after processing forward substitution
    // z is the solution for w after backward substitution
    Eigen::VectorXd w = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd z = Eigen::VectorXd::Zero(n);
    // count the number of iteration
    int counter = 0;


    // Iterative refinement
    while(true)
    {
        residual = b - A * x;
        counter++;

        // Solve the linear system (A-diag)z = residual with LU 
        // Forward Substitution
        er = 0; 
        for(auto i = 0; i < n; i++)
        {
            er = 0;
            if(i == 0)
            {
                w(i) = residual(i)/l(i, i);
            }
            else
            {
                for(auto j = 0; j < i; j++)
                {
                    er = er + l(i, j) * w(j);
                }
                w(i) = (residual(i) - er)/l(i, i);
            }
        }



        // Backward Substitution
        for(auto m = 0; m < n; m++)
        {
            er = 0;
            auto i = n - m - 1;
            if(i == n-1)
            {
                z(i) = w(i)/r(i, i);
            }
            else
            {
                for(auto j = i+1; j < n; j++)
                {
                    er = er + r(i, j) * z(j);
                }
                z(i) = (w(i) - er)/r(i, i);
            }
        }



        // Update the solution x = x + z
        for(auto i = 0; i < n; i++)
        {
            x(i) = x(i) + z(i);
        }
        
        // if the norm is below the tolerance, stop iterative refinement 
        // if the iteration number exceeds the maximum, exit
        if((z.norm() / x.norm()) < _rtol) 
        {
            break;
        }
        if (counter > _maxiter)
        {
            std::cout << "exceeded maximum number of iteration" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    
    // relocate the solution vector x to a matrix T with the regard of boundary values
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(Nx+2,Ny+2);
    for(auto j = 1; j < Ny+1; j++)
    {
        for(auto i = 1; i < Nx+1; i++)
        {
            T(j,i) = x((j-1)*Nx+i-1);
        }
    }

    return T;
}