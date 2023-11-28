#include "Coordinate.h"

/*
    _mat : the matrix of the linear system
    _vec : the right hand value of the linear system
*/

Coordinate::Coordinate()
{
    // User input
    std::cout << "choose coordinate : Cartesian(c) or Polar(p)" << std::endl;
    std::cin >> _type;

    // User decide the number of nodes in x and y directions
    if(_type == "c")
    {
        std::cout << "number of nodes in x direction" << std::endl;
        std::cin >> _m;
        std::cout << "number of nodes in y direction" << std::endl;
        std::cin >> _n;
        // Call functions that create the matrix and right hand side vector of the linear system
        cartesianMat();
        cartesianVec();
    }
    else if(_type == "p")
    {
        std::cout << "number of nodes in radius direction" << std::endl;
        std::cin >> _m;
        std::cout << "number of nodes in angular direction" << std::endl;
        std::cin >> _n;
        // Call functions that create the matrix and right hand side vector of the linear system
        polarMat();
        polarVec();
    }
    else
    {
        std::cout << "unrecognized name" << std::endl;
        exit(EXIT_FAILURE);
    }
}


Coordinate::Coordinate(std::string c_type, int m, int n)
    : _m(m), _n(n), _type(c_type)
{
    if(_type == "c")
    {
        // Call functions that create the matrix and right hand side vector of the linear system
        cartesianMat();
        cartesianVec();
    }
    else if(_type == "p")
    {
        // Call functions that create the matrix and right hand side vector of the linear system
        polarMat();
        polarVec();
    }
    else
    {
        std::cout << "unrecognized name" << std::endl;
        exit(EXIT_FAILURE);
    }
}


// Create the matrix of the linear system for the Cartesian coordinate
void Coordinate::cartesianMat()
{
    // Initialize the matrix
    _mat = Eigen::MatrixXd::Zero(_m*_n,_m*_n);

    /*
    Second derivative of the finite difference
    Txx = (T(i-1,j) - 2T(i,j) + T(i+1,j)) / hx      hx = i/(Nx+1)
    Tyy = (T(i,j-1) - 2T(i,j) + T(i,j+1)) / hy      hy = j/(Ny+1)
    Assume that the boundary is all zero
    */
    for(auto j = 0; j < _n; j++)
    {
        for(auto i = 0; i < _m; i++)
        {
            _mat(j*_m+i,j*_m+i) = -2*(pow(_m+1,2)+pow(_n+1,2));
            if(j != 0)
            {
                _mat(j*_m+i,(j-1)*_m+i) = pow(_n+1,2);
            }
            if(i != 0)
            {
                _mat(j*_m+i-1,j*_m+i) = pow(_m+1,2);
            }
            if(j != _n-1)
            {
                _mat(j*_m+i,(j+1)*_m+i) = pow(_n+1,2);
            }
            if(i != _m-1)
            {
                _mat(j*_m+i+1,j*_m+i) = pow(_m+1,2);
            }

        }
    }
    /*
    std::cout << "Matrix A of the linear system" << std::endl;
    std::cout << mat << std::endl;
    std::cout << std::endl;
    */
}

// Create the vector of the linear system for Cartesian coordinate
void Coordinate::cartesianVec()
{
    // Initialize the vector
    _vec = Eigen::VectorXd::Zero(_m*_n);

    // Approximate the steady state value with sin functions
    // b = -2*(pi^2)*sin(pi*x)*sin(pi*y)
    for(auto j = 0; j < _n; j++)
    {
        for(auto i = 0; i < _m; i++)
        {
            _vec(j*_m+i) = -2*pow(M_PI,2)*sin(M_PI*(i+1)/(_m+1))*sin(M_PI*(j+1)/(_n+1));
        }
    }
    /*
    std::cout << "Vector b of the linear system" << std::endl;
    std::cout << vec << std::endl;
    std::cout << std::endl;
    */
}



// Create the matrix of the linear system for the Polar coordinate
void Coordinate::polarMat()
{
    // Initialize the Matrix
    _mat = Eigen::MatrixXd::Zero(_m*_n,_m*_n);

    /*
    Second derivative of the finite difference
    Trr = (T(i-1,j) - 2T(i,j) + T(i+1,j)) / hr      hr = i/(Nr+1)
    Taa = (T(i,j-1) - 2T(i,j) + T(i,j+1)) / ha      ha = j/(Na+1)*i/(Nr+1)
    Assume that the boundary is all zero
    */
    for(auto j = 0; j < _n; j++)
    {
        for(auto i = 0; i < _m; i++)
        {
            _mat(j*_m+i,j*_m+i) = -2*(pow(_m+1,2)*(1+pow(_n+1,2)));
            if(j != 0)
            {
                _mat(j*_m+i,(j-1)*_m+i) = pow(_m+1,2)*pow(_n+1,2);
            }
            if(i != 0)
            {
                _mat(j*_m+i-1,j*_m+i) = pow(_m+1,2);
            }
            if(j != _n-1)
            {
                _mat(j*_m+i,(j+1)*_m+i) = pow(_m+1,2)*pow(_n+1,2);
            }
            if(i != _m-1)
            {
                _mat(j*_m+i+1,j*_m+i) = pow(_m+1,2);
            }
        }
    }

    /*
    std::cout << "Matrix A of the linear system" << std::endl;
    std::cout << mat << std::endl;
    std::cout << std::endl;
    */
}

// Create the vector of the linear system for Polar coordinate
void Coordinate::polarVec()
{
    // Initialize the b vector 
    _vec = Eigen::VectorXd::Zero(_m*_n);

    // Approximate the steady state value with cos functions
    // b = -2*(pi^2)*cos(pi*r/2)
    for(auto j = 0; j < _n; j++)
    {
        for(auto i = 0; i < _m; i++)
        {
            _vec(j*_m+i) = -2*pow(M_PI,2)*cos(M_PI*(i+1)/(_m+1)/2);
        }
    }
    /*
    std::cout << "Vector b of the linear system" << std::endl;
    std::cout << vec << std::endl;
    std::cout << std::endl;
    */
}