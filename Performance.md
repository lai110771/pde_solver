# Performance Optimization


**In order to optimize the previous version of pde_solver, several perfromance tests were executed via [Quick C++ Benchmarks](https://quick-bench.com/)**


1. **If_Branching**

    Reducing branching could be a way to optimize, and thus we start the for loop index from 1 to eliminate the condition statement.

    ```
    #include <iostream>
    #include <cmath>

    static void Option1(benchmark::State& state) {
    for (auto _ : state){
        int _m = 10, _n = 10;
        double _mat[_m*_n][_m*_n];
        for(auto j = 1; j < _n-1; j++)
        {
            for(auto i = 1; i < _m-1; i++)
            {
                _mat[j*_m+i][j*_m+i] = -2*(pow(_m+1,2)+pow(_n+1,2));
                _mat[j*_m+i][(j-1)*_m+i] = pow(_n+1,2);
                _mat[j*_m+i-1][j*_m+i] = pow(_m+1,2);
                _mat[j*_m+i][(j+1)*_m+i] = pow(_n+1,2);
                _mat[j*_m+i+1][j*_m+i] = pow(_m+1,2);
            }
        }
    }
    }
    BENCHMARK(Option1);

    static void Option2 (benchmark::State& state) {
    for (auto _ : state){
        int _m = 10, _n = 10;
        double _mat[_m*_n][_m*_n];
        for(auto j = 0; j < _n; j++)
        {
            for(auto i = 0; i < _m; i++)
            {
                _mat[j*_m+i][j*_m+i] = -2*(pow(_m+1,2)+pow(_n+1,2));
                if(j != 0)
                {
                    _mat[j*_m+i][(j-1)*_m+i] = pow(_n+1,2);
                }
                if(i != 0)
                {
                    _mat[j*_m+i-1][j*_m+i] = pow(_m+1,2);
                }
                if(j != _n-1)
                {
                    _mat[j*_m+i][(j+1)*_m+i] = pow(_n+1,2);
                }
                if(i != _m-1)
                {
                    _mat[j*_m+i+1][j*_m+i] = pow(_m+1,2);
                }
            }
        }
    }
    }
    BENCHMARK(Option2);
    ```

    [Sceenshot of Benchmark](https://gitlab.lrz.de/00000000014AE8CE/pde_solver/-/blob/performance/branching.png)

    However, in this case, those two options seemed equivalent regarding the execution time.
    

2. **Calculation of Square Root**

    Noticed that calculation of square root is somehow costly, in the comparison statement, the square root calculation is replaced by the square of both sides of the statement.
    
    ```
    #include <iostream>
    #include <cmath>

    static void Option1(benchmark::State& state) {
    for (auto _ : state){
        double residual = 10^(-5);
        double _rtol = 10^(-2);
        if(residual < _rtol*_rtol){}
    }
    }
    BENCHMARK(Option1);

    static void Option2 (benchmark::State& state) {
    for (auto _ : state){
        double residual = 10^(-5);
        double _rtol = 10^(-2);
        if(sqrt(residual) < _rtol){}
    }
    }
    BENCHMARK(Option2);
    ```

    [Sceenshot of Benchmark](https://gitlab.lrz.de/00000000014AE8CE/pde_solver/-/blob/performance/sqrt.png)

    The consequence shows that it is optimized when the square root calculation is replaced.

    figures in the graph (cpu_time):
    | Option1 | Option2 |
    | ------ | ------ |
    | 0.00000402 | 22.22567308 |

3. **Operator Divide**

    Since division is somehow expensive than multiplication regarding executing, we make a trail on substitution division by multiplication.

    ```
    #include <iostream>
    #include <cmath>

    static void Option1(benchmark::State& state) {
    for (auto _ : state){
        int _m = 100, _n = 100;
        double _vec[_m*_n];
        for(auto j = 0; j < _n; j++)
        {
            for(auto i = 0; i < _m; i++)
            {
                _vec[j*_m+i] = -2*pow(M_PI,2)*cos(M_PI*(i+1)/(_m+1)*0.5);
            }
        }
    }
    }
    BENCHMARK(Option1);

    static void Option2 (benchmark::State& state) {
    for (auto _ : state){
        int _m = 100, _n = 100;
        double _vec[_m*_n];
        for(auto j = 0; j < _n; j++)
        {
            for(auto i = 0; i < _m; i++)
            {
                _vec[j*_m+i] = -2*pow(M_PI,2)*cos(M_PI*(i+1)/(_m+1)/2);
            }
        }
    }
    }
    BENCHMARK(Option2);
    ```

    [Sceenshot of Benchmark](https://gitlab.lrz.de/00000000014AE8CE/pde_solver/-/blob/performance/division.png)

    However, in this case, those two options seemed equivalent regarding the execution time.

4. **Forward Substitution**
    
    In the forward substitution algorithm, again, the conditional statement in the for loop was encountered. In this case, trying to move the inital conditon out of the for loop to avoids the branch executed in the loop.

    ```
    #include <iostream>
    #include <vector>

    static void Option1(benchmark::State& state) {
    for (auto _ : state){
        int n = 10;
        std::vector<double> y(n,0), b(n,2);
        std::vector<std::vector<double>> l(n, std::vector<double>(n,1));
        double er = 0;
        y[0] = b[0]/l[0][0];
        for(auto i = 1; i < n; i++)
        {
            er = 0;
            for(auto j = 0; j < i; j++)
            {
                er = er + l[i][j] * y[j];
            }
            y[i] = (b[i] - er)/l[i][i];
        }
    }
    }
    BENCHMARK(Option1);

    static void Option2 (benchmark::State& state) {
    for (auto _ : state){
        int n = 10;
        std::vector<double> y(n,0), b(n,2);
        std::vector<std::vector<double>> l(n, std::vector<double>(n,1));
        double er = 0;
        for(auto i = 0; i < n; i++)
        {
            er = 0;
            if(i == 0)
            {
                y[i] = b[i]/l[i][i];
            }
            else
            {
                for(auto j = 0; j < i; j++)
                {
                    er = er + l[i][j] * y[j];
                }
                y[i] = (b[i] - er)/l[i][i];
            }
        }
    }
    }
    BENCHMARK(Option2);
    ```

    [Sceenshot of Benchmark](https://gitlab.lrz.de/00000000014AE8CE/pde_solver/-/blob/performance/forward_substitution.png)

    The consequence shows that it is optimized when the branching moved out of the looping.

    figures in the graph (cpu_time):
    | Option1 | Option2 |
    | ------ | ------ |
    | 974.9234 | 1028.1174 |


5. **Optimization Level**

    When building the Makefile, setting up the `CMAKE_CXX_FLAGS` as `-O3` to optimize the compilation.

    Here are the runtimes according to different level: 
    | -o0 | -o1 | -o2 | -o3 |
    | ------ | ------ | ------ | ------ |
    | 891ms | 667ms | 1072ms | 513ms |
