# PDE SOLVER

optional project in Advanced Programming

**2D steady state heat equation on C++**

Determines the temperature distribution of a board in steady state. 

- **Physical model** : 2 dimension board has temperature distribution profile represented with sin curves.
    `Txx + Tyy = -2pi^2 * sin(pi * x) * sin(pi * y)`
The boundaries are set to 0.


- **Analytical Solution** : `T(x,y) = sin(pi * x) * sin(pi * y)`

- **Simulation Concept**: It allows different coordinates and enables users to compute with different methods. 

- **How to Execute the Code**: It accepts both cartesian coordinate and polar coordinate. Users can choose either. For the cartesian coordinate, user must input the nodes number of x and y direction. In the same way, user must input the nodes number of radial and angular direction for polar coordinate. Then, user can choose either Gauss-Seidel method or Jacobi method. Both methods returns the solution as a matrix. user can output the results in csv file.

- **Example**:In the example, when the code executed, using the configuration file(c) to do the simulation. And then take Cartesian(c) as the coordinate system and assign 10 nodes in both x and y directions. Last, choose Gauss-Seidel(g) to solver the configuration. In the end of the process, user could save the consequence as a csv. file to visualize the outcome.

In this case, the consequence shown as below.



|  |  |  |  |  |  |  |  |  |  |  |  |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 0 | 0 |0 |  0 | 0 | 0 |0 | 0 | 0 | 0| 0 |0|
| 0 | 0.079914 |0.153354 |  0.21437 | 0.258019 | 0.280765  |0.280765  | 0.25802 | 0.214371 | 0.153354| 0.0799143 |0|
| 0 | 0.153354 |0.294284  |0.411373|  0.495136|  0.538785|  0.538785|  0.495136|  0.411374|  0.294285|  0.153354 |0|
|0  | 0.21437 | 0.411373 |  0.57505 | 0.692139 | 0.753156 | 0.753156 |  0.69214 | 0.575051 | 0.411374|  0.214371 |0|
|0 | 0.258019 | 0.495136 | 0.692139 |  0.83307  | 0.90651 | 0.906511 | 0.833071  | 0.69214 | 0.495137 |  0.25802|         0|
|0  |0.280765  |0.538785  |0.753156   |0.90651  |0.986425 | 0.986425 | 0.906511 | 0.753157 | 0.538786 | 0.280766    |     0|
|0  |0.280765  |0.538785  |0.753156  |0.906511  |0.986425 | 0.986426 |0.906512  |0.753157  |0.538786 | 0.280766     |    0|
|0   |0.25802  |0.495136  | 0.69214  |0.833071  |0.906511  |0.906512  |0.833072  |0.692141  |0.495137   |0.25802         |0|
|0  |0.214371  |0.411374  |0.575051   |0.69214  |0.753157  |0.753157  |0.692141  |0.575052  |0.411375  |0.214371         |0|
|0  |0.153354  |0.294285  |0.411374  |0.495137  |0.538786  |0.538786  |0.495137  |0.411375  |0.294286  |0.153355         |0|
|0 |0.0799143  |0.153354  |0.214371   |0.25802  |0.280766  |0.280766   |0.25802  |0.214371  |0.153355 |0.0799145         |0|
|0       |  0         |0       |  0    |     0  |       0   |      0   |      0    |     0   |      0    |   0  |0         |0|
