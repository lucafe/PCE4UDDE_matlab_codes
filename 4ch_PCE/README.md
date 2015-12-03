# Polynomial Chaos Theory

In this folder, you can find the codes that are used in the fourth chapter of the Master thesis: 
*“Polynomial Chaos Theory”*.
In particular, you can find 4 sub-folders:

- **Examples\_FIELD\_&\_GRIGORIU**  In this folder there are the codes (the new codes are **Field\_Grigoriu\_y1**[*script*] and **Field_Grigoriu_y2**[*script*]) to evaluate the Hermite PCE of the functions:
 Y\_1=\alpha+ exp(\beta(Z)) and 
 Y\_2=abs(Z)
 with Z standard normal distribution. 
The chaos coefficients are evaluated exactly in the article:
R.V. Field and M. Grigoriu, "On the accuracy of the polynomial chaos approximation", *Probabilistic Engineering Mechanics*, 19 (2004) pages 65 - 80. 
The codes are tested to see if the Gaussian Quadrature can fit the exact results.

- **Grids** Construction of the sparse and the tensor grids and corresponding weights. **tensorQUAD** construction of the tensor grid, **sparseQUAD** construction of the sparse grids given the stochastic dimension **D** and order of the construction **l**. The function **sparseQUAD** was inspired by the file [spquad](http://it.mathworks.com/matlabcentral/fileexchange/19063-sparse-grid-quadrature/content/spquad.m)


- **Hermite_PCE** Hermite PCE of standard normal distribution (**Hermite\_PCE\_of\_Normal\_Distribution**[*script*]) and exponential distribution (**Hermite\_PCE\_of\_Exponential\_Distribution**[*script*]). In particular for both distributions we evaluate the Chaos Coefficients, for the Exponential Distribution we evaluate also the pdf estimated by the PCE up to 4th order. 


- **Stochastic_ODE** The main file is **MAIN\_FILE** [*script*] and computes the solution of the stochastic process, solution of the stochastic Maltusian ODE with random parameter normally distributed. It uses the built in function **ode23** with the function **stochastic_ode_PCE**[*function*]

In these folders, there are also the Matlab codes
-**GaussHermite**
-**GaussLegendre_sh**
-**hermite**
Described in [3ch_Approximation](https://github.com/lucafe/PCE4UDDE_matlab_codes/3ch_Approximation).

[Return to the main folder](https://github.com/lucafe/PCE4UDDE_matlab_codes).