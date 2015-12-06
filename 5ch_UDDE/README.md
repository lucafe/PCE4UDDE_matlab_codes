# Uncertain Delay Differential Equation

In this folder, you can find the codes that are used in the fifth chapter of the Master thesis: 
*“Uncertain Delay Differential Equation”*.
In particular, you can find:

 - **Coefficients\_Legendre**[*script*] evaluates the different inner products used for the Galerkin projection (it must be run once in order to evaluate them).

- **Hayes\_UDDE** [*function*] describes the Hayes Uncertain Delay Differential Equations and uses the coefficients evaluated in **Coefficients\_Legendre**.

-**Indexes** [*function*] evaluates the indexes to construct a bi dimensional polynomial basis

- **MAIN\_FILE** [*script*] Solve the Hayes UDDE defined in the function **Hayes\_UDDE** and require to run at least once the script **Coefficients\_Legendre**.

In these folder, there are also the Matlab codes
-**GaussLegendre_sh**
-**legendre_sh**
Described in [3ch_Approximation](https://github.com/lucafe/PCE4UDDE_matlab_codes/3ch_Approximation).


[Return to the main folder](https://github.com/lucafe/PCE4UDDE_matlab_codes).