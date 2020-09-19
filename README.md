# Ellptic-or-Hyperbolic-Localization-Using-Minimum-Measurement-Solutions
MATLAB Processing Code for "Elliptic or Hyperbolic Localization Using Minimum Measurement Solutions"

If you use any of the following codes in your research, please cite the corresponding paper as a reference in your publication. Thank you!

# Project Abstract

Localization of an object using a number of sensors is often challenged by outlier observations and solution finding. This project aims at deriving a new algebraic positioning solution using a minimum number of measurements, and from which to develop an outlier detector and an object location estimator. Two measurements are sufficient in 2-D and three in 3-D to yield a solution if they are consistent. The derived minimum measurement solution is exact and reduces the computation to the roots of a quadratic equation. The solution derivation leads to simple criteria to ascertain if the line of positions from two measurements intersect. The intersection condition enables us to establish an outlier detector based on graph processing. By partitioning the overdetermined set of measurements first to obtain the individual minimum measurement solutions, we propose a best linear unbiased estimator (BLUE) to form the final location estimate. Analysis supports the proposed estimator in reaching the CRLB accuracy under Gaussian noise. A measurement partitioning scheme is developed to improve performance when the noise level becomes large. We mainly use elliptic time delay measurements for presentation, and the derived results are applicable to the hyperbolic time difference measurements as well. Both the 2-D and 3-D scenarios are considered.

# Code Description

IndvLocSol.m: Individual (MinMsr) Solution

Config_Comb.m: Measurement Combinations for Individual Solutions

SolDetect.m: Individual Solution Ambiguity Elimination

IndvSolSel.m: Individual Solution Selection for BLUE

BLUEest.m: Combining Individual Solutions by BLUE

Example_Fig12Fig13.m: Example Reproduction for Ellpitic Localization

Example_Fig15Fig16.m: Example Reproduction for Hyperbolic Localization

# Reference

Sanaa S. A. Al-Samahi, Yang Zhang, K. C. Ho, "Elliptic and hyperbolic localizations using minimum measurement solutions," Elsevier Signal Process., vol. 167, Feb. 2020.

https://cisp.ece.missouri.edu/index.html
