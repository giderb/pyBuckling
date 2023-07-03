# pyBuckling
This class will perform buckling calculations based on a generalized Ritz Energy solution for doubly curved anisotropic (orthotropic and isotropic as well) rectangular plates under generalized boundary conditions.
Inputs:
- xls_path : Excel path to the input file.
- a : Panel dimension - a. See input file more details.
- b : Panel dimension - b. See input file more details.
- thickness : Panel total thickness (mm).
- xbc : Boundary condition for edges along x axis. Options: "SS", "SF", "CC", "CS", "CF", "FF".
- ybc : Boundary condition for edges along y axis. Options: "SS", "SF", "CC", "CS", "CF", "FF".
-         First letter represents the first edge condition, whilst the second parallel edge is represented by second letter:
         Letter codes: "S" for simply supported, "C" for clamped, "F" for free.
         For instance, "SF" is one edge simply supported and the other parallel edge is unconstrained (free)
- Rx : Radius of curvature along x axis (mm).
- Ry : Radius of curvature along y axis (mm).
- percent_fix : Percent fixity for specifying fixity ratio between all edges simply supported and all edges clamped boundary conditions
-         Should be between 0 and 100. It specifies the panel buckling boundary condition in between simple and fixed. 
               If set to 0, the percent fixity is not considered and the specified boundary conditions are used to obtain buckling margin of safety. 
               For instance, in order to achieve 50% fixity, the buckling program is run with SS conditions and again with CC boundary conditions and the two critical buckling loads are averaged.
- num_terms : Buckling half wave number and the number of terms for Ritz solution. 12 is generally sufficient to cover all possible combinations.
             Can be lowered down to 6 or even 4, if the number of half waves are known.

This module reads the inputs from an Excel file, which includes "ABD_Matrix" and "Loads" sheets. 6x6 ABD matrix should be entered in the upperleft corner of the "ABD_Matrix" sheet. Several load cases can be inputted in the "Loads" page, with respect to the instructions provided in the "Conventions" sheet. An example input file is provided in the distribution - "Inputs.xlsx".
