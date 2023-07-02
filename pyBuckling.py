import numpy as np
import pandas as pd
from scipy.linalg import eigvals
from scipy.integrate import trapezoid
import os
   
class ritz_buckling:
    """
    ritz_buckling class:
        This class will perform buckling calculations based on a generalized Ritz Energy solution.
        Inputs:
            xls_path: Excel path to the input file.
            a: Panel dimension - a. See input file more details.
            b: Panel dimension - b. See input file more details.
            thickness: Panel total thickness (mm).
            xbc: Boundary condition for edges along x axis. Options: "SS", "SF", "CC", "CS", "CF", "FF".
            ybc: Boundary condition for edges along y axis. Options: "SS", "SF", "CC", "CS", "CF", "FF".
                First letter represents the first edge condition, 
                whilst the second parallel edge is represented by second letter:
                    Letter codes: "S" for simply supported, "C" for clamped, "F" for free.
                    For instance, "SF" is one edge simply supported and the other parallel edge is unconstrained (free)
            Rx: Radius of curvature along x axis (mm).
            Ry: Radius of curvature along y axis (mm).
            percent_fix: Percent fixity for specifying fixity ratio between all edges simply supported and all edges clamped boundary conditions
                * Should be between 0 and 100. It specifies the panel buckling boundary condition in between simple and fixed. 
                 If set to 0, the percent fixity is not considered and the specified boundary conditions are used to obtain buckling margin of safety. 
                 For instance, in order to achieve 50% fixity, the buckling program is run with SS conditions and again with CC boundary conditions and the two critical buckling loads are averaged.
            num_terms: Buckling half wave number and the number of terms for Ritz solution. 12 is generally sufficient to cover all possible combinations.
                Can be lowered down to 6 or even 4, if the number of half waves are known.
    """
    def __init__(self, xls_path, a, b, thickness, xbc, ybc, Rx, Ry, percent_fix, num_terms) :
      self.xls_path = xls_path
      self.a = a
      self.b = b
      self.thickness = thickness
      self.xbc = xbc
      self.ybc = ybc
      self.Rx = Rx
      self.Ry = Ry
      self.percent_fix = percent_fix
      m_num_terms = num_terms
      n_num_terms = num_terms
      xbc_orig = xbc
      ybc_orig = ybc
      
      
      # Read ABD matrix from the Excel input file
      ABD = pd.read_excel(xls_path, sheet_name='ABD_Matrix', header=None)
      A11 = ABD.iloc[0,0]
      A12 = ABD.iloc[0,1]
      A16 = ABD.iloc[0,2]
      A22 = ABD.iloc[1,1]
      A26 = ABD.iloc[1,2]
      A66 = ABD.iloc[2,2]
      
      D11 = ABD.iloc[3,3]
      D12 = ABD.iloc[3,4]
      D16 = ABD.iloc[3,5]
      D22 = ABD.iloc[4,4]
      D26 = ABD.iloc[4,5]
      D66 = ABD.iloc[5,5]
      
      # Read loads from the Excel input file
      Loads = pd.read_excel(xls_path, sheet_name='Loads')
      
      # Convert stiffness parameters into equivalent parameters to have integrals ranging from 0 to 1    
      alpha,beta,gamma,delta,vD,A11_p,A12_p,A22_p,A16_p,A26_p,A66_p,Zx,Zy = self.equiv_param(A11,A12,A22,A16,A26,A66,D11,D12,D22,D16,D26,D66,thickness,a,b,Rx,Ry)
      
      # Begin calculations
      U = np.zeros((3 * m_num_terms * n_num_terms, 3 * m_num_terms * n_num_terms))
      
      # If percent fixity is not specified, do the following
      if percent_fix == 0:
      
          # Calculate the strain energy of the laminate
          for ma in range(0, m_num_terms):
              for na in range(0, n_num_terms):
                  for mb in range(0, m_num_terms):
                      for nb in range(0, n_num_terms):
                          i = ma + na * m_num_terms
                          j = mb + nb * m_num_terms
                          if i>=j:
                              func = self.u(ma+1, na+1, mb+1, nb+1, xbc, ybc)
                              # U terms
                              U[i, j] = 12 * (A11_p / alpha**2 * func[6] + A66_p * func[7] + A16_p / alpha * (func[8] + func[9]))
                              U[i, j + m_num_terms * n_num_terms] = 12 * (A12_p * func[9] + A16_p / alpha * func[6] + alpha * A26_p * func[7] + A66_p * func[8])
                              U[i, j + 2 * m_num_terms * n_num_terms] = 12 * ((A11_p / alpha**2 * Zx + A12_p * Zy) * func[13] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[14])
                              # V terms
                              U[i + m_num_terms * n_num_terms, j] = 12 * (A12_p * func[8] + A16_p / alpha * func[6] + alpha * A26_p * func[7] + A66_p * func[9])
                              U[i + m_num_terms * n_num_terms, j + m_num_terms * n_num_terms] = 12 * (alpha**2 * A22_p * func[7] + A66_p * func[6] + alpha * A26_p * (func[8] + func[9]))
                              U[i + m_num_terms * n_num_terms, j + 2 * m_num_terms * n_num_terms] = 12 * ((A12_p * Zx + alpha**2 * A22_p * Zy) * func[14] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[13])
                              # W terms
                              U[i + 2 * m_num_terms * n_num_terms, j] = 12 * ((A11_p / alpha**2 * Zx + A12_p * Zy) * func[11] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[12])
                              U[i + 2 * m_num_terms * n_num_terms, j + m_num_terms * n_num_terms] = 12 * ((A12_p * Zx + alpha**2 * A22_p * Zy) * func[12] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[11])
                              U[i + 2 * m_num_terms * n_num_terms, j + 2 * m_num_terms * n_num_terms] = (12 * (A11_p / alpha**2 * Zx**2 + 2 * A12_p * Zx * Zy + alpha**2 * A22_p * Zy**2)) * func[10] \
                                          + 1 / alpha**2 * func[1] + alpha**2 * func[2] + vD * (func[3] + func[4]) + 2 * (beta - vD) * func[5] \
                                          + 2 * gamma / alpha * (func[15] + func[16]) + 2 * alpha * delta * (func[17] + func[18])
          
          U = self.symmetrize(U)
      
      # If percent fixity is specified, do the following
      else:
          # Do the calculations for the simply supported BC first
          xbc, ybc = "SS", "SS"
          # Calculate the strain energy of the laminate
          for ma in range(0, m_num_terms):
              for na in range(0, n_num_terms):
                  for mb in range(0, m_num_terms):
                      for nb in range(0, n_num_terms):
                          i = ma + na * m_num_terms
                          j = mb + nb * m_num_terms
                          if i>=j:
                              func = self.u(ma+1, na+1, mb+1, nb+1, xbc, ybc)
                              # U terms
                              U[i, j] = 12 * (A11_p / alpha**2 * func[6] + A66_p * func[7] + A16_p / alpha * (func[8] + func[9]))
                              U[i, j + m_num_terms * n_num_terms] = 12 * (A12_p * func[9] + A16_p / alpha * func[6] + alpha * A26_p * func[7] + A66_p * func[8])
                              U[i, j + 2 * m_num_terms * n_num_terms] = 12 * ((A11_p / alpha**2 * Zx + A12_p * Zy) * func[13] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[14])
                              # V terms
                              U[i + m_num_terms * n_num_terms, j] = 12 * (A12_p * func[8] + A16_p / alpha * func[6] + alpha * A26_p * func[7] + A66_p * func[9])
                              U[i + m_num_terms * n_num_terms, j + m_num_terms * n_num_terms] = 12 * (alpha**2 * A22_p * func[7] + A66_p * func[6] + alpha * A26_p * (func[8] + func[9]))
                              U[i + m_num_terms * n_num_terms, j + 2 * m_num_terms * n_num_terms] = 12 * ((A12_p * Zx + alpha**2 * A22_p * Zy) * func[14] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[13])
                              # W terms
                              U[i + 2 * m_num_terms * n_num_terms, j] = 12 * ((A11_p / alpha**2 * Zx + A12_p * Zy) * func[11] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[12])
                              U[i + 2 * m_num_terms * n_num_terms, j + m_num_terms * n_num_terms] = 12 * ((A12_p * Zx + alpha**2 * A22_p * Zy) * func[12] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[11])
                              U[i + 2 * m_num_terms * n_num_terms, j + 2 * m_num_terms * n_num_terms] = (12 * (A11_p / alpha**2 * Zx**2 + 2 * A12_p * Zx * Zy + alpha**2 * A22_p * Zy**2)) * func[10] \
                                          + 1 / alpha**2 * func[1] + alpha**2 * func[2] + vD * (func[3] + func[4]) + 2 * (beta - vD) * func[5] \
                                          + 2 * gamma / alpha * (func[15] + func[16]) + 2 * alpha * delta * (func[17] + func[18])
          
          U1 = self.symmetrize(U)
          # Do the calculations for clamped BC then
          U = np.zeros((3 * m_num_terms * n_num_terms, 3 * m_num_terms * n_num_terms))
          xbc, ybc = "CC", "CC"
          # Calculate the strain energy of the laminate
          for ma in range(0, m_num_terms):
              for na in range(0, n_num_terms):
                  for mb in range(0, m_num_terms):
                      for nb in range(0, n_num_terms):
                          i = ma + na * m_num_terms
                          j = mb + nb * m_num_terms
                          if i>=j:
                              func = self.u(ma+1, na+1, mb+1, nb+1, xbc, ybc)
                              # U terms
                              U[i, j] = 12 * (A11_p / alpha**2 * func[6] + A66_p * func[7] + A16_p / alpha * (func[8] + func[9]))
                              U[i, j + m_num_terms * n_num_terms] = 12 * (A12_p * func[9] + A16_p / alpha * func[6] + alpha * A26_p * func[7] + A66_p * func[8])
                              U[i, j + 2 * m_num_terms * n_num_terms] = 12 * ((A11_p / alpha**2 * Zx + A12_p * Zy) * func[13] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[14])
                              # V terms
                              U[i + m_num_terms * n_num_terms, j] = 12 * (A12_p * func[8] + A16_p / alpha * func[6] + alpha * A26_p * func[7] + A66_p * func[9])
                              U[i + m_num_terms * n_num_terms, j + m_num_terms * n_num_terms] = 12 * (alpha**2 * A22_p * func[7] + A66_p * func[6] + alpha * A26_p * (func[8] + func[9]))
                              U[i + m_num_terms * n_num_terms, j + 2 * m_num_terms * n_num_terms] = 12 * ((A12_p * Zx + alpha**2 * A22_p * Zy) * func[14] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[13])
                              # W terms
                              U[i + 2 * m_num_terms * n_num_terms, j] = 12 * ((A11_p / alpha**2 * Zx + A12_p * Zy) * func[11] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[12])
                              U[i + 2 * m_num_terms * n_num_terms, j + m_num_terms * n_num_terms] = 12 * ((A12_p * Zx + alpha**2 * A22_p * Zy) * func[12] + (A16_p / alpha * Zx + alpha * A26_p * Zy) * func[11])
                              U[i + 2 * m_num_terms * n_num_terms, j + 2 * m_num_terms * n_num_terms] = (12 * (A11_p / alpha**2 * Zx**2 + 2 * A12_p * Zx * Zy + alpha**2 * A22_p * Zy**2)) * func[10] \
                                          + 1 / alpha**2 * func[1] + alpha**2 * func[2] + vD * (func[3] + func[4]) + 2 * (beta - vD) * func[5] \
                                          + 2 * gamma / alpha * (func[15] + func[16]) + 2 * alpha * delta * (func[17] + func[18])
          
          U2 = self.symmetrize(U)
      
      # Set the headers for the output CSV file
      headers =  ["Load Case Name", "a (mm)", "b (mm)", "Thickness (mm)", "X Boundary Condition", "Y Boundary Condition", "Percent Fixity", "Radius of Curvature X (mm)",
                      "Radius of Curvature Y (mm)", "RF", "NxU Applied (N/mm)", "NxL Applied (N/mm)", "NyU Applied (N/mm)", "NyL Applied (N/mm)", "Nxy Applied (N/mm)", 
                      "NxU Critical (N/mm)", "NxL Critical (N/mm)", "NyU Critical (N/mm)", "NyL Critical (N/mm)", "Nxy Critical (N/mm)",
                      "A11", "A12", "A16", "A22", "A26", "A66", 
                      "D11", "D12", "D16", "D22", "D26", "D66" ]
      N_array = np.empty((0, len(headers)), float)
      N_array_DF = pd.DataFrame(N_array)
      N_array_DF.columns = headers
      
      # Set the output location for the output CSV file
      csv_filename = os.path.join(os.path.dirname(xls_path),''.join(['Ritz_Buckling.csv']))
      
      if os.path.exists(csv_filename):
              os.remove(csv_filename)
              
      N_array_DF.to_csv(csv_filename, mode='a', header=True, index=False)
              
      N_array=[]
      
      # Calculate the potential energy of the applied loads
      for index, row in Loads.iterrows():
          W = np.zeros((3 * m_num_terms * n_num_terms, 3 * m_num_terms * n_num_terms))
          LC_Name = row["Load Case Name"]
          
          # Multiply all loads by -1, as the Ritz solution requires compressive loads to be positive
          NxU = row["NxU"] * -1
          NxL = row["NxL"] * -1
          
          if NxU < NxL or NxL == 0:
              NxU = row["NxL"] * -1
              NxL = row["NxU"] * -1
              
          NyU = row["NyU"] * -1
          NyL = row["NyL"] * -1
          
          if NyU < NxL or NyL == 0:
              NyU = row["NyL"] * -1
              NyL = row["NyU"] * -1    
          
          if NxL !=0:
              sigma_x = (NxL - NxU) / NxL
          else:
              sigma_x = 0
          
          if NyL !=0:
              sigma_y = (NyL - NyU) / NyL
          else:
              sigma_y = 0
          
          Nxy = row["Nxy"] * -1
          
          # If percent fixity is not specified, do the following
          if percent_fix == 0: 
              for ma in range(0, m_num_terms):
                  for na in range(0, n_num_terms):
                      for mb in range(0, m_num_terms):
                          for nb in range(0, n_num_terms):
                              i = ma + na * m_num_terms
                              j = mb + nb * m_num_terms
                              if i>=j:
                                  func = self.w(ma+1, na+1, mb+1, nb+1, xbc, ybc, sigma_x, sigma_y)
                                  W[i + 2 * m_num_terms * n_num_terms, j + 2 * m_num_terms * n_num_terms] = np.pi**2 * (NxL * func[1] +  alpha**2 * NyL * func[2] - alpha * Nxy * (func[3] + func[4]))
              W = self.symmetrize(W)
              
              # Calculate eigenvalues for the current load case and find the minimum positive eigenvalue
              eigenvalues = eigvals(U, W)
              eig_real=eigenvalues.real
              eig_real=min(eig_real[eig_real>0])
          # If percent fixity is specified, do the following
          else:
              # Do the calculations for the simply supported BC first
              xbc, ybc = "SS", "SS"
              for ma in range(0, m_num_terms):
                  for na in range(0, n_num_terms):
                      for mb in range(0, m_num_terms):
                          for nb in range(0, n_num_terms):
                              i = ma + na * m_num_terms
                              j = mb + nb * m_num_terms
                              if i>=j:
                                  func = self.w(ma+1, na+1, mb+1, nb+1, xbc, ybc, sigma_x, sigma_y)
                                  W[i + 2 * m_num_terms * n_num_terms, j + 2 * m_num_terms * n_num_terms] = np.pi**2 * (NxL * func[1] +  alpha**2 * NyL * func[2] - alpha * Nxy * (func[3] + func[4]))
              W1 = self.symmetrize(W)
              # Do the calculations for clamped BC then
              xbc, ybc = "CC", "CC"
              W = np.zeros((3 * m_num_terms * n_num_terms, 3 * m_num_terms * n_num_terms))
              for ma in range(0, m_num_terms):
                  for na in range(0, n_num_terms):
                      for mb in range(0, m_num_terms):
                          for nb in range(0, n_num_terms):
                              i = ma + na * m_num_terms
                              j = mb + nb * m_num_terms
                              if i>=j:
                                  func = self.w(ma+1, na+1, mb+1, nb+1, xbc, ybc, sigma_x, sigma_y)
                                  W[i + 2 * m_num_terms * n_num_terms, j + 2 * m_num_terms * n_num_terms] = np.pi**2 * (NxL * func[1] +  alpha**2 * NyL * func[2] - alpha * Nxy * (func[3] + func[4]))
              W2 = self.symmetrize(W)
              # Calculate eigenvalues for the current load case and find the minimum positive eigenvalue for simply supported boundary condition
              eigenvalues1 = eigvals(U1, W1)
              eig_real1=eigenvalues1.real
              eig_real1=min(eig_real1[eig_real1>0])
              # Calculate eigenvalues for the current load case and find the minimum positive eigenvalue for clamped boundary condition
              eigenvalues2 = eigvals(U2, W2)
              eig_real2=eigenvalues2.real
              eig_real2=min(eig_real2[eig_real2>0])
              eig_real = (eig_real1 * (100 - percent_fix) + eig_real2 * percent_fix) / 100
          
          # Find critical buckling load in each direction using the minimum eigenvalue
          Nx_cr = eig_real * NxL * (np.pi / b)**2 * np.sqrt(D11 * D22)
          Ny_cr = eig_real * NyL * (np.pi / b)**2 * D22
          Nxy_cr = eig_real * Nxy * (np.pi / b)**2 * (D11 * D22**3)**0.25
          
          # Calculate minimum reserve factors
          RF_x = Nx_cr / NxL if NxL != 0 else 1e20
          RF_y = Ny_cr / NyL  if NyL != 0 else 1e20
          RF_xy = Nxy_cr / Nxy  if Nxy != 0 else 1e20
          RFmin = min(RF_x, RF_y, RF_xy)
          
          # Revert back to original sign of the loading
          NxcrL = RFmin * -NxL
          NxcrU = RFmin * -NxU
          NycrL = RFmin * -NyL
          NycrU = RFmin * -NyU
          Nxycr = RFmin * -Nxy  
          
          # Append output to the array 
          N_array.append([LC_Name, a, b, thickness, xbc_orig, ybc_orig, percent_fix, Rx, Ry, RFmin, row["NxU"], row["NxL"], row["NyU"], row["NyL"], row["Nxy"], NxcrU, NxcrL, NycrU, NycrL, Nxycr,
                        A11, A12, A16, A22, A26, A66, D11, D12, D16, D22, D26, D66])         
      
      # Output output array to CSV file
      N_array_DF = pd.DataFrame(N_array)
      N_array_DF.to_csv(csv_filename, mode='a', header=False, index=False)
      
    def symmetrize(self, a):
        return a + a.T - np.diag(a.diagonal())
    
    def equiv_param(self,A11,A12,A22,A16,A26,A66,D11,D12,D22,D16,D26,D66,thickness,a,b,Rx,Ry):
        alpha = a / b * (D22 / D11)**0.25
        beta = (D12 + 2 * D66) / (D11 * D22)**0.5
        gamma = D16 / (D11**3 * D22)**0.25
        delta = D26 / (D11 * D22**3)**0.25
        vD = D12 / (D11 * D22)**0.5
        A11_p = A11 / D11 * thickness**2 / 12
        A12_p = A12 / (D11 * D22)**0.5 * thickness**2 / 12
        A22_p = A22 / D22 * thickness**2 / 12
        A16_p = A16 / (D11**3 * D22)**0.25 * thickness**2 / 12
        A26_p = A26 / (D11 * D22**3)**0.25 * thickness**2 / 12
        A66_p = A66 / (D11 * D22)**0.5 * thickness**2 / 12
        if Rx != 0:
            Zx = a**2 / (thickness * Rx) 
        else :
            Zx = 0
        
        if Ry != 0:
            Zy = b**2 / (thickness * Ry)
        else:
            Zy = 0
            
        return alpha,beta,gamma,delta,vD,A11_p,A12_p,A22_p,A16_p,A26_p,A66_p,Zx,Zy
    
    def bc(self, m, n, xbc, ybc):
        
        # X boundary condition
        if xbc == "SS":
            # simply supported 
            B1 = (m) * np.pi
            C1, B2, C2 = 0, 0, 0
        elif xbc == "CC":
            # clamped - clamped
            B1 = (m + 1) * np.pi
            C1 = -np.pi / 2
            B2 = (m - 1) * np.pi
            C2 = np.pi / 2
        elif xbc == "CS":
            # clamped - simply supported
            B1 = (m + 1 / 2) * np.pi
            C1 = -np.pi / 2
            B2 = (m - 1 / 2) * np.pi
            C2 = np.pi / 2
        elif xbc == "SF":
            # simply supported - free
            B1 = (m / 2 - 1 / 4) * np.pi
            C1, B2, C2 = 0, 0, 0
        elif xbc == "CF":
            # clamped - free
            B1 = (m / 2 + 1 / 4) * np.pi
            C1 = -np.pi / 2
            B2 = (m / 2 - 3 / 4) * np.pi
            C2 = np.pi / 2
        elif xbc == "FF":
            # free - free
            if (m % 2) == 0:
                B1 = ((m - 1) / 2) * np.pi
                C1 = -((m + 1) / 4) * np.pi
                B2 = 0
                C2 = 0
            else:
                B1 = (m / 2) * np.pi
                C1 = -(m / 4) * np.pi
                B2 = 0
                C2 = 0
            
        # Y boundary condition    
        if ybc == "SS":
            # simply supported 
            D1 = (n) * np.pi
            E1, D2, E2 = 0, 0, 0
        elif ybc == "CC":
            # clamped - clamped
            D1 = (n + 1) * np.pi
            E1 = -np.pi / 2
            D2 = (n - 1) * np.pi
            E2 = np.pi / 2
        elif ybc == "CS":
            # clamped - simply supported
            D1 = (n + 1 / 2) * np.pi
            E1 = -np.pi / 2
            D2 = (n - 1 / 2) * np.pi
            E2 = np.pi / 2
        elif ybc == "SF":
            # simply supported - free
            D1 = (n / 2 - 1 / 4) * np.pi
            E1, D2, E2 = 0, 0, 0
        elif ybc == "CF":
            # clamped - free
            D1 = (n / 2 + 1 / 4) * np.pi
            E1 = -np.pi / 2
            D2 = (n / 2 - 3 / 4) * np.pi
            E2 = np.pi / 2
        elif ybc == "FF":
            # free - free
            if (n % 2) == 0:
                D1 = ((n - 1) / 2) * np.pi
                E1 = -((n + 1) / 4) * np.pi
                D2 = 0
                E2 = 0
            else:
                D1 = (n / 2) * np.pi
                E1 = -(n / 4) * np.pi
                D2 = 0
                E2 = 0
            
        return B1, C1, B2, C2, D1, E1, D2, E2
    
    def u(self, ma, na, mb, nb, xbc, ybc):
        func = {}
        
        B1a, C1a, B2a, C2a, D1a, E1a, D2a, E2a = self.bc(ma, na, xbc, ybc)
        B1b, C1b, B2b, C2b, D1b, E1b, D2b, E2b = self.bc(mb, nb, xbc, ybc)
            
        tt = np.linspace(0,1,50)
        
    #   SS d2_fwi/dx^2 * d2_fwj/dx^2 
        intx = trapezoid( (-B1a**2 * np.sin(B1a*tt+C1a) - B2a**2 * np.sin(B2a*tt+C2a)) * (-B1b**2 * np.sin(B1b*tt+C1b) - B2b**2 * np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (np.sin(D1a*tt+E1a) + np.sin(D2a*tt+E2a)) * (np.sin(D1b*tt+E1b) + np.sin(D2b*tt+E2b)) , tt)    
        func[1] = intx*inty
        
    #   SS d2_fwi/dy^2 * d2_fwj/dy^2 
        intx = trapezoid( (np.sin(B1a*tt+C1a) + np.sin(B2a*tt+C2a)) * (np.sin(B1b*tt+C1b) + np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (-D1a**2 * np.sin(D1a*tt+E1a) - D2a**2 * np.sin(D2a*tt+E2a)) * (-D1b**2 * np.sin(D1b*tt+E1b) - D2b**2 * np.sin(D2b*tt+E2b)) , tt)    
        func[2] = intx*inty
        
    #   SS d2_fwi/dy^2 * d2_fwj/dx^2    
        intx = trapezoid( (np.sin(B1a*tt+C1a) + np.sin(B2a*tt+C2a)) * (-B1b**2 * np.sin(B1b*tt+C1b) - B2b**2 * np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (-D1a**2 * np.sin(D1a*tt+E1a) - D2a**2 * np.sin(D2a*tt+E2a)) * (np.sin(D1b*tt+E1b) + np.sin(D2b*tt+E2b)) , tt)    
        func[3] = intx*inty
        
    #   SS d2_fwi/dx^2 * d2_fwj/dy^2         
        intx = trapezoid( (-B1a**2 * np.sin(B1a*tt+C1a) - B2a**2 * np.sin(B2a*tt+C2a)) * (np.sin(B1b*tt+C1b) + np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (np.sin(D1a*tt+E1a) + np.sin(D2a*tt+E2a)) * (-D1b**2 * np.sin(D1b*tt+E1b) - D2b**2 * np.sin(D2b*tt+E2b)) , tt)    
        func[4] = intx*inty
        
    #   SS d2_fwi/dx.dy * d2_fwj/dx.dy     
        intx = trapezoid( (B1a * np.cos(B1a*tt+C1a) + B2a * np.cos(B2a*tt+C2a)) * (B1b * np.cos(B1b*tt+C1b) + B2b * np.cos(B2b*tt+C2b))  , tt) 
        inty = trapezoid( (D1a * np.cos(D1a*tt+E1a) + D2a * np.cos(D2a*tt+E2a))  *  (D1b * np.cos(D1b*tt+E1b) + D2b * np.cos(D2b*tt+E2b)) , tt)    
        func[5] = intx*inty
        
    #   SS d_fwi/dx  * d_fwj/dx
        intx = trapezoid( (B1a * np.cos(B1a*tt+C1a) + B2a * np.cos(B2a*tt+C2a)) *(B1b * np.cos(B1b*tt+C1b) + B2b * np.cos(B2b*tt+C2b))  , tt) 
        inty = trapezoid( (np.sin(D1a*tt+E1a) + np.sin(D2a*tt+E2a)) * (np.sin(D1b*tt+E1b) + np.sin(D2b*tt+E2b)) , tt)    
        func[6] = intx*inty
    
    #   SS d_fwi/dy * d_fwj/dy
        intx = trapezoid( (np.sin(B1a*tt+C1a) + np.sin(B2a*tt+C2a)) * (np.sin(B1b*tt+C1b) + np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (D1a * np.cos(D1a*tt+E1a) + D2a * np.cos(D2a*tt+E2a)) * (D1b * np.cos(D1b*tt+E1b) + D2b * np.cos(D2b*tt+E2b)) , tt)    
        func[7] = intx*inty
    
    #   SS d_fwi/dx * d_fwj/dy
        intx = trapezoid( (B1a * np.cos(B1a*tt+C1a) + B2a * np.cos(B2a*tt+C2a)) * (np.sin(B1b*tt+C1b) + np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (np.sin(D1a*tt+E1a) + np.sin(D2a*tt+E2a)) * (D1b * np.cos(D1b*tt+E1b) + D2b * np.cos(D2b*tt+E2b)) , tt)    
        func[8] = intx*inty
    
    #   SS d_fwi/dy * d_fwj/dx   
        intx = trapezoid( (np.sin(B1a*tt+C1a) + np.sin(B2a*tt+C2a)) * (B1b * np.cos(B1b*tt+C1b) + B2b * np.cos(B2b*tt+C2b))  , tt) 
        inty = trapezoid( (D1a * np.cos(D1a*tt+E1a) + D2a * np.cos(D2a*tt+E2a)) * (np.sin(D1b*tt+E1b) + np.sin(D2b*tt+E2b)) , tt)    
        func[9] = intx*inty
    
    #   SS fwi * fwj   
        intx = trapezoid( (np.sin(B1a*tt+C1a) + np.sin(B2a*tt+C2a)) * (np.sin(B1b*tt+C1b) + np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (np.sin(D1a*tt+E1a) + np.sin(D2a*tt+E2a)) * (np.sin(D1b*tt+E1b) + np.sin(D2b*tt+E2b)) , tt)    
        func[10] = intx*inty
        
    #   SS d_fui/dx * fwj, d_fvi/dx * fwj  
        intx = trapezoid( (B1a * np.cos(B1a*tt+C1a) + B2a * np.cos(B2a*tt+C2a)) * (np.sin(B1b*tt+C1b) + np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (np.sin(D1a*tt+E1a) + np.sin(D2a*tt+E2a)) * (np.sin(D1b*tt+E1b) + np.sin(D2b*tt+E2b)) , tt)    
        func[11] = intx*inty
    
    #   SS d_fui/dy * fwj, d_fvi/dy * fwj
        intx = trapezoid( (np.sin(B1a*tt+C1a) + np.sin(B2a*tt+C2a)) * (np.sin(B1b*tt+C1b) + np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (D1a * np.cos(D1a*tt+E1a) + D2a * np.cos(D2a*tt+E2a)) * (np.sin(D1b*tt+E1b) + np.sin(D2b*tt+E2b)) , tt)    
        func[12] = intx*inty
    
    #   SS fwi * d_fuj/dx, fwi * d_fvj/dx
        intx = trapezoid( (np.sin(B1a*tt+C1a) + np.sin(B2a*tt+C2a)) * (B1b * np.cos(B1b*tt+C1b) + B2b * np.cos(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (np.sin(D1a*tt+E1a) + np.sin(D2a*tt+E2a)) * (np.sin(D1b*tt+E1b) + np.sin(D2b*tt+E2b)) , tt)    
        func[13] = intx*inty
    
    #   SS fwi * d_fuj/dy, fwi * d_fvj/dy
        intx = trapezoid( (np.sin(B1a*tt+C1a) + np.sin(B2a*tt+C2a)) * (np.sin(B1b*tt+C1b) + np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (np.sin(D1a*tt+E1a) + np.sin(D2a*tt+E2a)) * (D1b * np.cos(D1b*tt+E1b) + D2b * np.cos(D2b*tt+E2b)) , tt)    
        func[14] = intx*inty
        
    #   SS d2_fwi/dx.dy * d2_fwj/dx^2 
        intx = trapezoid( (B1a * np.cos(B1a*tt+C1a) + B2a * np.cos(B2a*tt+C2a)) * (-B1b**2 * np.sin(B1b*tt+C1b) - B2b**2 * np.sin(B2b*tt+C2b))  , tt) 
        inty = trapezoid( (D1a * np.cos(D1a*tt+E1a) + D2a * np.cos(D2a*tt+E2a))  *  (np.sin(D1b*tt+E1b) + np.sin(D2b*tt+E2b)) , tt)   
        func[15] = intx*inty
    
    #   SS d2_fwi/dx^2 * d2_fwj/dx.dy 
        intx = trapezoid( (-B1a**2 * np.sin(B1a*tt+C1a) - B2a**2 * np.sin(B2a*tt+C2a)) * (B1b * np.cos(B1b*tt+C1b) + B2b * np.cos(B2b*tt+C2b)) , tt)
        inty = trapezoid( (np.sin(D1a*tt+E1a) + np.sin(D2a*tt+E2a)) * (D1b * np.cos(D1b*tt+E1b) + D2b * np.cos(D2b*tt+E2b)) , tt)    
        func[16] = intx*inty
        
    #   SS d2_fwi/dx.dy * d2_fwj/dy^2  
        intx = trapezoid( (B1a * np.cos(B1a*tt+C1a) + B2a * np.cos(B2a*tt+C2a)) * (np.sin(B1b*tt+C1b) + np.sin(B2b*tt+C2b))  , tt) 
        inty = trapezoid( (D1a * np.cos(D1a*tt+E1a) + D2a * np.cos(D2a*tt+E2a))  *  (-D1b**2 * np.sin(D1b*tt+E1b) - D2b**2 * np.sin(D2b*tt+E2b)) , tt)   
        func[17] = intx*inty
    
    #   SS d2_fwi/dy^2 * d2_fwj/dx.dy  
        intx = trapezoid( (np.sin(B1a*tt+C1a) + np.sin(B2a*tt+C2a)) * (B1b * np.cos(B1b*tt+C1b) + B2b * np.cos(B2b*tt+C2b)) , tt)
        inty = trapezoid( (-D1a**2 * np.sin(D1a*tt+E1a) - D2a**2 * np.sin(D2a*tt+E2a)) * (D1b * np.cos(D1b*tt+E1b) + D2b * np.cos(D2b*tt+E2b)) , tt)    
        func[18] = intx*inty
        
        return func
    
    
    def w(self, ma, na, mb, nb, xbc, ybc, sigma_x, sigma_y):
        func = {}
        
        B1a, C1a, B2a, C2a, D1a, E1a, D2a, E2a = self.bc(ma, na, xbc, ybc)
        B1b, C1b, B2b, C2b, D1b, E1b, D2b, E2b = self.bc(mb, nb, xbc, ybc)
            
        tt = np.linspace(0,1,50)
            
    #   SS (1 - sigma_x * y) * d_fwi/dx  * d_fwj/dx
        intx = trapezoid( (B1a * np.cos(B1a*tt+C1a) + B2a * np.cos(B2a*tt+C2a)) *(B1b * np.cos(B1b*tt+C1b) + B2b * np.cos(B2b*tt+C2b))  , tt) 
        inty = trapezoid( (1 - sigma_x * tt) * (np.sin(D1a*tt+E1a) + np.sin(D2a*tt+E2a)) * (np.sin(D1b*tt+E1b) + np.sin(D2b*tt+E2b)) , tt)    
        func[1] = intx*inty
    
    #   SS (1 - sigma_y * x) * d_fwi/dy * d_fwj/dy
        intx = trapezoid( (1 - sigma_y * tt) * (np.sin(B1a*tt+C1a) + np.sin(B2a*tt+C2a)) * (np.sin(B1b*tt+C1b) + np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (D1a * np.cos(D1a*tt+E1a) + D2a * np.cos(D2a*tt+E2a)) * (D1b * np.cos(D1b*tt+E1b) + D2b * np.cos(D2b*tt+E2b)) , tt)    
        func[2] = intx*inty
    
    #   SS d_fwi/dx * d_fwj/dy
        intx = trapezoid( (B1a * np.cos(B1a*tt+C1a) + B2a * np.cos(B2a*tt+C2a)) * (np.sin(B1b*tt+C1b) + np.sin(B2b*tt+C2b)) , tt) 
        inty = trapezoid( (np.sin(D1a*tt+E1a) + np.sin(D2a*tt+E2a)) * (D1b * np.cos(D1b*tt+E1b) + D2b * np.cos(D2b*tt+E2b)) , tt)    
        func[3] = intx*inty
    
    #   SS d_fwi/dy * d_fwj/dx   
        intx = trapezoid( (np.sin(B1a*tt+C1a) + np.sin(B2a*tt+C2a)) * (B1b * np.cos(B1b*tt+C1b) + B2b * np.cos(B2b*tt+C2b))  , tt) 
        inty = trapezoid( (D1a * np.cos(D1a*tt+E1a) + D2a * np.cos(D2a*tt+E2a)) * (np.sin(D1b*tt+E1b) + np.sin(D2b*tt+E2b)) , tt)    
        func[4] = intx*inty
        
        return func