import numpy as np

def get_cheb_Dmat(x_vec):
    
    N_pts = len(x_vec)
    
    N = N_pts - 1
    
    C_mat = np.zeros([N_pts, N_pts])
    C_mat = np.asmatrix(C_mat)
    
    for i_idx in range(0, N_pts):
      for j_idx in range(0, N_pts):
          
          if (((i_idx > 0) and (i_idx < N)) \
                                        and \
              ((j_idx > 0) and (j_idx < N))):
              
              if (i_idx != j_idx):
                  
                  C_mat[i_idx, j_idx] = ((-1)**(i_idx + j_idx)) / (x_vec[i_idx] - x_vec[j_idx])

              else:
                  
                  C_mat[i_idx, j_idx] = -x_vec[j_idx] / (2 * (1 - (x_vec[j_idx])**2))
                  
          elif ((i_idx == 0) and ((j_idx > 0) and (j_idx < N))):

                  C_mat[i_idx, j_idx] =  2 * (((-1)**j_idx) / (1 - x_vec[j_idx]))
                  
          elif ((i_idx == N) and ((j_idx > 0) and (j_idx < N))):
              
                  C_mat[i_idx, j_idx] = -2 * (((-1)**(N + j_idx)) / (1 + x_vec[j_idx]))
                  
          elif ((j_idx == 0) and ((i_idx > 0) and (i_idx < N))):
              
                  C_mat[i_idx, j_idx] = -(1/2) * (((-1)**i_idx) / (1 - x_vec[i_idx]))
                  
          elif ((j_idx == N) and ((i_idx > 0) and (i_idx < N))):
              
                  C_mat[i_idx, j_idx] =  (1/2) * (((-1)**(N + i_idx)) / (1 + x_vec[i_idx])) 
             
          elif ((i_idx == 0) and (j_idx == 0)):
              
                  C_mat[i_idx, j_idx] = (2 * (N**2) + 1) / 6
                  
          elif ((i_idx == N) and (j_idx == 0)):
              
                  C_mat[i_idx, j_idx] = -(1/2) * (-1)**N
                  
          elif ((i_idx == 0) and (j_idx == N)):
              
                  C_mat[i_idx, j_idx] = (1/2) * (-1)**N
                  
          elif ((i_idx == N) and (j_idx == N)):
              
                  C_mat[i_idx, j_idx] = -(2 * (N**2) + 1) / 6
                
             
          else: 
              
                  print('Error! All possible if--else--if conditions are exhausted.')
                  
    return C_mat
             
                
             
                
             
                
             
                
             
                
             
                
             
                
             
                
             
                
             
                
             