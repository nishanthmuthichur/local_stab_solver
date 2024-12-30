import numpy as np

import ls_utils as lsu

PI = np.pi

DUMMY = -10000

class loc_cls():
    
    def __init__(self):
        
        self.nu     = DUMMY
        
        self.r_co   = DUMMY        
        self.U_mean = DUMMY

        self.m_azi  = DUMMY
        self.k_axi  = DUMMY
   
def cons_A_mat(Ublk):
    
    N_rad = len(Ublk.r_co)
    
    m_azi = Ublk.m_azi
    k_axi = Ublk.k_axi
    nu    = Ublk.nu
    
    A_mat =      np.zeros([4 * N_rad, 4 * N_rad]) + \
            1j * np.zeros([4 * N_rad, 4 * N_rad])
            
    Dr_by_r_mat = np.zeros([N_rad, N_rad])
    
    Dr_mat  = lsu.get_D_mat(N_rad)
    Drr_mat = lsu.get_DD_mat(N_rad)
    
    for idx in range(0, N_rad):
        
        Dr_by_r_mat[idx, :] = Dr_mat[idx, :] / Ublk.r_co[idx]
        
#   Row 1    
    
    A_mat[1 : (N_rad - 1), 0 : N_rad] = Dr_mat[1 : (N_rad - 1), :]
    
    for idx in range(1, (N_rad - 1)):
    
        A_mat[idx, idx]             = A_mat[idx, idx] + 1 / Ublk.r_co[idx]
        
        A_mat[idx, N_rad + idx]     = A_mat[idx, N_rad + idx] + \
                                               ((1j * m_azi) / Ublk.r_co[idx])
        
        A_mat[idx, 2 * N_rad + idx] = A_mat[idx, 2 * N_rad + idx] + \
                                                (1j * k_axi)
        
#   Row 2
    
    A_mat[N_rad + 1 : (2 * N_rad - 1), 0 : N_rad] =  \
                                (1j * nu) *     Drr_mat[1 : (N_rad - 1), :] + \
                                (1j * nu) * Dr_by_r_mat[1 : (N_rad - 1), :]  
    
    A_mat[N_rad + 1 : (2 * N_rad - 1), 3 * N_rad : 4 * N_rad] = -1j * Dr_mat[1 : (N_rad - 1), :]
    
    for idx in range(1, (N_rad - 1)):
        
        A_mat[N_rad + idx, idx] = A_mat[N_rad + idx, idx]                      + \
                                  k_axi * Ublk.uz_mean[idx]                    - \
                                  ((1j * nu * (m_azi**2))/(Ublk.r_co[idx]**2)) - \
                                  ((1j * nu) / (Ublk.r_co[idx]**2))            - \
                                  (1j * nu * (k_axi**2))
                    
        A_mat[N_rad + idx, N_rad + idx] = A_mat[N_rad + idx, N_rad + idx] + \
                                          ((2 * m_azi * nu)/(Ublk.r_co[idx]**2))
                                          
#   Row 3

    A_mat[(2 * N_rad + 1) : (3 * N_rad - 1), N_rad : (2 * N_rad)] = \
                                (1j * nu) *     Drr_mat[1 : (N_rad - 1), :] + \
                                (1j * nu) * Dr_by_r_mat[1 : (N_rad - 1), :]       
                            
    for idx in range(1, (N_rad - 1)):
        
        A_mat[2 * N_rad + idx, idx] = -2 * nu * m_azi / (Ublk.r_co[idx]**2)
        
        A_mat[2 * N_rad + idx, N_rad + idx] = \
                                        A_mat[2 * N_rad + idx, N_rad + idx] + \
                                                (k_axi * Ublk.uz_mean[idx]) - \
                                 ((1j * nu * m_azi**2)/(Ublk.r_co[idx]**2)) - \
                                 ((1j * nu)           /(Ublk.r_co[idx]**2)) - \
                                                     (1j * nu * (k_axi**2))
                                  
        A_mat[2 * N_rad + idx, 3 * N_rad + idx] = m_azi / (Ublk.r_co[idx])
        
        
#   Row 4

    A_mat[(3 * N_rad + 1) : (4 * N_rad - 1), (2 * N_rad) : (3 * N_rad)] = \
                                (1j * nu) *     Drr_mat[1 : (N_rad - 1), :] + \
                                (1j * nu) * Dr_by_r_mat[1 : (N_rad - 1), :]       

    for idx in range(1, (N_rad - 1)):
        
        A_mat[(3 * N_rad + idx), idx] = -1j * Ublk.uz_mean_r
        
        A_mat[(3 * N_rad + idx), (2 * N_rad + idx)] = \
                                A_mat[(3 * N_rad + idx), (2 * N_rad + idx)] + \
                                                (k_axi * Ublk.uz_mean[idx]) - \
                                 ((1j * nu * m_azi**2)/(Ublk.r_co[idx]**2)) - \
                                                     (1j * nu * (k_axi**2))
                                  
        A_mat[3 * N_rad + idx, 3 * N_rad + idx] = k_axi
    
    return A_mat
        
        
