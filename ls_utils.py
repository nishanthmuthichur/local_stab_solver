import numpy as np
import scipy as scp

import ls_utils as lsu
import num_sch as ns

PI = np.pi

DUMMY = -10000

class loc_cls():
    
    def __init__(self):
        
        self.nu        = DUMMY

        self.m_azi  = DUMMY
        self.k_axi  = DUMMY
        
        self.r_co      = DUMMY        
        self.uz_mean   = DUMMY
        self.uz_mean_r = DUMMY

        self.Dr_mat      = DUMMY
        self.Drr_mat     = DUMMY        
        self.Dr_by_r_mat = DUMMY

        self.A_mat = DUMMY
        self.B_mat = DUMMY
        
        self.EW = DUMMY
        self.EV = DUMMY



def cons_deriv_mat(Ublk):

    N_rad = len(Ublk.r_co)

    r_co = Ublk.r_co

    Dr_by_r_mat = np.zeros([N_rad, N_rad])
    
    Dr_mat  = ns.get_cheb_Dmat(r_co)
    
    Drr_mat = Dr_mat @ Dr_mat

    for idx in range(0, N_rad):
        
        Dr_by_r_mat[idx, :] = Dr_mat[idx, :] / r_co[idx]

    Ublk.Dr_mat  = Dr_mat
    Ublk.Drr_mat = Drr_mat    
    Ublk.Dr_by_r_mat = Dr_by_r_mat

    return Ublk
    
    
def cons_A_mat(Ublk):
    
    N_rad = len(Ublk.r_co)
    
    m_azi = Ublk.m_azi
    k_axi = Ublk.k_axi
    nu    = Ublk.nu
    
    r_co = Ublk.r_co
    uz_mean = Ublk.uz_mean
    
    A_mat =      np.zeros([4 * N_rad, 4 * N_rad]) + \
            1j * np.zeros([4 * N_rad, 4 * N_rad])

    Dr_mat      = Ublk.Dr_mat
    Drr_mat     = Ublk.Drr_mat
    Dr_by_r_mat = Ublk.Dr_by_r_mat

    uz_mean_r = Dr_mat @ r_co
    uz_mean_r = np.asarray(uz_mean_r)
    uz_mean_r = uz_mean_r.flatten()    
       
#   Row 1    
    
    A_mat[1 : (N_rad - 1), 0 : N_rad] = Dr_mat[1 : (N_rad - 1), :]
    
    for idx in range(1, (N_rad - 1)):
    
        A_mat[idx, idx]             = A_mat[idx, idx] + 1 / r_co[idx]
        
        A_mat[idx, N_rad + idx]     = A_mat[idx, N_rad + idx] + \
                                               ((1j * m_azi) / r_co[idx])
        
        A_mat[idx, 2 * N_rad + idx] = A_mat[idx, 2 * N_rad + idx] + \
                                                (1j * k_axi)
        
#   Row 2
    
    A_mat[N_rad + 1 : (2 * N_rad - 1), 0 : N_rad] =  \
                                (1j * nu) *     Drr_mat[1 : (N_rad - 1), :] + \
                                (1j * nu) * Dr_by_r_mat[1 : (N_rad - 1), :]  
    
    A_mat[N_rad + 1 : (2 * N_rad - 1), 3 * N_rad : 4 * N_rad] = -1j * Dr_mat[1 : (N_rad - 1), :]
    
    for idx in range(1, (N_rad - 1)):
        
        A_mat[N_rad + idx, idx] = A_mat[N_rad + idx, idx]                      + \
                                  k_axi * uz_mean[idx]                    - \
                                  ((1j * nu * (m_azi**2))/(r_co[idx]**2)) - \
                                  ((1j * nu) / (r_co[idx]**2))            - \
                                  (1j * nu * (k_axi**2))
                    
        A_mat[N_rad + idx, N_rad + idx] = A_mat[N_rad + idx, N_rad + idx] + \
                                          ((2 * m_azi * nu)/(r_co[idx]**2))
                                          
#   Row 3

    A_mat[(2 * N_rad + 1) : (3 * N_rad - 1), N_rad : (2 * N_rad)] = \
                                (1j * nu) *     Drr_mat[1 : (N_rad - 1), :] + \
                                (1j * nu) * Dr_by_r_mat[1 : (N_rad - 1), :]       
                            
    for idx in range(1, (N_rad - 1)):
        
        A_mat[2 * N_rad + idx, idx] = -2 * nu * m_azi / (r_co[idx]**2)
        
        A_mat[2 * N_rad + idx, N_rad + idx] = \
                                        A_mat[2 * N_rad + idx, N_rad + idx] + \
                                                (k_axi * uz_mean[idx]) - \
                                 ((1j * nu * m_azi**2)/(r_co[idx]**2)) - \
                                 ((1j * nu)           /(r_co[idx]**2)) - \
                                                     (1j * nu * (k_axi**2))
                                  
        A_mat[2 * N_rad + idx, 3 * N_rad + idx] = m_azi / (r_co[idx])
        
        
#   Row 4

    A_mat[(3 * N_rad + 1) : (4 * N_rad - 1), (2 * N_rad) : (3 * N_rad)] = \
                                (1j * nu) *     Drr_mat[1 : (N_rad - 1), :] + \
                                (1j * nu) * Dr_by_r_mat[1 : (N_rad - 1), :]       

    for idx in range(1, (N_rad - 1)):
        
        A_mat[(3 * N_rad + idx), idx] = -1j * uz_mean_r[idx]
        
        A_mat[(3 * N_rad + idx), (2 * N_rad + idx)] = \
                                A_mat[(3 * N_rad + idx), (2 * N_rad + idx)] + \
                                                (k_axi * uz_mean[idx]) - \
                                 ((1j * nu * m_azi**2)/(r_co[idx]**2)) - \
                                                     (1j * nu * (k_axi**2))
                                  
        A_mat[3 * N_rad + idx, 3 * N_rad + idx] = k_axi
    
    Ublk.A_mat = A_mat
    
    return Ublk
        
        
def cons_B_mat(Ublk):
    
    N_rad = len(Ublk.r_co)
    
    B_mat = np.zeros([4 * N_rad, 4 * N_rad]) 
    
    for idx in range(1, (N_rad - 1)):

        B_mat[    N_rad + idx,             idx] = 1        
        B_mat[2 * N_rad + idx,     N_rad + idx] = 1        
        B_mat[3 * N_rad + idx, 2 * N_rad + idx] = 1                

    Ublk.B_mat = B_mat

    return Ublk


def set_incomp_BC(Ublk):
    
    N_rad = len(Ublk.r_co)    
    
    A_mat = Ublk.A_mat
    
    Dr_mat = Ublk.Dr_mat
    
    m_azi  = Ublk.m_azi
    
    # Wall boundary condition     
    A_mat[        0,         0] = 1
    A_mat[    N_rad,     N_rad] = 1    
    A_mat[2 * N_rad, 2 * N_rad] = 1    
    A_mat[3 * N_rad, (3 * N_rad) : (4 * N_rad)] = Dr_mat[0, :]
    
    
    # r=0 boundary condition
    A_mat[(N_rad - 1),  0 : N_rad] = Dr_mat[(N_rad - 1), :]
    
    A_mat[(2 * N_rad - 1), (    N_rad - 1)] = 1
    A_mat[(2 * N_rad - 1), (2 * N_rad - 1)] = 1j * m_azi
    
    A_mat[(3 * N_rad - 1), (3 * N_rad - 1)] = 0
    
    A_mat[(4 * N_rad - 1), (4 * N_rad - 1)] = 0

    Ublk.A_mat = A_mat

    return Ublk

def comp_gen_eig(Ublk):
    
    A_mat = Ublk.A_mat
    B_mat = Ublk.B_mat
    
    EW, EV = scp.linalg.eig(A_mat, B_mat)
    
    Ublk.EW = EW
    Ublk.EV = EV
    
    return Ublk