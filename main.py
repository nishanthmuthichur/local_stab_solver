import numpy as np

import matplotlib.pyplot as plt

import ls_utils as lsu

def comp_loc_stab_main(ip_obj):
    
# Construct the derivative matrices
    ip_obj = lsu.cons_deriv_mat(ip_obj)    
    
# Construct the A matrix
    ip_obj = lsu.cons_A_mat(ip_obj)

    print('Local stability: A matrix has been constructed.')
# Construct the B matrix
    ip_obj = lsu.cons_B_mat(ip_obj)
    print('Local stability: B matrix has been constructed.')    

# Setup the BC
    ip_obj = lsu.set_incomp_BC(ip_obj)

# Compute the eigenvalues. 
    ip_obj = lsu.comp_gen_eig(ip_obj)

    EW = ip_obj.EW[-300 : -1]

    fig, ax = plt.subplots(1, 1)
    fig.set_dpi(300)

    ax.plot(np.real(EW), np.imag(EW), 'rx')
    ax.grid()

    plt.show()