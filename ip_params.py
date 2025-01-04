import numpy as np
import ls_utils as lsu
import main as mn

PI = np.pi

ip_obj = lsu.loc_cls()

ip_obj.N_rad = 201

ip_obj.nu = 1 / 5000
ip_obj.U0 = 1
ip_obj.r0 = 1

ip_obj.m_azi = 1
ip_obj.k_axi = 1

ip_obj.N = ip_obj.N_rad - 1

ip_obj.r_co = np.cos((np.linspace(0, ip_obj.N, ip_obj.N_rad) * PI) / ip_obj.N)

ip_obj.uz_mean = (1 - (ip_obj.r_co)**2)

mn.comp_loc_stab_main(ip_obj)
