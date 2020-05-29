 #furnishing\
#Verified by hand calculation

def furnishing_gd(N_fdc, N_pax, N_cc, MTOW, P_c):
    K_lav = 0.31 #short range, lavatory
    K_buf = 1.02 #short range
    # P_c = design cabin pressure
    #N_fdc = flight deck crew
    #N_cc = cabin crew
    return 55*N_fdc + 32*N_pax + 15*N_cc + K_lav*(N_pax**1.33) + K_buf *(N_pax**1.12) \
        +109*((N_pax*(1+P_c)/100)**0.505) + 0.771*(MTOW/1000)
        

# def furnishing_torenbeek(MTOW, W_F):
#     return 0.211*((MTOW-W_F)**0.91)
