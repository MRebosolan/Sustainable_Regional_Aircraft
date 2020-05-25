#Class I Weight estimation
def CLASS1WEIGHT(hydro):
    e=2.71828182846
    n_pax=75
    W_pax=93 #includes luggage
    W_cargo=1000
    n_crew=4
    W_payload=n_pax*W_pax+W_cargo
    Design_range=2000#[km]
    g=9.81
    
    LD_c=15
    LD_c2=17
    LD_loiter=17
    
    V_c=230.3
    V_c2=0.8*V_c
    V_loiter=0.6*V_c
    
    R_c=Design_range-100 #km
    
    cj_c=1.98291*10**(-5)
    cj_c2=1.84128*10**(-5)
    cj_loiter=1.41637*10**(-5)
    
    np_c=0.82
    np_loiter=0.77
    
    t_loiter=1800#s
    R_loiter=t_loiter*V_loiter
    
    R_c2=463/2.5
    if hydro:
        cj_c=cj_c/3
        
    trapped=0.001 #Trapped fuel as fraction of MTOW
    
    ###FUEL FRACTIONS---ROSKAM---VERIFIED
    end1=0.99+hydro*0.005
    Mff1=end1
    
    end2=0.99+hydro*0.005
    Mff2=Mff1*end2
    
    end3=0.995+hydro*0.0025
    Mff3=Mff2*end3
    
    end4=0.98+hydro*0.01
    Mff4=Mff3*end4
    
    end5=1/e**(R_c*1000*g*cj_c/LD_c/V_c)
    Mff5=Mff4*end5
    
    end6=0.99+hydro*0.005
    Mff6=Mff5*end6
    
    end7=0.99+hydro*0.005
    Mff7=Mff6*end7
    
    end8=1/e**(R_c2*1000*g*cj_c2/LD_c2/V_c2)
    Mff8=Mff7*end8
    
    end9=1/e**(t_loiter*g*cj_loiter/LD_loiter)
    Mff9=Mff8*end9
    
    end10=0.99
    Mff10=Mff9*end10
    
    end11=0.992+hydro*0.004
    Mff11=Mff10*end11
    
    FUELFRACMTOW=1-Mff11
    ###REGRESSION DATA OEW=a MTOW+b
    a_reg=0.5753
    b_reg=393.89
    
    MTOW=(b_reg+W_payload)/(1-a_reg-FUELFRACMTOW)
    FUEL=FUELFRACMTOW*MTOW
    OEW=MTOW-W_payload-FUEL
    MZFW=MTOW-FUEL
    return(MTOW,OEW,FUEL,W_payload,MZFW)
