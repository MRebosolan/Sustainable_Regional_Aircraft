import input 
from CarbonFootprint import cf
from cabindesign import cabin_design

def CLASS1WEIGHTHYBRID(H_to_ker_ratio = input.H_to_ker_ratio,OEWINPUT = 1, top_selecter = 0):
    W_hydrosys=H_to_ker_ratio*1500 #initial guess for hydro system weight
    e=2.71828182846
    n_pax= input.Npax
    W_pax= input.W_pax
    W_cargo= input.W_cargo
    n_crew= input.n_crew
    W_payload=input.W_payload
    Design_range= input.Design_range
    g=9.81
    HYDROGEN_DENSITY= input.rho_hydrogen
    CABIN_LENGTH= input.lpax

    LD_c=input.LD_c
    LD_c2=input.LD_c2
    LD_loiter=input.LD_loiter

    V_c= input.V_C*0.514444
    V_c2= input.V_C2*0.514444
    V_loiter= input.V_loiter*0.514444

    R_c=Design_range-100 #km, correct for take off and landing covered distance

    cj_ck=input.cj_ck
    cj_c=(H_to_ker_ratio-1)*0.349*cj_ck/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio)+0.349*(cj_ck-(H_to_ker_ratio-1)*0.349*cj_ck/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio))
    

    cj_ck2=input.cj_ck2
    cj_c2=(H_to_ker_ratio-1)*0.349*cj_ck2/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio)+0.349*(cj_ck2-(H_to_ker_ratio-1)*0.349*cj_ck2/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio))
    
    cj_kloiter=input.cj_kloiter
    cj_loiter=(H_to_ker_ratio-1)*0.349*cj_kloiter/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio)+0.349*(cj_kloiter-(H_to_ker_ratio-1)*0.349*cj_kloiter/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio))

    np_c=0.82
    np_loiter=0.77

    t_loiter=input.t_loiter#s
    R_loiter=t_loiter*V_loiter

    R_c2=200


    trapped=0.001 #Trapped fuel as fraction of MTOW

    ###FUEL FRACTIONS---ROSKAM---VERIFIED
    end1=0.99+(1-0.99)*H_to_ker_ratio*(1-0.349)
    Mff1=end1

    end2=0.99+(1-0.99)*H_to_ker_ratio*(1-0.349)
    Mff2=Mff1*end2

    end3=0.995+(1-0.995)*H_to_ker_ratio*(1-0.349)
    Mff3=Mff2*end3

    end4=0.98+(1-0.98)*H_to_ker_ratio*(1-0.349)
    Mff4=Mff3*end4

    end5=1/e**(R_c*1000*g*cj_c/LD_c/V_c)
    Mff5=Mff4*end5

    end6=0.99+(1-0.99)*H_to_ker_ratio*(1-0.349)
    Mff6=Mff5*end6

    end7=0.99+(1-0.99)*H_to_ker_ratio*(1-0.349)
    Mff7=Mff6*end7

    end8=1/e**(R_c2*1000*g*cj_c2/LD_c2/V_c2)
    Mff8=Mff7*end8


    end9=1/e**(t_loiter*g*cj_loiter/LD_loiter)
    Mff9=Mff8*end9

    end10=0.99+(1-0.99)*H_to_ker_ratio*(1-0.349)
    Mff10=Mff9*end10

    end11=0.992+(1-0.992)*H_to_ker_ratio*(1-0.349)
    Mff11=Mff10*end11

    FUELFRACMTOW=1-Mff11
    ###REGRESSION DATA OEW=a MTOW+b
    a_reg=0.5753
    b_reg=393.89+W_hydrosys

    MTOW=(b_reg+W_payload)/(1-a_reg-FUELFRACMTOW)

    if OEWINPUT!=1:
        MTOW=(OEWINPUT+W_payload)/(1-FUELFRACMTOW) #no trapped fuel as hydrogen aircraft
    
    FUEL=FUELFRACMTOW*MTOW
    OEW=MTOW-W_payload-FUEL
    MZFW=MTOW-FUEL
    KEROSENE=(1-H_to_ker_ratio)*FUEL
    HYDROGEN=FUEL-KEROSENE
    HYDROGENVOLUME=1.1*1.072*HYDROGEN/HYDROGEN_DENSITY #NASA PAPER
    print(HYDROGENVOLUME)
    
    if HYDROGENVOLUME!=0:
        t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail\
           ,t_top,m_top,tm_top,d_top,l_top,t_pod,m_pod,tm_pod,d_pod,l_pod,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,V_tank_pod,\
           tm_tanksystem,CGtank,CGfuelfull,CGcomb,totdrag,fuselage_weight,CDzerofus,FFbody,Cfturb,fuselage_area,CDzeropods,fusdrag,poddrag=cabin_design(0,0,HYDROGENVOLUME, top_selecter)
    else:
        d_top=0
        tm_tanksystem=0
        


    INFO=[MTOW,OEW,FUEL,W_payload,(MZFW),(KEROSENE),(HYDROGEN),HYDROGENVOLUME,d_top,tm_tanksystem]
    return INFO

