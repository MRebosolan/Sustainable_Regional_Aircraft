import input

##Class I Weight estimation
#def CLASS1WEIGHT(hydro):
#    W_hydrosys=1200
#    e=2.71828182846
#    n_pax=75
#    W_pax=93 #includes luggage
#    W_cargo=1000
#    n_crew=4
#    W_payload=n_pax*W_pax+W_cargo
#    Design_range=2000#[km]
#    g=9.81
#    
#    LD_c=15
#    LD_c2=17
#    LD_loiter=17
#    HYDROGEN_DENSITY=70.8
#    V_c=230.3
#    V_c2=0.8*V_c
#    V_loiter=0.6*V_c
#    
#    R_c=Design_range-100 #km, correct for take off and landing covered distance
#    
#    cj_ck=1.98291*10**(-5) #kerosene cj
#    cj_ch=1.98291*10**(-5)*0.349 #hydrogen cj
#    
#    cj_c2=1.84128*10**(-5)
#    cj_loiter=1.41637*10**(-5)
#    
#    np_c=0.82
#    np_loiter=0.77
#    
#    t_loiter=1800#s
#    R_loiter=t_loiter*V_loiter
#    
#    R_c2=200
#    
#        
#    trapped=0.001 #Trapped fuel as fraction of MTOW
#    
#    ###FUEL FRACTIONS---ROSKAM---VERIFIED
#    end1=0.99+hydro*0.0065
#    Mff1=end1
#    
#    end2=0.99+hydro*0.0065
#    Mff2=Mff1*end2
#    
#    end3=0.995+hydro*0.00325
#    Mff3=Mff2*end3
#    
#    end4=0.98+hydro*0.013
#    Mff4=Mff3*end4
#    
#    end5=1/e**(R_c*1000*g*cj_ch/LD_c/V_c)
#    Mff5=Mff4*end5
#    
#    end6=0.99+hydro*0.0065
#    Mff6=Mff5*end6
#    
#    end7=0.99+hydro*0.0065
#    Mff7=Mff6*end7
#    
#    end8=1/e**(R_c2*1000*g*cj_c2/LD_c2/V_c2)
#    Mff8=Mff7*end8
#    
#    
#    end9=1/e**(t_loiter*g*cj_loiter/LD_loiter)
#    Mff9=Mff8*end9
#    
#    end10=0.99
#    Mff10=Mff9*end10
#    
#    end11=0.992+hydro*0.005
#    Mff11=Mff10*end11
#    
#    FUELFRACMTOW=1-Mff11
#    ###REGRESSION DATA OEW=a MTOW+b
#    a_reg=0.5753
#    b_reg=393.89+W_hydrosys
#    
#    MTOW=(b_reg+W_payload)/(1-a_reg-FUELFRACMTOW)
#    FUEL=FUELFRACMTOW*MTOW
#    OEW=MTOW-W_payload-FUEL
#    MZFW=MTOW-FUEL
#    KEROSENE=MTOW*(Mff8-Mff10)
#    HYDROGEN=FUEL-KEROSENE
#    HYDROGENVOLUME=HYDROGEN/HYDROGEN_DENSITY
#
#    INFO=[int(MTOW),int(OEW),int(FUEL),int(W_payload),int(MZFW),int(KEROSENE),int(HYDROGEN),int(HYDROGENVOLUME)]
#    return(INFO)
    
    

def CLASS1WEIGHTHYBRID(H_to_ker_ratio = input.H_to_ker_ratio,OEWINPUT = 1):
    W_hydrosys=H_to_ker_ratio*1200
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

    V_c= input.V_c
    V_c2= input.V_c2
    V_loiter= input.V_loiter
    
    R_c=Design_range-100 #km, correct for take off and landing covered distance
    
    cj_ck=input.cj_ck
    cj_c = cj_ck * 0.349 * H_to_ker_ratio + cj_ck * (1 - H_to_ker_ratio)
    
    cj_ck2=input.cj_ck2
    cj_c2 = cj_ck2 * 0.349 * H_to_ker_ratio + cj_ck2 * (1 - H_to_ker_ratio)

    cj_kloiter=input.cj_kloiter
    cj_loiter = cj_kloiter * 0.349 * H_to_ker_ratio + cj_kloiter * (1 - H_to_ker_ratio)


    np_c=0.82
    np_loiter=0.77
    
    t_loiter=1800#s
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
        MTOW=(OEWINPUT+W_payload)/(1-FUELFRACMTOW)
        
    FUEL=FUELFRACMTOW*MTOW
    OEW=MTOW-W_payload-FUEL
    MZFW=MTOW-FUEL
    KEROSENE=(1-H_to_ker_ratio)*FUEL
    HYDROGEN=FUEL-KEROSENE
    HYDROGENVOLUME=1.072*HYDROGEN/HYDROGEN_DENSITY
    
    ###TANK
    R=0
    result=0
    while result<HYDROGENVOLUME:
        R+=0.001
        result=4/3*3.14159*R**3+3.14159*CABIN_LENGTH*R**2
        
    TANK_DIAMETER=R*2
    TANK_SURFACE_AREA = CABIN_LENGTH*3.14159*TANK_DIAMETER + 3.14159*TANK_DIAMETER**2
    TANK_MATERIAL_DENSITY = 2825 #MONOLITHIC METAL Aluminium alloy 2219 KG/M3
    TANK_THICKNESS = 3 # thickness in mm
    STRUCTURAL_TANK_MASS = TANK_SURFACE_AREA*TANK_MATERIAL_DENSITY*TANK_THICKNESS*10**(-3) #tank mass exluding insulation + other systems required
    
    ###TANK
    INFO=[MTOW,OEW,FUEL,W_payload,(MZFW),(KEROSENE),(HYDROGEN),HYDROGENVOLUME,TANK_DIAMETER,STRUCTURAL_TANK_MASS]
    return(INFO)    


mtowlist=[]
xlist=[]
oewlist=[]
kerosenelist=[]
hydrogenlist=[]
tfuellist=[]
tmasslist=[]
hvollist=[]
tdiameterlist=[]

import matplotlib.pyplot as plt

for i in range(0,101):
    outputc1h=CLASS1WEIGHTHYBRID(i/100,1)
    mtowlist.append(outputc1h[0])
    oewlist.append(outputc1h[1])
    kerosenelist.append(outputc1h[5])
    hydrogenlist.append(outputc1h[6])
    tfuellist.append(outputc1h[2])
    hvollist.append(outputc1h[7])
    tdiameterlist.append(outputc1h[8])
    tmasslist.append(outputc1h[9])
    xlist.append(i)

plt.subplot(2,3,1)
plt.plot(xlist,mtowlist,label='MTOW')  
plt.plot(xlist,oewlist,label='OEW')
plt.plot([0,100],[28992,28992],label='MTOW kerosene reserve frac')
plt.plot([0,100],[18273,18273],label='OEW kerosene reserve frac')

plt.ylabel('WEIGHT [kg]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(2,3,2)
plt.plot(xlist,kerosenelist,label='kerosene mass')
plt.plot(xlist,hydrogenlist,label='hydrogen mass')
plt.plot(xlist,tfuellist,label='total fuel mass')
plt.plot([0,100],[660,660],label='kerosene in kerosene reserve frac')
plt.plot([0,100],[2084,2084],label='hydrogen in kerosene reserve frac')
plt.plot([0,100],[2744,2744],label='total fuel in kerosene reserve frac')

plt.ylabel('WEIGHT [kg]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()
plt.subplot(2,3,3)
plt.plot(xlist,hvollist,label='Hydrogen volume')
plt.ylabel('Volume [CUBIC METERS]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(2,3,4)
plt.plot(xlist,tdiameterlist,label='Tank diameter')
plt.ylabel('Meters')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()
plt.subplot(2,3,5)
plt.plot(xlist,tmasslist,label='Tank mass')
plt.ylabel('kgs')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()
plt.show() 

