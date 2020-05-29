#Input in this function is H to ker ratio, which is 1 (put in as integer). OEWINPUT is 1 for first iteration, anything else is for further iterations.
#Output is INFO=[MTOW,OEW,FUEL,W_payload,(MZFW),(KEROSENE),(HYDROGEN),HYDROGENVOLUME,TANK_DIAMETER,TOTAL_STRUCTURAL_TANK_MASS]
#Do not change anything in this function please.
import input 
from CarbonFootprint import cf
from hydrogen_tank_sizing import tank_sizing

def CLASS1WEIGHTHYBRID(H_to_ker_ratio = input.H_to_ker_ratio,OEWINPUT = 1):
    W_hydrosys=H_to_ker_ratio*2000
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

    V_c= input.V_C
    V_c2= input.V_C2
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
        MTOW=(OEWINPUT+W_payload)/(1-FUELFRACMTOW)
    
    FUEL=FUELFRACMTOW*MTOW
    OEW=MTOW-W_payload-FUEL
    MZFW=MTOW-FUEL
    KEROSENE=(1-H_to_ker_ratio)*FUEL
    HYDROGEN=FUEL-KEROSENE
    HYDROGENVOLUME=1.1*1.072*HYDROGEN/HYDROGEN_DENSITY #NASA PAPER
    TANK_THICKNESS,STRUCTURAL_TANK_MASS, TOTAL_STRUCTURAL_TANK_MASS, TANK_DIAMETER=tank_sizing(HYDROGENVOLUME+0.01,CABIN_LENGTH,1)
    


    INFO=[MTOW,OEW,FUEL,W_payload,(MZFW),(KEROSENE),(HYDROGEN),HYDROGENVOLUME,TANK_DIAMETER,TOTAL_STRUCTURAL_TANK_MASS]
    return INFO


mtowlist=[]
xlist=[]
oewlist=[]
kerosenelist=[]
hydrogenlist=[]
tfuellist=[]
tmasslist=[]
hvollist=[]
tdiameterlist=[]
energylist=[]
cjclist=[]
emissionslist=[]
emissionsratiolist=[]


import matplotlib.pyplot as plt

for i in range(99,101):
    outputc1h=CLASS1WEIGHTHYBRID(i/100,1)
    mtowlist.append(outputc1h[0])
    oewlist.append(outputc1h[1])
    kerosenelist.append(outputc1h[5])
    hydrogenlist.append(outputc1h[6])
    tfuellist.append(outputc1h[2])
    hvollist.append(outputc1h[7])
    tdiameterlist.append(outputc1h[8])
    tmasslist.append(outputc1h[9])
    energylist.append(kerosenelist[-1]*42.8+hydrogenlist[-1]*122.8)
    emissionslist.append(cf(tfuellist[-1], i/100, 1-i/100, input.NOx_H2, input.GWP)[0])
    emissionsratiolist.append(cf(tfuellist[-1], i/100, 1-i/100, input.NOx_H2, input.GWP)[1])
    xlist.append(i)


plt.subplot(3,3,1)
plt.plot(xlist,mtowlist,label='MTOW')
plt.plot(xlist,oewlist,label='OEW')
#plt.plot([0,100],[28992,28992],label='MTOW kerosene reserve frac')
#plt.plot([0,100],[18273,18273],label='OEW kerosene reserve frac')

plt.ylabel('WEIGHT [kg]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(3,3,2)
plt.plot(xlist,kerosenelist,label='kerosene mass')
plt.plot(xlist,hydrogenlist,label='hydrogen mass')
plt.plot(xlist,tfuellist,label='total fuel mass')
#plt.plot([0,100],[660,660],label='kerosene in kerosene reserve frac')
#plt.plot([0,100],[2084,2084],label='hydrogen in kerosene reserve frac')
#plt.plot([0,100],[2744,2744],label='total fuel in kerosene reserve frac')

plt.ylabel('WEIGHT [kg]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()
plt.subplot(3,3,3)
plt.plot(xlist,hvollist,label='Hydrogen volume')
plt.ylabel('Volume [CUBIC METERS]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(3,3,4)
plt.plot(xlist,tdiameterlist,label='Tank diameter')
plt.ylabel('Meters')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(3,3,5)
plt.plot(xlist,tmasslist,label='Tank mass')
plt.ylabel('kgs')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(3,3,6)
plt.plot(xlist,energylist,label='total E carried')
plt.plot(xlist, [i*122.8 for i in hydrogenlist],label='E hydrogen')
plt.plot(xlist, [i*42.8 for i in kerosenelist],label='E kerosene')
plt.ylabel('MJ')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(3,3,7)
plt.plot(xlist,[i*input.hydrogen_cost for i in hydrogenlist],label='Hydrogen cost')
plt.plot(xlist,[i*0.6/0.81 for i in kerosenelist],label='kerosene cost')
plt.plot(xlist,[sum(x) for x in zip([i*input.hydrogen_cost for i in hydrogenlist], [i*0.6/0.81 for i in kerosenelist])],label='total cost')
plt.legend()

plt.ylabel('US DOLLARS')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')

plt.subplot(3,3,8)
plt.plot(xlist,emissionslist,label='emissions')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(3,3,9)
plt.plot(xlist,emissionsratiolist,label='CF RATIO wrt CRJ700')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()
plt.show()

