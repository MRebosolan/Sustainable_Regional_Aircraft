#Class I Weight estimation


n_pax=75
W_pax=93 #includes luggage
W_cargo=1000
n_crew=4
Design_range=2000#[km]

LD_c=15
LD_c2=17
LD_loiter=17
g=9.81
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
R_c2=463

###FUEL FRACTIONS###
end1=0.99
end2=0.99
end3=0.995
end4=0.98
end5=
end6=0.99
end7=0.99
end8=
end9=
end10=0.99
end11=0.992

