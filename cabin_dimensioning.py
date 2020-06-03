# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 22:16:55 2020

@author: malfl
"""

#FUSELAGE SIZING

#ROSKAM

npax=75

w_seat=18.5*0.0254
w_armrest=1.9*0.0254
w_aisle=0.49
w_clearance=0.4*0.0254


h_headroom=1.27 #m
h_aisle=1.9 #m

three=3*w_seat+4*w_armrest+w_clearance
two=2*w_seat+3*w_armrest+w_clearance

totalwidth=two+three+w_aisle
seat_pitch=32*0.0254
aft_galley=1.8
Vol=30
R=1.500
fuel_length=Vol/3.14159/R**2
cabin_length=seat_pitch*15+aft_galley+fuel_length

print(three,two,totalwidth)
print(fuel_length,aft_galley,cabin_length)