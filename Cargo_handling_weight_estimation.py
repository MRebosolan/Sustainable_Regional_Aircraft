''' 
to calculate apu weight based on cargo floor area in class 2
outputs cargo handling equipment weight
responsible: matteo
'''

def cargo_handling_weight(Sff):
    #Sff= freight floor area in square ft
    Wch=3*Sff
    return Wch