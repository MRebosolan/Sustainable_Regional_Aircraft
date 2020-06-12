"""
input: number of pax
output: oxygen sys weight for class 2
responsible: jorn
"""
def oxygen_system_weight(Npax):
    #Npax=number of passengers
    Wox = 30 + 1.2*Npax
    return Wox
