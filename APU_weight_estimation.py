''' 
to calculate apu weight based on take off weight
outputs apu weight
responsible: matteo'''
def APU_weight_estimation(Wto):
    #Wto = take-off weight
    Wapu = 0.0075 * Wto
    return Wapu