#%%

def pitch_ratio_fun(D_o, Layout_angle_deg):

    if Layout_angle_deg == 45 or Layout_angle_deg == 0: # Square arrangement 
        if D_o == 1/8:
            Pitch_ratio = 1.25  
        elif D_o == 1/4:
            Pitch_ratio = 1.25  
        elif D_o == 3/8:
            Pitch_ratio = 1.33    
        elif D_o == 1/2:
            Pitch_ratio = 1.25
        elif D_o == 5/8:
            Pitch_ratio = 1.25
        elif D_o == 3/4:
            Pitch_ratio = (1)/D_o
        elif D_o == 1:
            Pitch_ratio = (1+1/4)/D_o
        elif D_o == 1+1/4:
            Pitch_ratio = (1+9/16)/D_o
        elif D_o == 1+1/2:
            Pitch_ratio = (1+7/8)/D_o
        else:
            print("This outer diameter is not considered.")
            return -1

    elif Layout_angle_deg == 30 or Layout_angle_deg == 60: # Square arrangement 
        if D_o == 1/8:
            Pitch_ratio = 1.25  
        elif D_o == 1/4:
            Pitch_ratio = 1.25  
        elif D_o == 3/8:
            Pitch_ratio = 1.33     
        elif D_o == 1/2:
            Pitch_ratio = 1.25
        elif D_o == 5/8:
            Pitch_ratio = 1.25
        elif D_o == 3/4:
            Pitch_ratio = (15/16)/D_o
        elif D_o == 1:
            Pitch_ratio = (1+1/4)/D_o
        elif D_o == 1+1/4:
            Pitch_ratio = (1+9/16)/D_o
        elif D_o == 1+1/2:
            Pitch_ratio = (1+7/8)/D_o
        else:
            print("This outer diameter is not considered.")
            return -1
    else:
        print("This tube arrangement is not considered.")
        return -1

    return Pitch_ratio
