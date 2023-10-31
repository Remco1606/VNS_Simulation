
from model import  init_model,  calcesI, SD_Curve
v,tvec = init_model()
x = 5000
y = 3500
calcesI(x,y)
DataSD = SD_Curve(v)
import matplotlib.pyplot as plt
plt.plot(DataSD['SD_Duration'],DataSD['SD_Threshold'])
plt.title("Strength Duration Curve")
plt.suptitle("Electrode Location x = " + str(x/1000)+ "mm " + "y = " + str(y/1000)+"mm")
plt.xlabel("Pulse Duration [ms]")
plt.ylabel("Threshold Amplitude [mV]")
plt.show()