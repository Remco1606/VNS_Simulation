Running the code
in python:

>>> from model import init_model, Stimulation, calcesI, SD_Curve
>>> v, tvec = init_model()
Init_model complete
>>> calcesI(5000,50)
calcesI complete
>>> DataSD = SD_Curve(v)
>>> import matplotlib.pyplot as plt
>>> plt.plot(DataSD['SD_Duration'],DataSD['SD_Threshold'])
[<matplotlib.lines.Line2D object at 0x000002AD22B56C50>]
>>> plt.show()
