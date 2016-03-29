import quantumoptics as qo
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set(style='whitegrid', context='paper')

om_ds =  np.linspace(10, 10.8, 5) 
cam = qo.JaynesCummingsSystem(4,
                              10,
                              10,
                              om_ds,
                              [1, 0.2],
                              10,
                              40)

cam.draw_qps(plottype='c', 
             colormap='viridis',
             save=True,
             fontdict={'size':16,
                       'horizontalalignment':'left'})
plt.show()
