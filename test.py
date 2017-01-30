import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import trapz

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

# Simple mouse click function to store coordinates
def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    # print 'x = %d, y = %d'%(
    #     ix, iy)

    # assign global variable to access outside of function
    global coords
    coords.append((ix, iy))

    # Disconnect after 2 clicks
    if len(coords) == 2:
        fig.canvas.mpl_disconnect(cid)
        plt.close(1)
    return


x = np.arange(-10,10)
y = x**2

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(x,y)

coords = []

# Call click func
cid = fig.canvas.mpl_connect('button_press_event', onclick)

plt.show(1)


# limits for integration
ch1 = np.where(x == (find_nearest(x, coords[0][0])))
ch2 = np.where(x == (find_nearest(x, coords[1][0])))

# Calculate integral
y_int = trapz(y[ch1[0][0]:ch2[0][0]], x = x[ch1[0][0]:ch2[0][0]])

print ''
print 'Integral between '+str(coords[0][0])+ ' & ' +str(coords[1][0])
print y_int
