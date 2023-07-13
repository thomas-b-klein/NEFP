from matplotlib import pyplot as plt
import numpy as np
from ifr import ISENTROPIC as ISEN

def plot_nozzle(aeat, prop, sub = False, stag_point = None):
    #converging section
    xc = np.linspace(-3,0,301)
    yc_wall = 1 + (xc**2) / 4.5

    #diverging section
    xd = np.linspace(0,10,1001)
    yd_wall = np.tanh(2*xd/10.0) / np.tanh(2)*(aeat-1)+1

    x = np.concatenate((xc,xd),axis=0)
    y_wall =  np.concatenate((yc_wall,yd_wall))
    y = np.linspace(0,max(y_wall),int(max(y_wall)*100+1))
    XX,YY = np.meshgrid(x,y)
    ZZ = np.ones(np.shape(XX))

    converging = ISEN.from_aas(yc_wall,subsonic=True)
    diverging = ISEN.from_aas(yd_wall,subsonic=sub)
    c_values = getattr(converging, prop)
    d_values = getattr(diverging, prop)
    values = np.concatenate((c_values, d_values),axis=0)
    if stag_point: values = values*stag_point
    for i,lst in enumerate(ZZ):
        ZZ[i] = values

    fig,ax = plt.subplots()
    fig.set_size_inches(15,2.5)
    
    #plot wall
    ax.plot(x,y_wall,color='k')

    #plot contours
    vmin = 0
    vmax = int(np.max(ZZ))
    contourf_ = ax.contourf(XX,YY,ZZ,levels=1000)
    cbar_ = fig.colorbar(contourf_)
    plt.show()

# plot_nozzle(5, "tt0", sub=False, stag_point=1300)




