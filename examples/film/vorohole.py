import numpy as np
import matplotlib.pyplot as pl
from scipy.spatial import Voronoi,voronoi_plot_2d

# 31 holes
xh=np.array([-48.993637046529095,-43.753369272237194,-43.17579606194892 ,-40.78616352201257 ,-30.849491348578386,-30.7672646156706  ,-29.372566636717583,-26.185054373083002,-23.979591836734695,-21.505414479595217,-20.87235597788844 ,-20.579072953029033,-17.59445383525972 ,-14.69227313566936 ,-6.933962264150935 ,-5.5424528301886795,-3.715184186882304 ,9.016607251543334  ,9.404407244526404  ,10.309973045822098 ,13.217173661917592 ,24.713907410034466 ,25.669362084456395 ,29.49011680143755  ,31.704851752021554 ,32.367475292003604 ,34.21170851657446  ,38.25780072182374  ,45.03732116939662  ,46.606019766397125 ,49.35983827493261])
yh=np.array([-36.4480247623806  ,-34.213807680788804,-27.373786637357966, 15.07214778146384 ,-20.785372063521777,-14.146690844397398, 19.179069688672126, 11.507796245719163, 49.42837193450884 ,-13.975486990941299, 34.93752459964086 ,-29.565907369514704, 11.955438263101364, 15.508605322874203,-26.686039982030536,-2.468913681219931 ,-24.400262342479337,-37.377870848465776,-13.767086211537539, 36.93420207530489 ,-42.92346040577172 , 11.88515518052511 ,-49.24684526689242 , 19.47046991268691 ,-11.1310531711475  , 44.81586145585473 ,-21.357137325049436,-34.61274344858849 ,-22.528444432925525,-46.57185548989793 , 26.800256871011605])
holes=np.transpose(np.array([xh,yh]))

# Periodicity
holes=np.vstack((holes,holes+[+100,0],holes+[-100,0],holes+[0,+100],holes+[0,-100],holes+[+100,+100],holes+[-100,-100],holes+[+100,-100],holes+[-100,+100]))

# Form Voronoi diagram
vor=Voronoi(holes)

# Now plot Voronoi
fig=voronoi_plot_2d(vor,show_vertices=True,line_colors='orange',line_width=2,line_alpha=0.6, point_size=2)
pl.xlabel('$x/h_{film}$',size=16)
pl.ylabel('$y/h_{film}$',size=16)
ticks=np.linspace(-50,+50,5)
pl.xticks(ticks)
pl.yticks(ticks)
pl.xlim(xmin=-50,xmax=+50)
pl.ylim(ymin=-50,ymax=+50)
pl.show()

# Plot holes
#fig,ax=pl.subplots()
#ax.scatter(x=xh,y=yh,marker='o')
#ax.set_xlabel('$x/h_{film}$',size=16)
#ax.set_ylabel('$y/h_{film}$',size=16)
#ax.set_title(r'Initial hole distribution',size=18)
#ax.set_xlim(xmin=-50,xmax=+50)
#ax.set_ylim(ymin=-50,ymax=+50)
#ticks=np.linspace(-50,+50,5)
#ax.set_xticks(ticks)
#ax.set_yticks(ticks)
#pl.rcParams["figure.figsize"]=(2,2)
#fig.tight_layout()
#pl.show()