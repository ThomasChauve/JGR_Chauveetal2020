### Identification of spheres on discontinuity
import gts
import numpy as np
#import math

# angle of the crack with the horizontal (degrees)
#sc=0.0125
#sp=6
#Rmean=sc*100000.**(1./3)
#alpha_patch=45.


alpha=math.radians(alpha_patch)
radius = sp*sc/2.
discret=50
center=np.array(PP)


# creation des vertices des triangles
pcenter=gts.Vertex(center[0]*X,center[1]*Y,center[2]*Z)
pt=[]

angle_step=2.*math.pi/(discret)
for i in list(range(discret)):
	x=radius*math.cos(i*angle_step)
	y=radius*math.sin(i*angle_step)
	z=0.
	pt.append(gts.Vertex(center[0]*X+x,center[1]*Y+math.cos(alpha)*y,center[2]*Z+math.sin(alpha)*y))
	 

# creation des edges des triangles
edge=[]
for i in list(range(discret)):
	if i==discret-1:
		edge.append([gts.Edge(pcenter,pt[i]),gts.Edge(pt[i],pt[0]),gts.Edge(pt[0],pcenter)])
	else:
		edge.append([gts.Edge(pcenter,pt[i]),gts.Edge(pt[i],pt[i+1]),gts.Edge(pt[i+1],pcenter)])


# creation des faces des triangles
face=[]
for i in list(range(discret)):
	face.append(gts.Face(edge[i][0],edge[i][1],edge[i][2]))


# creation de la surface et importation dans la simulation
s = gts.Surface()
for i in list(range(discret)):
	s.add(face[i])

#print s.area()
#print math.pi*radius**2.

facet = gtsSurface2Facets(s,wire = False,material=wallMat)
O.bodies.append(facet)
