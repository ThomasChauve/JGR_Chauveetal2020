# -*- coding: utf-8 -*-
# encoding: utf-8
from yade import ymport, utils, plot,export
import math
import numpy as np
#from pylab import *
import warnings
warnings.filterwarnings("ignore")

#### SIMULATIONS DEFINED HERE
utils.readParamsFromTable(fPACK="Packing/",PACK='111_5k',Y1=4.30e+10,A1=0.01,T1=8e6,C1=30e6,SXX=-11e6,SYY=-30e6,SZZ=-11e6,sp=0.1,alpha_patch=90.,posP=[[0.5,0.5,0.35],[0.5,0.5,0.65],[0.5,0.214285,0.35],[0.5,0.214285,0.65],[0.5,0.785714,0.35],[0.5,0.785714,0.65]],sc=0.01,injectionf=1e-7,timeMax=0.000001,FOLDER='Test/',NAME='-test',factor=1,noTableOk=True)
from yade.params.table import *

#### parameters of the simulation (geometry, boundary conditions, output, etc...)
os.mkdir(FOLDER)

### Material microproperties (set of parameters for simulating Colton sandstone with K=10)
intR=1.06
DENS=4000
YOUNG=Y1
FRICT=2
RFRICT=0
### Mechanical loading (state of stress before the injection)
Sxx=SXX # Sigmaxx
Syy=factor*SYY # Sigmayy
Szz=factor*SZZ # Sigmazz

### Fluid properties
KFluid=2.2e11 # bulk modulus of the fluid (1/compressibility)
visc=1 # viscosity of the fluid
pFactor=1.8e-11 # to scale the permeability of the rock matrix: useless if lines 133-136 are not commented (impermeable matrix) -> cf. permeametre.py: 1.8e-11 gives a permeability of 1e-16 m2 for 111_10k
slotAperture=1e-3 # initial aperture of pre-existing fracture where the injection is done

### hydraulic loading
#flowRate=injectionf # injection flow rate

### Simulation Control
saveData=10 # data record interval
#iterMax=20e3
saveVTK=20 # number of Vtk files
OUT=PACK+NAME

### preprocessing to get dimensions
O.bodies.append(ymport.text(fPACK+PACK+'.spheres',scale=sc))

dim=utils.aabbExtrema()
xinf=dim[0][0]
xsup=dim[1][0]
X=xsup-xinf
yinf=dim[0][1]
ysup=dim[1][1]
Y=ysup-yinf
zinf=dim[0][2]
zsup=dim[1][2]
Z=zsup-zinf

R=0
Rmax=0
numSpheres=0.
for o in O.bodies:
 if isinstance(o.shape,Sphere):
   numSpheres+=1
   R+=o.shape.radius
   if o.shape.radius>Rmax:
     Rmax=o.shape.radius
Rmean=R/numSpheres

print 'X=',X,' | Y=',Y,' | Z=',Z,' || nbSpheres=',numSpheres,' | Rmean=',Rmean

###
O.reset() # all previous lines were for getting dimensions of the packing to create walls at the right positions (below) because walls have to be generated after spheres for FlowEngine
###

#### here we reconstruct the scene with right dimensions (because walls have to be imported before spheres for flow engine)

### material definition
def sphereMat1(): return JCFpmMat(type=1,density=DENS,young=Y1,poisson=A1,tensileStrength=T1,cohesion=C1,frictionAngle=radians(FRICT),jointNormalStiffness=Y1/(pi*Rmean),jointShearStiffness=A1*Y1/(pi*Rmean),jointTensileStrength=0.,jointCohesion=0.,jointFrictionAngle=radians(FRICT),jointDilationAngle=radians(0))


def wallMat(): return JCFpmMat(type=0,density=DENS,young=YOUNG,frictionAngle=radians(0))

### walls
mn,mx=Vector3(xinf+0.1*Rmean,yinf+0.1*Rmean,zinf+0.1*Rmean),Vector3(xsup-0.1*Rmean,ysup-0.1*Rmean,zsup-0.1*Rmean)

walls=utils.aabbWalls(oversizeFactor=1.5,extrema=(mn,mx),thickness=0.1*min(X,Y,Z),material=wallMat)
wallIds=O.bodies.append(walls)

### packing
O.bodies.append(ymport.text(fPACK+PACK+'.spheres',scale=sc,material=sphereMat1))

### DFN
for i in list(range(len(posP))):
	PP=posP[i]
	execfile('patch.py')

execfile('identifyInitialFractures.py')


#### engines
### Triaxial Engine
triax=TriaxialStressController(
	internalCompaction=False
	,stressMask=7
	,goal1=Sxx
	,goal2=Syy
	,goal3=Szz
	,max_vel=0.01
)

### Flow Engine
flow=DFNFlowEngine(	
        isActivated=False
        ,useSolver=3
        ,boundaryUseMaxMin = [0,0,0,0,0,0]
        ,bndCondIsPressure = [1,1,1,1,1,1]
        ,bndCondValue=[0,0,0,0,0,0]
        ,permeabilityFactor=pFactor
        ,viscosity=visc
        ,fluidBulkModulus=KFluid
        ### DFN related
        ,clampKValues=False
        ,jointsResidualAperture=slotAperture
        
        ## segfault with following
        #,updatePositions=True
        
        ## segfault with following
        #,defTolerance=-1
        #,updatePositions=True
        #,meshUpdateInterval=1000000000
        
)

### with DFNFlow, we can block every cells not concerned with fractures with the following function: if these lines are commented (permeable matrix), you will get warning about cholmod: is it an issue? I am not sure yet but it is annoying...
def blockStuff():
	for k in range(flow.nCells()): flow.blockCell(k,True)
flow.blockHook="blockStuff()"

### Simulation is defined here (DEM loop, interaction law, servo control, recording, etc...)
O.engines=[
        ForceResetter(),
        InsertionSortCollider([Bo1_Box_Aabb(),Bo1_Sphere_Aabb(aabbEnlargeFactor=intR,label='Saabb')]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=intR,label='SSgeom'),Ig2_Box_Sphere_ScGeom()],
		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1,label='interactionPhys')],
		[Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(smoothJoint=True,neverErase=1,recordCracks=True,Key=OUT,label='interactionLaw')]
	),
        GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=10,timestepSafetyCoefficient=0.8,defaultDt=0.1*utils.PWaveTimeStep()),
        triax,
        flow,
        NewtonIntegrator(damping=0.4,label="newton"),
        PyRunner(iterPeriod=int(1),initRun=True,command='crackCheck()',label='check'),
        PyRunner(iterPeriod=int(saveData),initRun=True,command='recorder()',label='recData'),
        PyRunner(iterPeriod=int(1),initRun=True,command='saveFlowVTK()',label='saveFlow',dead=1),
        #PyRunner(iterPeriod=int(1),initRun=True,command='saveAperture()',label='saveAperture',dead=1),
        VTKRecorder(iterPeriod=int(1),initRun=True,fileName=OUT+'-',recorders=['spheres','bstresses','cracks'],Key=OUT,label='saveSolid',dead=1)
]
 
# these lines can be a problem depending on the configuration of your computer
#from yade import qt
#v=qt.Controller()
#v=qt.View()

#### custom functions
os.chdir(FOLDER)
#### check if new cracks are created to update "flow mesh permeability"
cks=cks0=0
def crackCheck():
  global tensCks, shearCks, cks, cks0
  cks=interactionLaw.nbTensCracks+interactionLaw.nbShearCracks
  if cks>(cks0): 
    #print 'new crack! Update triangulation!'
    flow.updateTriangulation=True
  cks0=cks

### save flow field (pressure and velocity)
def saveFlowVTK():
 flow.saveVtk(folder='VTK')
 
#### save cracks aperture
#from yade import export
#vtkExporter = export.VTKExporter('cracks')
#def saveAperture():
  #vtkExporter.exportContactPoints(what=[('b','i.phys.isBroken'),('n','i.geom.normal'),('s','i.phys.crossSection'),('a','i.phys.crackJointAperture')])

### save macroscopic data
ex0=ey0=ez0=0
def recorder():
    global ex0,ey0,ez0
    crackVolume=crackSurface=0
    for i in O.interactions:
        if i.phys.isBroken:
	    crackVolume+=i.phys.crossSection*i.phys.crackJointAperture
	    crackSurface+=i.phys.crossSection
    yade.plot.addData( t=O.time
			,i=O.iter
			,ex=triax.strain[0]-ex0
			,ey=triax.strain[1]-ey0
			,ez=triax.strain[2]-ez0
			,sx=triax.stress(triax.wall_right_id)[0]
			,sy=triax.stress(triax.wall_top_id)[1]
			,sz=triax.stress(triax.wall_front_id)[2]
			,p=flow.getPorePressure((xinf+X/2.,yinf+Y/2.,zinf+Z/2.))
			,tc=interactionLaw.nbTensCracks
			,sc=interactionLaw.nbShearCracks
			,p32=crackSurface
			,p33=crackVolume
			,unbF=utils.unbalancedForce() 
			,momEnergy=kineticEnergy
    )
    plot.saveDataTxt(OUT)

def microstresses():
	s=utils.bodyStressTensors()
	TW=TesselationWrapper()
	TW.triangulate()
	TW.computeVolumes()
	for b in O.bodies:
		if type(b.shape)==Sphere:
			b.microstresses = s[b.id]*4.*pi/3.*b.shape.radius**3/TW.volume(b.id)
			b.TWvolume=TW.volume(b.id)
			b.cc=b.shape.color[0]+0.5*b.shape.color[1]+2.*b.shape.color[2]
	

	vtkExporter = export.VTKExporter(OUT+'_microstresses_'+str(O.iter))
	vtkExporter.exportSpheres(what=[('TWvolume','b.TWvolume'),('microstresses','b.microstresses'),('young','b.mat.young'),('color','b.cc')])


O.engines=O.engines+[PyRunner(iterPeriod=int(1),initRun=True,command='microstresses()',label='savemicrostresses',dead=1)]


#### Simulation starts here

### manage interaction detection factor during the first timestep (near neighbour bonds are created at first timestep)
O.step()
## initializes the interaction detection factor to default value (new contacts, frictional, between strictly contacting particles)
SSgeom.interactionDetectionFactor=-1.
Saabb.aabbEnlargeFactor=-1.
saveSolid.dead=1

### mechanical loading
while 1:
  O.run(100, True)
  print 'unbalanced force=',unbalancedForce()
  if ( unbalancedForce()<0.005 and ((abs(abs(triax.stress(triax.wall_right_id)[0])-abs(Sxx))/abs(Sxx))<0.001) and ((abs(abs(triax.stress(triax.wall_top_id)[1])-abs(Syy))/abs(Syy))<0.001) and ((abs(abs(triax.stress(triax.wall_front_id)[2])-abs(Szz))/abs(Szz))<0.001) ):
    print 'stabilizing || iteration=', O.iter
    O.run(100,True) # to further stabilize the system
    print 'confined state || Sxx=',triax.stress(triax.wall_right_id)[0],' | Syy=',triax.stress(triax.wall_top_id)[1],' | Szz=',triax.stress(triax.wall_front_id)[2]
    ex0=triax.strain[0]
    ey0=triax.strain[1]
    ez0=triax.strain[2]
    O.save(OUT+'_confined.yade')
    break
    
### hydraulic loading
#print 'activate flow engine now || iteration=', O.iter
#triax.max_vel=1
#flow.isActivated=1

saveFlow.dead=0
saveSolid.dead=0
savemicrostresses.dead=0
#saveAperture.dead=1
O.step() # needed to avoid segfault?
saveFlow.dead=1
saveSolid.dead=1
savemicrostresses.dead=1

for i in list(range(10)):
    triax.goal2=Syy-((i+1)*10**6)
    print('New goal='+str(triax.goal2))
    while 1:
        O.run(100, True)
        print 'unbalanced force=',unbalancedForce()
        if ( unbalancedForce()<0.005 and ((abs(abs(triax.stress(triax.wall_right_id)[0])-abs(Sxx))/abs(Sxx))<0.001) and ((abs(abs(triax.stress(triax.wall_top_id)[1])-abs(triax.goal2))/abs(triax.goal2))<0.001) and ((abs(abs(triax.stress(triax.wall_front_id)[2])-abs(Szz))/abs(Szz))<0.001) ):
            print 'stabilizing || iteration=', O.iter
            O.run(100,True) # to further stabilize the system
            print 'confined state || Sxx=',triax.stress(triax.wall_right_id)[0],' | Syy=',triax.stress(triax.wall_top_id)[1],' | Szz=',triax.stress(triax.wall_front_id)[2]
            #ex0=triax.strain[0]
            #ey0=triax.strain[1]
            #ez0=triax.strain[2]
            #O.save(OUT+'_confined.yade')
            break

    saveFlow.dead=0
    saveSolid.dead=0
    savemicrostresses.dead=0
    #saveAperture.dead=1
    O.step() # needed to avoid segfault?
    saveFlow.dead=1
    saveSolid.dead=1
    savemicrostresses.dead=1

