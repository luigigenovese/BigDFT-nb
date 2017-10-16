import extract_fragments as ef, numpy as np, Logfiles as lf

#load the system with PEN as reference fragment
#first read the logfile multipooles
log=lf.document_quantities(lf.get_log('OUT-big-rigid.yaml'),{"mp":['Multipole coefficients']})
#mp=lf.document_quantities(lf.get_log('log-one-1-11.yaml'),{"mp":['Multipole coefficients']})
#then load the system
system=ef.System(log['mp']['values'],nat_reference=36,units='A')

#also load the two pentacenes which have to be taken as references for the system
two=ef.System(xyz='two.xyz',nat_reference=36,units='A')

#then decompose the system into fragments, taking the two fragments as references
system.decompose(two.fragments)

#this prints out the nunmber of atoms of the respective systems
print len(two),len(system) 

#given this decomposition we then might inspect the values of the atomic dipole and monopoles
#system.xyz()
cxyz=system.centroid()
limit=20.0
#create a list of fragments which are internal and that correspond to each of the reference fragments
flst=[[],[]]
for i,f in enumerate(system.fragments):
    print i,np.linalg.norm(system.CMs[i]-cxyz),system.decomposition[i]['J']
    if np.linalg.norm(system.CMs[i]-cxyz) > limit: continue
    ival=system.decomposition[i]['id']
    flst[ival].append(f)

print len(flst[0])
print len(flst[1])

#build the average amond these fragments for the multipoles
ref=[]
for ival,f in enumerate(two.fragments):
    ref.append(ef.frag_average(f,flst[ival]))

print ref[0].d0(),ref[0].Q()
print ref[1].d0(),ref[1].Q()
#recompose the system from these two averaged fragments such as to eliminate the surface effects
MMsystem=ef.System(transformations=system.decomposition,units='A',reference_fragments=ref)

#then build the QM system as some internal fragments starting from the recomposed system
oxyz=MMsystem.CMs[MMsystem.central_fragment()]

#create the list of the QM fragments
QMlist=[]
limit=11.0
for i,f in enumerate(MMsystem.fragments):
    dist=np.linalg.norm(MMsystem.CMs[i]-oxyz)
    #print len(f),f.Q(),f.d0(),np.linalg.norm(f.d1()),dist,MMsystem.CMs[i][2]
    if dist*ef.AU_to_A < limit: QMlist.append(i)

QMlist.sort()
QMlist.reverse()
print len(MMsystem)
print QMlist

#now we should write the input file of BigDFT. Construct a system with the QM list fragments and a 
#external potential region with the exterior system pulled in
QMsystem=ef.System(units='A')
for ifrag in QMlist:
    QMsystem.append(MMsystem.pop(ifrag))

#then print the two systems
print len(QMsystem)
print len(MMsystem)

#suffix of the filename
sx=str(int(limit))

#print the xyz files of the QM and the MM system for cross-check
QMsystem.xyz('QM'+sx+'.xyz',units='angstroem')
MMsystem.xyz('MM'+sx+'.xyz',units='angstroem')
#also print the positions of the center of mass o the MM fragments and  the corresponding dipole norm
cents=ef.XYZfile('MM'+sx+'centroids.xyz')#,units='angstroem')
cents.append(MMsystem.CMs,basename='Cen')
cents.dump()

cents=ef.XYZfile('MM'+sx+'d1.dat')#,units='angstroem')
cents.append([[np.linalg.norm(f.d1())] for f in MMsystem.fragments],basename='d1')
cents.dump()

import yaml
f=open('MM'+sx+'.yaml','w')
f.write(yaml.dump(MMsystem.dict()))
f.close()
