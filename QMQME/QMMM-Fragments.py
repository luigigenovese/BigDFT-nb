
# coding: utf-8

# # Creation of External Potential from QM multipoles

# In this notebook we will illustrate how a (list of) instance(s) of a Logfile class of the BigDFT module can be used to create a (list of) input files where a QM active region is submitted to an external potential.
# In other terms, from a QM simulation of a system, we will create two distinct regions: the active QM region, made of a subset of the full systems, and the environment region, modelled by an external potential, defined in terms of electrostatic multipoles placed in points which may lie outside the simulation domain of the active QM region.

# Let us start by defining the function which writes the inputfiles after we performed the partition the system into QM and MM regions:

# In[2]:

from futile import Yaml
def append_inputfile(qm,mm,name,tar=None):
    rad=str(name)
    xyzname=rad+'.xyz'
    xyznamemm=rad+'-mm.xyz'
    yamlname=rad+'.yaml'
    out={'dft': {'external_potential': mm.dict()},'posinp': xyzname}
    qm.xyz(xyzname)
    mm.xyz(xyznamemm)
    Yaml.dump(out,yamlname,raw=True,tar=tar)
    return rad


# Let us now define the function that deals with a System instance to provide the two regions:

# In[3]:

def partition_system(system,frags):
    """
    Defines a QM and an MM system from a partition of the original one.
    the parameter frags contains a list of fragment ids that have to be used for the definition of
    the QM region.
    """
    import copy
    from BigDFT import Fragments as F
    Environment=copy.deepcopy(system)
    #let us exclude one fragment from the environment, that will be promoted to be the active QM region
    QM=F.System(units=system.units)
    #reverse the list such as the fragments might be pop-ed without reorder
    toqm=frags
    toqm.sort()
    toqm.reverse()
    for f in toqm:
        QM.append(Environment.pop(f))
    MM=Environment
    return QM,MM


# Let us now define the function which reads a Logfile instance to create a System instance

# In[4]:

from BigDFT import Fragments as F
def fragmentation(full,nat_frag):
    units=full.log['Multipole coefficients']['units']
    Environment=F.System(mp_dict=full.electrostatic_multipoles,units=units,nat_reference=nat_frag)
    return Environment


# And now we need a function which provide a System instance together with the name of the run

# In[5]:

from BigDFT import Logfiles as lf
def QM_snapshot(filename,fragment_size):
    FullQM=lf.Logfile(filename)
    s=filename
    name=s[s.find('log-')+4:].rstrip('.yaml')
    return name,fragmentation(FullQM,fragment_size)


# For a given file, iterate over a number of fragments to create the list of the QM snapshots

# In[6]:

def log_to_QMMM(filename,fragment_size,nfrag=None,tar=None,center=True,qm_list=[]):
    #extract the system out of the logfile
    name,system=QM_snapshot(filename,fragment_size)
    print 'Read run name "',name,'" from file "',filename,'", No. of fragments',len(system.fragments)
    list_posinp=[] #create the list of the run names that have to be performed
    limit=len(system.fragments) if nfrag is None else nfrag
    if center: 
        todo=[[system.central_fragment()]]
    elif len(qm_list)>0 :
        todo=[qm_list]
    else:
        todo=[[i] for i in range(limit)]
    for i in todo:
        #for each chosen fragment partition the system in QM and MM region
        qm,mm=partition_system(system,i)
        print i,name,limit
        # Now we can create an input file which is associated to the corresponding run
        list_posinp.append(append_inputfile(qm,mm,name+'-'+str(i),tar=tar))
    return list_posinp


# Set the arguments that have to be used when using this script from the command line.
# First, determine if we are in a notebook of not:

# In[7]:

try:
    __IPYTHON__
    innb=True
except:
    innb=False


# Then perform the test run or the full run in the script case

# In[8]:

import yaml
if innb:
    #perform a test with the test file
    log_to_QMMM('logs-fullQM/log-snap02000-fullQM.yaml',3,1)
else:        
    import UniParse
    args=UniParse.UniParser('Fragment extraction from homogeneous QM/MM creation')
    args.option('-s','--size',help='fragment size for the QM region',default=3)
    args.option('-f','--files',remainder=True,help='files to be used for the snapshots',default='test')
    args.option('-t','--tarfile',default=None,
                help='archive on which to write the data, will be extended by .tar.bz2')
    args.option('-m','--mode',
                help="""QMMM mode: if "center", 
                only the central fragment will be converted, otherwise all the fragments.""",default='center')
    args.option('-e','--environment',
               help='A list might be passed in order to indicate which fragments have to be promoted as environmental',
               default=[])
    arg=args.args()
    #open the tarfile archive if required to do so
    if arg.tarfile: 
        import tarfile
        tar=tarfile.open(arg.tarfile+'.tar.bz2',mode='w:bz2')
    else:
        tar=None
    list_posinp=[] #store the names of the different runs to be treated in the list 
    center=arg.mode=='center'
    qm_list=[]
    if not center: qm_list=yaml.load(arg.mode)
    for f in arg.files:
        list_posinp+=log_to_QMMM(f,int(arg.size),tar=tar,center=center,qm_list=qm_list)
        print 'File:"',f,'" treated, tot snapshots',len(list_posinp)
    if tar: tar.close()
    Yaml.dump(list_posinp,'list_posinp')


# In[ ]:



