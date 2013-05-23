from msmbuilder import Trajectory, metrics
from msmbuilder.geometry import dihedral as _dihedralcalc
import numpy


# function to get mass weights based on atom name
def get_weights(names):
    atomweights=dict()
    atomweights['C']=12
    atomweights['O']=16
    atomweights['N']=14
    atomweights['H']=1
    atomweights['S']=32
    weights=[]
    for x in names:
        try:
            int(x[0])
            weights.append(atomweights[x[1]])
        except ValueError:
            weights.append(atomweights[x[0]])
    return numpy.array(weights)


# load in trajectory, here I have just a pdb with one frame
# you can loop over frames in a trajectory
conf=Trajectory.load_from_pdb('new-aln-state847.cent.pdb')

# here i organize the groups in dictionaries
coms=dict()
groups=dict()
group_wts=dict()
residues=dict()
numbers=[1,2,3,4]
start=1
for num in numbers:
    groups[num]=dict() # group ids
    group_wts[num]=dict() # group mass weight ids
    end=start+1
    residues[num]=[start,end] # adding consecutive resids to groups in this example
    start=end+1

# compute the COM for each group of 2 residues
for group in groups.keys():
    start=residues[group][0]
    for resid in residues[group]:
        array=numpy.where(conf['ResidueID']==resid)[0]
        types=conf['AtomNames'][array]
        if resid==start:
            groups[group]=conf['AtomID'][array]
            group_wts[group]=get_weights(types)
        else:
            groups[group]=numpy.hstack((groups[group], conf['AtomID'][array]))
            group_wts[group]=numpy.hstack((group_wts[group], get_weights(types)))
    indices=[i-1 for i in groups[group]]
    x=[coor[0] for coor in conf['XYZList'][0][indices]]
    y=[coor[1] for coor in conf['XYZList'][0][indices]]
    z=[coor[2] for coor in conf['XYZList'][0][indices]] # here just using one frame
    coors=numpy.vstack((x,y,z))
    coms[group]=numpy.array([numpy.average(numpy.array(x), weights=group_wts[group]), numpy.average(numpy.array(y), weights=group_wts[group]), numpy.average(numpy.array(z), weights=group_wts[group])])


# msmbuilder _dihedralcalc needs a trajectory-like object of coordinates, so stack the
# com coordinates
com_traj=numpy.zeros((1,4, 3)) # here the size is 1 because one frame, 4 groups for the dihedral
indices=numpy.zeros((1,4)) # here i just need to tell _dihedralcalc that we have one set of 4 coordinates
for group in groups.keys():
    index=group-1
    com_traj[0][index]=coms[group]
    indices[0][index]=index

pseudodihed=_dihedralcalc.compute_dihedrals(com_traj, indices, degrees=True)
print residues
print pseudodihed

