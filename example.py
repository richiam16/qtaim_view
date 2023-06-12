from QTAIM_p_3 import main

bencene = main("examples/C6H6_6dp.wfn")

# bencene.graph(0,"xy",-4,4,1000,10e-5,0,200) # A pop out window will show with a surface graph of the electron density (We choose bencene as a testing molecule for being a planar molecule)
# bencene.graph(0,"xy",-4,4,1000,10e-5,1,200) # A pop out window will show with a isodensities graph of the electron density
 # Note: the method object.graph() can graph both surface plots and isodensities plot. it is explanatory to talk about the meaning of each of the statment in above example:

 # 0: stands for the first molecular orbital (it can plot from 0 to the total number of orbitals that the molecule has)
 # "xy": stands for the plane in which the graph will show, this method accepts the following commands: "x", "y", "z", "xy", "xz" and "yz"
 # -4 the distance where the plot starts in a.u
 # -4 the distance where the plot ends in a.u
 # the number of points for the surface plot
 # the maximun value of density to plot, electron density has cusps and therefore sometimes the graph will show only the cusps of some atoms. This part of the method is to help visualize other elementrs that might not be visible
 # 0 or 1: 0 for a surface plot and 1 for a isodensities plot
 # number of isodensities to be plotted

 # Note: originally the object.graph() method was created a way to visualize the electron density and intended to test the reconstruction of the electron density. In the future a better approach to plor surface and isodensities will be developed.

 
