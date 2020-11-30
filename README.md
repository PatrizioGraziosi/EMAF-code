# EMAF-code

The Effective Mass Finder (EMAF) code aims to compute the meaningful effective mass parameters for a given, generic and arbitrary, three-dimensional (3D) bandstructure of semiconductors. In particular the EMAF code extracts the density of states (DOS) effective mass, and the conductivity effective mass. We introduce the meaning of the effective masses that the EMAF code computes, and the utilization of the code. 

The code can be used by editing a text file or via a graphical user interface (GUI). The former is contained in a folder named “EMAF_v2 release” while the latter can be downloaded as an app anmed “EMAP_app”. The functionality is the same. In the folder “EMAF_v2” there is a script named “EMAF_run”; this is the only file that needs to be modified as explained in section 4 and launched. When using the EMAF app, these instruction will be edited in the GUI.


# The core idea

The EMAF code extracts the effective masses values which are meaningful for the charge transport, for a 3D semiconductor electronic structure. These are the DOS effective mass and the conductivity efective mass. The EMAF objective is to dleiver effective mass values which are representative of the band structure as a whole, not of a single specific valley. 
Thus, the DOS effective mass is the effective mass of an isotropic parabolic band that gives the same carrier density of the actual bandstructure under investigation. 
The conductivity effective mass is the effective mass of an isotropic parabolic band that gives the same injection velocity as a ballistic field effect transistor (FET) having a channel with the investigated bandstructure and operated in the subthreshold regime. Thus, the computed conductivity effective mass grasps the essence of the bandstructure ballistic transport properties. It is also separately computed for the three space directions.
