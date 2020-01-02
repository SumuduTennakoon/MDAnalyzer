#CALCULATING GEOMETRIC RELATIONS FOR MINERALS USING THE METRICAL MATRIX"
#Sumudu Tennakoon (Created: 2017-10-15)
#Reference: G. V Gibbs, The Metrical Matrix in Teaching Mineralogy (1990).

from metrical import *

print("CALCULATING GEOMETRIC RELATIONS FOR MINERALS USING THE METRICAL MATRIX")
print("Reference: G. V Gibbs, The Metrical Matrix in Teaching Mineralogy (1990).")

###########################################################
#EXAMPLE 1: Constructing a metrical matrix
print("#EXAMPLE 1: Constructing a metrical matrix")
triRhodUC = UnitCell(10.497, 9.797, 12.185, 103.0, 108.51, 82.5)
print("G=\n%s" % str(triRhodUC.G()))
print("-------------------------------------------------------\n")

#EXAMPLE 2: A calculation of the angle between two zones
print("EXAMPLE 2: A calculation of the angle between two zones")
r1 = vector(-1, 1, 2) #(h k l)
r2 = vector(2, 1,3) #(h k l)
print("r1=\n%s"%str(r1))
print("r2=\n%s"%str(r2))
print("Angle = %.2f degrees"%angle_delta(r1,r2,triRhodUC)) 
print("-------------------------------------------------------\n")

#EXAMPLE 3: A calculation of a d(hkl)-spacing
print("#EXAMPLE 3: A calculation of a d(hkl)-spacing")
kyaUC = UnitCell(7.126, 7.852, 5.572, 89.99, 101.11, 106.03)
s = vector(2, -3, 1) #(h k l)
print("s=\n%s"%str(s))
print("G=\n%s" % str(kyaUC.G()))
print("G^-1=\n%s" % str(kyaUC.G().I))
print("d_(2-31)= %.5f angstroms"%d_hkl(s,kyaUC))
print("-------------------------------------------------------\n")

#EXAMPLE 4: Calculation of the interfacial angle between two planes
print("#EXAMPLE 4: Calculation of the interfacial angle between two planes")
s1 = vector(-1, 2, 2) #(h k l)
s2 = vector(3, -2, 0) #(h k l)
print("s1=\n%s"%str(s1))
print("s2=]n%s"%str(s2))
print("Angle = %.2f degrees"%angle_DELTA(s1,s2,kyaUC))
print("-------------------------------------------------------\n")

#EXAMPLE 5: Calculation of the angle between a zone and a face pole
print("#EXAMPLE 5: Calculation of the angle between a zone and a face pole")
s = vector(-1, 2, 2) #(h k l)
r = vector(8, 7 ,2) #[u, v, w]
print("s=\n%s"%str(s))
print("r=\n%s"%str(r))
print("Angle = %.2f degrees"%angle_eps(s,r,kyaUC))
print("-------------------------------------------------------\n")

#EXAMPLE 6: Calculation of the volume of the unit cell
#monoclinic silica
print("#EXAMPLE 6: Calculation of the volume of the unit cell")
mcSiOUC = UnitCell(7.135, 12.372, 7.173, 90.0, 90.0, 120.36)
print("G=\n%s" % str(mcSiOUC.G()))
print("Volume= %.4f angstrom^3"%mcSiOUC.volume())    
print("-------------------------------------------------------\n")

#EXAMPLE 7: Calculation of SiO bond lengths
print("# EXAMPLE 7: Calculation of SiO bond lengths")
SiOUC = UnitCell(7.988, 7.040, 7.025, 90.51, 95.18, 102.47)
v1 = vector(0.3955, 0.9092, 0.2746) #0-O8
v2 = vector(0.2150, 0.9544, 0.3440) #0-Si2
v3 = vector(0.4505, 0.7353, 0.1447) #0-Si3
v4 = v2-v1 #O8-Si2
v5 = v3-v1 #O8-Si3
v4_length = bond_length(v2,v1, SiOUC) 
v5_length = bond_length(v3,v1, SiOUC)
print("Unit Cell:\na =%.5f angstroms\nb =%.5f angstroms\nc =%.5f angstroms\nalpha =%.5f degrees\nbeta =%.5f degrees\ngamma =%.5f degrees"%(7.988, 7.040, 7.025, 90.51, 95.18, 102.47))
print("Volume= %.4f angstrom^3"%SiOUC.volume())
print("G=\n%s" % str(SiOUC.G()))
print("v_O8=\n"+str(v1))
print("v_Si2=\n"+str(v2))
print("v_Si3=\n"+str(v3))
print("v_O8-Si2=\n"+str(v4))
print("v_O8-Si3=\n"+str(v5))
print("Length(O8-Si2)= %.6f  angstroms" % v4_length)
print("Length(O8-Si3)= %.6f  angstroms" % v5_length)
print("-------------------------------------------------------\n")

#EXAMPLE 8: A calculation of the SiOSi angle  
print("# EXAMPLE 8: A calculation of the SiOSi angle")
print("Angle(Si3O8Si2) = %.2f degrees" % angle_theta(v4, v5, SiOUC))
