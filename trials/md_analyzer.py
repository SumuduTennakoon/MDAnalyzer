# MOLEDULAR DYNAMICS RESULTS ANALYSER
# PROCESS XDATCAR, PLOT ATOMIC POSITIONS AND MOTION.
# COMPUTES VDATCAR AND VELOCITY AUTOCORELATION, GENERATE POWER SPECTRUM.
# PERFOM BOND LEGTH AND ANGLE STATISTICS
# Sumudu Tennakoon (Created: 2017-10-22)

from math import *
import numpy as np
from metrical import * 
import time
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
FN =2000 #Number of frams to be processed

###############################################################################
# Wrapper functions(use untions)
def xdatcar2vdatcar(XDATCARFILE='XDATCAR', VDATCARFILE='VDATCAR', OutValues=False, OutFile=True):
    '''
    Function computes velocity profile (VDATCAR) given an XDATCAR file from vasp output.
    '''
    if OutFile == True:
        # Outputs VDATCAR file
    if OutValues == True:
        # return velocity vs. time
        return VDATCAR

def velocitycorr(VDATCARFILE='VDATCAR', VELCORRFILE = 'VCORR', dt=None, N=None, interval=None, OutValues=True, OutFile=False):
    '''
    Function computes power spectrum (PSPectra) given velocity profile (VDATCAR).
    '''
    if OutFile == True:
        # Outputs VCORR file
    if OutValues == True:
        # return power spectrum
        return VDATCAR    
###############################################################################    
def element_constant(element):
    if element=="H":
        mass = 1.00794 #u
        radius = 25 #pm
        color = "g"
    elif element=="O":
        mass = 15.9994
        radius = 60    
        color="r"
    elif element=="Mg":
        mass = 24.3050 
        radius = 150 
        color="b"
    else:
        mass = 1.0
        radius = 1.0   
        color = "k"       
        
    return mass,radius,color
        
def cords2vector(x, y, z):
    return np.matrix([[x], [y], [z]])

class FRAME(): 
    def __init__(self,frame_id=0, atom_id=[], atom_uvw=[]):
        self.frame_id=frame_id
        self.atom_id = atom_id
        self.atom_uvw = atom_uvw
        
    def __str__(self): 
        out = "# u,v,w in fractional cooridnates\n"
        out = out + "#\tAtom_ID\tu        \tv        \tw       \n"
        for i in range(len(self.atom_id)):
            out = out+"%3d\t%s\t%.6f\t%.6f\t%.6f\n"%((i+1), self.atom_id[i], self.atom_uvw[i][0], self.atom_uvw[i][1], self.atom_uvw[i][2])
        return out
    
    def parse_frame(self,frame_lines, elements):
        uvw=[]
        for line in frame_lines:
            uvw.append(line.split())
        uvw=np.array(uvw,dtype=float)
        self.add_atom
    
    def positions(self):        
        return self.atom_uvw
        
class XDATCAR():
    def __init__(self, filename='XDATCAR'):
        self.filename=filename
        self.frames=[]
        self.V=[]
        self.atom_element=[]
        self.atom_mass=[]
        self.atom_radius=[] 
        self.atom_id=[]
        self.atom_color=[] #color used to plot atom

    # calculate displacement of atoms
    def displacements(self,X1,X2):
	return X2-X1
    
    # calculate velocity of atoms		
    def velocities(self,frame1,frame2,UC,dt=1.0):
        a=UC.a
        b=UC.b
        c=UC.c
        R1=frame1.positions()
        R2=frame2.positions()
        R=self.displacements(R1,R2) #abc2xyz(X1)-abc2xyz(X2)
        
        # applying the periodic boundary conditions
        for i in range(len(R)):
            while(R[i][0]<-0.5):
                R[i][0]=R[i][0]+1.0
            while(R[i][1]<-0.5):
                R[i][1]=R[i][1]+1.0
            while(R[i][2]<-0.5):
                R[i][2]=R[i][2]+1.0
            while(R[i][0]>0.5):
                R[i][0]=R[i][0]-1.0
            while(R[i][1]>0.5):
                R[i][1]=R[i][1]-1.0
            while(R[i][2]>0.5):
                R[i][2]=R[i][2]-1.0 
                
        # convert fracional coordinates to cartecian cooridnates [use of matrix algebra (A*B).T = (B.T)(A.T)]
        R=np.array(np.matrix(R)*(self.UC.M.T))
        V=R/dt #Velocity unit: A/fs 
        
        return V

    # Generate velocity data (VDATCAR file)
    def vdatcar(self,file_name='VDATCAR'):
        # generate vdatcar
        vdatcarfile=""
        vdatcarfile=vdatcarfile+"%s\n"%(self.header)
        vdatcarfile=vdatcarfile+"%10s\t%10s\t%10s\n"%(self.elements[0],self.elements[1],self.elements[2])
        vdatcarfile=vdatcarfile+"%10d\t%10d\t%10d\n"%(self.atoms_count[0],self.atoms_count[1],self.atoms_count[2])
        
        for i in range(1,len(self.frames)):
            V=self.velocities(self.frames[i-1],self.frames[i],self.UC)
            self.V.append(V)
            vdatcarfile=vdatcarfile+"Frame: %d\n"%(i)
            for j in range(len(V)):
                vdatcarfile=vdatcarfile+"  % .8f\t% .8f\t% .8f\n"%(V[j][0],V[j][1],V[j][2])
        self.V=np.array(self.V)
        vac(self.V,self.atom_mass)
        fv = open(file_name,'w')
        fv.write(vdatcarfile)
        fv.close()
        del(vdatcarfile)
    
    # Calculate O-H bond length and andles
    def OH_bond_stats(self):
        OH_bond_length=[]
        OH_bond_angle=[]
        OH_bond_angle2D=[]
        for frame in self.frames:
            lstO=[]
            lstH=[]
            a_ax=np.matrix([[1],[0],[0]])
            b_ax=np.matrix([[0],[1],[0]])
            c_ax=np.matrix([[0],[0],[1]])
            for i in range(len(self.atom_element)):
                if self.atom_element[i]=="H":
                    lstH.append(i)
                elif self.atom_element[i]=="O":
                    lstO.append(i)
            OH_bond_info=""
            for i in lstO:
                rO=np.matrix([[frame.atom_uvw[i][0]],[frame.atom_uvw[i][1]],[frame.atom_uvw[i][2]]])
                for j in lstH:
                    rH=np.matrix([[frame.atom_uvw[j][0]],[frame.atom_uvw[j][1]],[frame.atom_uvw[j][2]]])
                    OH_bond=bond_length(rO, rH,self.UC)
                    if OH_bond < 1.05:
                        OH_bond_length.append(OH_bond)
                        rOH2D=np.matrix([[(rH-rO).item(0)],[(rH-rO).item(1)],[0.0]])
                        OH_bond_angle2D.append(angle_delta(rOH2D,a_ax,self.UC))                        
                        OH_bond_angle.append(angle_delta((rH-rO),c_ax,self.UC))   
                        OH_bond_info=OH_bond_info+"%.6f\t%.6f\t%.6f"%(0,0,0) 
                        
        fb = open("OHBONDS",'w')
        fb.write(OH_bond_info)
        fb.close()
        
        plt.figure(100)
        n, bins, patches = plt.hist(OH_bond_length, bins=1000, histtype='step')
        
        plt.figure(200)
        n, bins, patches = plt.hist(OH_bond_angle2D, bins=1000, histtype='step')
        
        plt.figure(300)
        n, bins, patches = plt.hist(OH_bond_angle, bins=1000, histtype='step')
        
        plt.show()       
            
    # load data from XDATCAR file        
    def load(self):
        # read XDATCAR file content to memory
        with open(self.filename, 'r') as f:
            xdatcar=f.readlines()
       
        # extract file header
        self.header = xdatcar[0].strip()
        # scale factor
        self.scale_factor= float(xdatcar[1].strip())

        # lattice basis vectors
        self.basis = []
        for i in [2,3,4]:
            self.basis.append(xdatcar[i].split())
        self.basis=np.array(self.basis,dtype=float)
        
        # elements and atoms count
        self.elements = xdatcar[5].split()
        self.atoms_count = np.array(xdatcar[6].split(),dtype=int)
        self.total_atoms = sum(self.atoms_count)
        
        # extracting atom data
        for i in range(len(self.elements)):
            for j in range(self.atoms_count[i]):
                elem = self.elements[i]
                mass,radius,color = element_constant(elem)
                self.atom_element.append(elem)  
                self.atom_mass.append(mass)
                self.atom_radius.append(radius)
                self.atom_id.append("%s%d"%(elem,j+1))     
                self.atom_color.append(color)
                    
        #Lattice vectors in xyz cooridinates
        a_xyz=cords2vector(self.basis[0][0], self.basis[0][1], self.basis[0][2])
        b_xyz=cords2vector(self.basis[1][0], self.basis[1][1], self.basis[1][2])
        c_xyz=cords2vector(self.basis[2][0], self.basis[2][1], self.basis[2][2])
        
        self.a_xyz=a_xyz
        self.b_xyz=b_xyz
        self.c_xyz=c_xyz
        
        #lattice constants
        a=sqrt((a_xyz.T*a_xyz).item(0))
        b=sqrt((b_xyz.T*b_xyz).item(0))
        c=sqrt((c_xyz.T*c_xyz).item(0))   

        #calculate angles between latttice vectors          
        cos_alpha = ((b_xyz.T*c_xyz).item(0))/(b*c)
        cos_beta= ((c_xyz.T*a_xyz).item(0))/(c*a)
        cos_gamma= ((a_xyz.T*b_xyz).item(0))/(a*b)

        alpha=degrees(acos(cos_alpha))
        beta=degrees(acos(cos_beta))
        gamma=degrees(acos(cos_gamma))       

        self.UC=UnitCell(a, b, c, alpha, beta, gamma)

        print("Description:\t"+str(self.header))
        print("Scale Factor:\t"+str(self.scale_factor))
        print("Lattice Vectors:\n"+str(self.basis))
        print("Elements:\t"+str(self.elements))
        print("Atoms:   \t"+str(self.atoms_count))
        print("Total Atoms:   \t"+str(self.total_atoms))      
        print("Lattice Parameters:")
        print("a \tb  \tc \talpha\tbeta\tgamma\n%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f"%(a,b,c,alpha,beta,gamma))
        print("unit-cell Volume: %.4f angstrom"%self.UC.volume())
        
        self.frames_count=0
        self.frames=[]

        self.read_frame(xdatcar[8:8+self.total_atoms],1,True) #dsiplay first frame
        
        for i in range(FN):
            # extracting frame data
            frame_start = 8+i*(self.total_atoms+1)
            frame_end = frame_start+self.total_atoms
            frame_lines = xdatcar[frame_start:frame_end] 
            self.frames.append(self.read_frame(frame_lines,i+1, False))  
            self.frames_count=i+1
            
        print "%d Frames Extracted"%self.frames_count
        del(xdatcar)     
    
    # plot unitcell boundaries
    def plot_unitcell(self,fig,ax):
        from matplotlib.path import Path
        import matplotlib.patches as patches
        verts = [(0., 0.),(self.a_xyz.item(0), self.a_xyz.item(1)),(self.a_xyz.item(0)+self.b_xyz.item(0), self.a_xyz.item(1)+self.b_xyz.item(1)),(self.b_xyz.item(0), self.b_xyz.item(1)),(0., 0.)] 
        path=Path(verts)
        patch = patches.PathPatch(path, fill=None, lw=2)
        ax.add_patch(patch)#    
    
    # Plot single frame in XDATCAR file        
    def plot_frame(self,frame_id,fig,x_max=1.0, x_min=0.0, y_max = 1.0, y_min=0.0, z_max=1.0, z_min = 0.0):  
        ax = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
        
        ax.set_title("a-b")
        ax2.set_title("b-c")
        ax3.set_title("a-c")
        ax4.set_title("a-b layer")
        self.plot_unitcell(fig,ax)
        self.plot_unitcell(fig,ax4)
        
        fig.suptitle("After %d frames"%frame_id)  
        
        frame = self.frames[frame_id-1]
        for i in range(len(self.atom_id)):
            r=np.matrix([frame.atom_uvw[i][0], frame.atom_uvw[i][1],frame.atom_uvw[i][2]]).T
            scale = self.atom_radius[i]
            color = self.atom_color[i]
            x=self.UC.M*r                                
            ax.scatter(x.item(0),x.item(1),c=color,s=scale*2.0)
            ax2.scatter(x.item(1),x.item(2),c=color,s=scale*2.0)
            ax3.scatter(x.item(0),x.item(2),c=color,s=scale*2.0)
            if frame.atom_uvw[i][2]<z_max and frame.atom_uvw[i][2]>z_min:
                ax4.scatter(x.item(0),x.item(1),c=color,s=scale*2.0)
                #ax4.annotate(frame.atom_id[i],  xy=(x.item(0),x.item(1)))
        del(frame)
    
    # Read single frame in XDATCAR file
    def read_frame(self, frame_lines,frame_id=1, print_frame=False):               
        atom_uvw=[]
        
        # extracting atom position data
        for n in range(len(frame_lines)):
            atom_uvw.append(frame_lines[n].split())  # atom positions

        atom_uvw=np.array(atom_uvw,dtype=float)
        # create new frame
        frame=FRAME(frame_id, self.atom_id, atom_uvw)
        
        if(print_frame==True):
            print(str(frame))
            
        return frame




