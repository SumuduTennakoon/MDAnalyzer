# MOLEDULAR DYNAMICS RESULTS ANALYSER
# PROCESS XDATCAR, PLOT ATOMIC POSITIONS AND MOTION.
# COMPUTES VDATCAR AND VELOCITY AUTOCORELATION, GENERATE POWER SPECTRUM.
# PERFOM BOND LEGTH AND ANGLE STATISTICS
# Sumudu Tennakoon (Created 2017-10-22, Last update: 2019-12-02)

from math import *
import numpy as np
from metrical import * 
import time
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
FN =2000 #Number of frams to be processed

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

def vac(vframes, filename="XDATCAR", atoms=None):
    atoms=None
    dt=1e-15
    frame_count=len(vframes)
    vacf=np.zeros(frame_count)
    #print("Frames "+str(frame_count))
    #print(atoms)
    for t in range(0,frame_count):
        num = 0.0
        denom = 0.0    
        for t0 in range(0,1): #frame_count): #Corelate with origin
            for i in range(0,len(vframes[t])):
                #print(i)
                if atoms != None:
                    if i not in atoms:
                        continue
                        
                vx0i = vframes[t0][i][0]
                vy0i = vframes[t0][i][1]
                vz0i = vframes[t0][i][2]            
                vxti = vframes[t][i][0]
                vyti = vframes[t][i][1]
                vzti = vframes[t][i][2]
                
                v0vt = vx0i*vxti + vy0i*vyti + vz0i*vzti
                v0v0 = vx0i*vx0i + vy0i*vy0i + vz0i*vz0i

                num = num + v0vt
                #denom = denom + v0v0
                denom = 1.0 # No normalization
                
        vacf[t] = num/denom
    
    # time steps
    n = frame_count*20
    #n = frame_count
    
    ts=np.linspace(0, dt*frame_count, frame_count)
    #freq = np.fft.fftfreq(frame_count, dt)/1e12*4.13567 #meV
    freq = np.fft.fftfreq(n, dt)/1e12*33.35641 #conversion to cm-1
    
    # WINDOW
    from scipy.signal import hanning
    window = np.hanning(len(vacf)) #Apply Hann window  
    #window = 1.0
    
    ps = np.fft.fft(vacf*window, n=n)
#    gaussian = np.exp(-(freq/20.0)**2/2)
#    ps = np.convolve(ps, gaussian, mode="full")   
    ps = ps.real**2+ps.imag**2
 

    vacf_File = ""
    ps_file=""
	
    #Writes Power Spectrum PSPEC file
    for i in range(len(ps)//2):
        ps_file=ps_file+"%.2f\t%.4f\n"%(freq[i],ps[i])	
    ff = open(filename+"_"+"PSPEC",'w')
    ff.write(ps_file)
    ff.close()    
    
    #Writes Velocity Autocorelation Function to VACF file
    for i in range(len(vframes)):
        vacf_File=vacf_File+"%e\t%.4f\n"%(ts[i],vacf[i])	
    ff = open(filename+"_"+"VACF",'w')
    ff.write(vacf_File)
    ff.close()   
    
    plt.figure(1)
    plt.plot(ts*1e12, vacf,label=filename)
    plt.xlabel("Time [ps]")
    plt.ylabel("VACF")
    plt.legend()
    plt.figure(2)
    #plt.psd(vacf, 512, 1.0/dt,pad_to=4096)
    plt.plot(freq,ps,"-",label=filename)  
    plt.xlabel("Frequency [cm-1]")
    plt.ylabel("Vibrational Density")
    plt.legend()
    plt.xlim(0,4500)
    
def vac0(vframes,filename="XDATCAR"): #was used on 2019-12-02
    nom=0.0
    denom=0.0
    dt=1e-15
    frame_count=len(vframes)
    vacf=np.zeros(frame_count)
    #print("Frames "+str(frame_count))

    for t in range(0,frame_count):
        vact = 0.0
        for i in range(0,len(vframes[t])):
            vx0i = vframes[0][i][0]
            vy0i = vframes[0][i][1]
            vz0i = vframes[0][i][2]            
            vxti = vframes[t][i][0]
            vyti = vframes[t][i][1]
            vzti = vframes[t][i][2]
            #Use to normalize
            vti = sqrt(vxti**2+vyti**2+vzti**2)
            v0i = sqrt(vx0i**2+vy0i**2+vz0i**2)
            vact = vact + (vx0i*vxti+vy0i*vyti+vz0i*vzti)/(vti*v0i)
        vacf[t] = vact/len(vframes[t])
    
    # time steps
    n = frame_count*10
    #n = frame_count
    
    ts=np.linspace(0, dt*frame_count, frame_count)
    #freq = np.fft.fftfreq(frame_count, dt)/1e12*4.13567 #meV
    freq = np.fft.fftfreq(n, dt)/1e12*33.35641
    
    # WINDOW
    from scipy.signal import hanning
    window = np.hanning(len(vacf)) #Apply Hann window  
    #window = 1.0
    
    ps = np.fft.fft(vacf*window, n=n)
#    gaussian = np.exp(-(freq/20.0)**2/2)
#    ps = np.convolve(ps, gaussian, mode="full")   
    ps = ps.real**2+ps.imag**2
 

    vacf_File = ""
    ps_file=""
	
    #Writes Power Spectrum PSPEC file
    for i in range(len(ps)//2):
        ps_file=ps_file+"%.2f\t%.4f\n"%(freq[i],ps[i])	
    ff = open(filename+"_"+"PSPEC",'w')
    ff.write(ps_file)
    ff.close()    
    
    #Writes Velocity Autocorelation Function to VACF file
    for i in range(len(vframes)):
        vacf_File=vacf_File+"%e\t%.4f\n"%(ts[i],vacf[i])	
    ff = open(filename+"_"+"VACF",'w')
    ff.write(vacf_File)
    ff.close()   
    
    plt.figure(1)
    plt.plot(ts*1e12, vacf,label=filename)
    plt.xlabel("Time [ps]")
    plt.ylabel("VACF")
    plt.legend()
    plt.figure(2)
    #plt.psd(vacf, 512, 1.0/dt,pad_to=4096)
    plt.plot(freq,ps,"-",label=filename)  
    plt.xlabel("Frequency [cm-1]")
    plt.ylabel("Vibrational Density")
    plt.legend()
    plt.xlim(0,4500)
      
def vac1(vframes,filename="XDATCAR"):
    nom=0.0
    denom=0.0
    dt=1e-15
    frame_count=len(vframes)
    vacf=np.zeros(frame_count)
    #print("Frames "+str(frame_count))
    tmax_frame=frame_count
    for t in range(0,frame_count): 
        vact = 0.0
        for t0 in range(0,tmax_frame-t):
            vact0 = 0.0
            for i in range(0,len(vframes[t0])):
                vx0i = vframes[0][i][0]
                vy0i = vframes[0][i][1]
                vz0i = vframes[0][i][2]            
                vxti = vframes[t0][i][0]
                vyti = vframes[t0][i][1]
                vzti = vframes[t0][i][2]
                #Use to normalize
                vti = sqrt(vxti**2+vyti**2+vzti**2)
                v0i = sqrt(vx0i**2+vy0i**2+vz0i**2)
                vact0 = vact0 + (vx0i*vxti+vy0i*vyti+vz0i*vzti)/(vti*v0i)
            vact = vact + vact0/(tmax_frame-t)
        vacf[t] = vact/len(vframes[t])

    ts=np.linspace(0, dt*frame_count, frame_count)

    from scipy.signal import hanning
    window = np.hanning(len(vacf)) #Apply Hann window  
    ps = np.fft.fft(vacf*window, n=20480)
    ps = ps.real**2+ps.imag**2
    #freq = np.fft.fftfreq(frame_count, dt)/1e12*4.13567 #meV
    freq = np.fft.fftfreq(20480, dt)/1e12*33.35641
#    gaussian = np.exp(-(freq/20.0)**2/2)
#    ps = np.convolve(ps, gaussian, mode="full")    
    vacf_File = ""
    ps_file=""
	
    #Writes Power Spectrum PSPEC file
    for i in range(len(ps)/2):
        ps_file=ps_file+"%.2f\t%.4f\n"%(freq[i],ps[i])	
    ff = open(filename+"_"+"PSPEC",'w')
    ff.write(ps_file)
    ff.close()    
    
    #Writes Velocity Autocorelation Function to VACF file
    for i in range(len(vframes)):
        vacf_File=vacf_File+"%e\t%.4f\n"%(ts[i],vacf[i])	
    ff = open(filename+"_"+"VACF",'w')
    ff.write(vacf_File)
    ff.close()   
    
    plt.figure(1)
    plt.plot(ts*1e12, vacf,label=filename)
    plt.xlabel("Time [ps]")
    plt.ylabel("VACF")
    plt.legend()
    plt.figure(2)
    #plt.psd(vacf, 512, 1.0/dt,pad_to=4096)
    plt.plot(freq,ps,"-",label=filename)  
    plt.xlabel("Frequency [cm-1]")
    plt.ylabel("Vibrational Density")
    plt.legend()
    plt.xlim(0,4500)


def cords2vector(x, y, z):
    return np.matrix([[x], [y], [z]])

class FRAME(): 
    def __init__(self,frame_id=1, atom_id=[], atom_uvw=[]):
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
            V=self.velocities(self.frames[i-1], self.frames[i], self.UC)
            self.V.append(V)
            vdatcarfile=vdatcarfile+"Frame: %d\n"%(i)
            for j in range(len(V)):
                vdatcarfile=vdatcarfile+"  % .8f\t% .8f\t% .8f\n"%(V[j][0],V[j][1],V[j][2])
        self.V=np.array(self.V)

        # Isolate Atoms
        lstH = []
        lstO = []
        for i in range(len(self.atom_element)):
            if self.atom_element[i]=="H":
                lstH.append(i)
            elif self.atom_element[i]=="O":
                lstO.append(i)
        atoms = lstO + lstH
        
        vac(self.V, self.filename, atoms=atoms)
        
        fv = open(self.filename+"_"+file_name,'w')
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
        n, bins, patches = plt.hist(OH_bond_length, bins=1000, histtype='step', label=self.filename)
        plt.xlabel('bond length')
        #plt.ylabel('count')
        plt.legend()
        
        plt.figure(200)
        n, bins, patches = plt.hist(OH_bond_angle2D, bins=1000, histtype='step', label=self.filename)
        plt.xlabel('bond angle')
        plt.legend()
        
        plt.figure(300)
        n, bins, patches = plt.hist(OH_bond_angle, bins=1000, histtype='step', label=self.filename)
        plt.xlabel('bond angle')
        plt.legend()
        
        #plt.show()       
            
    # load data from XDATCAR file        
    def load(self,use_frames=1):
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
        
        for i in range(use_frames):
            # extracting frame data
            frame_start = 8+i*(self.total_atoms+1)
            frame_end = frame_start+self.total_atoms
            frame_lines = xdatcar[frame_start:frame_end] 
            self.frames.append(self.read_frame(frame_lines,i+1, False))  
            self.frames_count=i+1
            
        print("%d Frames Extracted"%self.frames_count)
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

# Analyse MD output given XDATCAR
def analyze_md_output(xdatcar_file="XDATCAR",fig_num=000, plot_structure=False, calc_velocities=False, calc_bond_stats=False, use_frames=1):
    start = time.time()
    X=XDATCAR(xdatcar_file)
    print("\nLoading XDARCAR...")
    X.load(use_frames)
    
    if plot_structure==True:
        print("\nPlotting atom positions...")
        fig = plt.figure(fig_num)
        for i in range(0,use_frames,10):
            X.plot_frame(i,fig,z_max=.250)
            
    if calc_velocities==True:
        print("\nCalculating velocities...")
        X.vdatcar()
    
    if calc_bond_stats==True:
        print("\nCalculating O-H bond length and angles...")
        X.OH_bond_stats()

    end=time.time()    
    print("Run time= %.3f seconds"%(end-start))



