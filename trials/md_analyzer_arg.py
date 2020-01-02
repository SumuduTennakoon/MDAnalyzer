# -*- coding: utf-8 -*-
# MOLEDULAR DYNAMICS RESULTS ANALYSER
# PROCESS XDATCAR, PLOT ATOMIC POSITIONS AND MOTION.
# COMPUTES VDATCAR AND VELOCITY AUTOCORELATION, GENERATE POWER SPECTRUM.
# PERFOM BOND LEGTH AND ANGLE STATISTICS
# Sumudu Tennakoon (Created 2017-10-22, Last update: 2018-06-13) Argument Parse Version

import argparse
import sys

import pandas as pd
import numpy as np
import metrical

def cords2vector(x, y, z):
    return np.matrix([[x], [y], [z]])

def basis2unitcell(basis):
    #Lattice vectors in xyz cooridinates
    a_xyz= cords2vector(basis[0][0], basis[0][1], basis[0][2])
    b_xyz= cords2vector(basis[1][0], basis[1][1], basis[1][2])
    c_xyz= cords2vector(basis[2][0], basis[2][1], basis[2][2])
      
    #lattice constants
    a= np.sqrt((a_xyz.T*a_xyz).item(0))
    b= np.sqrt((b_xyz.T*b_xyz).item(0))
    c= np.sqrt((c_xyz.T*c_xyz).item(0))   

    #calculate angles between latttice vectors          
    cos_alpha = ((b_xyz.T*c_xyz).item(0))/(b*c)
    cos_beta= ((c_xyz.T*a_xyz).item(0))/(c*a)
    cos_gamma= ((a_xyz.T*b_xyz).item(0))/(a*b)

    alpha= np.degrees(np.arccos(cos_alpha))
    beta= np.degrees(np.arccos(cos_beta))
    gamma= np.degrees(np.arccos(cos_gamma))    
    
    return a, b, c, alpha, beta, gamma      

def displacements(self, Frame0=None, Frame1=None):    
    X0 = Frame0.positions()
    X1 = Frame1.positions()
    return X1-X0

def velocity(self, Frame0=None, Frame1=None, dt=1):        
    return self.atom_uvw
    
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
       
    def positions(self):        
        return self.atom_uvw


    
class XDATCAR():
    def __init__(self, file_name='XDATCAR', file_type=None, version=None):
        self.file_name=file_name
        self.frames=[] #Stores an array of frames
        self.elements=[] # Atom element Symbol
        self.atoms_count=[] #Count of atome elements
        self.atoms=[]
        self.atom_ids = []
        self.total_atoms = 0
        self.a = 1.0
        self.b = 1.0
        self.c = 1.0
        self.UC = None

    # calculate displacement of atoms
    def displacements(self,X1,X2):
        return X2-X1
    
    # Read single frame in XDATCAR file
    def read_frame(self, frame_lines, frame_id=1, print_frame=False):               
        atom_uvw=[]
        
        # extracting atom position data
        for n in range(len(frame_lines)):
            atom_uvw.append(frame_lines[n].split())  # atom positions
        atom_uvw=np.array(atom_uvw,dtype=float)
        
        # create new frame
        frame=FRAME(frame_id, self.atom_ids, atom_uvw)
        
        if(print_frame==True):
            print(str(frame))
            
        return frame
    
    # load data from XDATCAR file        
    def read_vasp(self, file_name, version=5, frame_limit=None):
        self.file_name=file_name
        # read XDATCAR file content to memory
        with open(file_name, 'r') as f:
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
        
        print(self.elements)
        print(self.atoms_count)
        # assign atom ids
        for i in range(0, len(self.elements)):
            for j in range(0, self.atoms_count[i]):
                self.atom_ids.append("%s%d"%(self.elements[i],j+1)) 
                
        # compute unit cell from basis
        a, b, c, alpha, beta, gamma = basis2unitcell(self.basis)
        self.UC = metrical.UnitCell(a, b, c, alpha, beta, gamma) 

        #Display File Header Information
        print("Description:\t"+str(self.header))
        print("Scale Factor:\t"+str(self.scale_factor))
        print("Lattice Vectors:\n"+str(self.basis))
        print("Elements:\t"+str(self.elements))
        print("Atoms Count:   \t"+str(self.atoms_count))
        print("Total Atoms:   \t"+str(self.total_atoms))      
        print("Lattice Parameters:")
        print("a \tb  \tc \talpha\tbeta\tgamma\n%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f"%(a,b,c,alpha,beta,gamma))
        print("unit-cell Volume: %.4f angstrom"%self.UC.volume())
        
        #Set number of frames to read
        if frame_limit==None:
            frame_limit = len(xdatcar[8:]/(self.total_atoms+1))
            
        #Start Reading Frames        
        self.frames_count=0
        self.frames=[]
        
        #Read and disiplay first frame
        frame_start = 8
        frame_end = frame_start + self.total_atoms
        frame_lines = xdatcar[frame_start:frame_end]
        self.read_frame(frame_lines=frame_lines, frame_id=1, print_frame=True) 
        self.frames_count=1
            
        for i in range(1,frame_limit):
            # Read frames
            frame_start = 8 + i*(self.total_atoms+1)
            frame_end = frame_start + self.total_atoms
            frame_lines = xdatcar[frame_start:frame_end] 
            self.frames.append(self.read_frame(frame_lines=frame_lines, frame_id=i+1, print_frame=False))  
            self.frames_count=i
            
        del(xdatcar)    
        print("%d Frames Extracted"%self.frames_count)
        
if __name__== "__main__":
    try:      
        parser = argparse.ArgumentParser(description='pyXDATCAR')
        parser.add_argument(dest='path', metavar='path', type=str, help='[File Path]/FileName')
        parser.add_argument(dest='frames', metavar='-frames', nargs='?', type=str, help='Frames to extract | optional')
        parser.add_argument(dest='ftype', metavar='-ftype', nargs='?', type=str, help='File type (E.g. vasp) | optional')
        parser.add_argument(dest='version', metavar='-version', nargs='?', type=str, help='Version | optional')
        args = parser.parse_args()
        
        print('Example usage: >python md_analysis_v2.py XDATCAR100 10')
        path = args.path
        if args.frames==None:
            frames=1
        else:
            frames = int(args.frames)
        ftype = args.ftype
        version = args.version
        
        print('File path: {}'.format(path))
        print('Frames to extract: {}'.format(frames))
        print('File type: {}'.format(ftype))
        print('File version: {}'.format(version))        
    except:
        print('Error parsing the input parameters {}'.format(sys.exc_info()[0]))
        sys.exit(1)   
    
    Displacements = XDATCAR()
    Displacements.read_vasp(path, version=version, frame_limit=frames)
    sys.exit(0)
    
    
    
        