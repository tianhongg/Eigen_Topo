import numpy as np
import math 
from scipy import linalg

class Topo:
    def __init__(self,N,dr,r00,Bz,By,r0,l, DenProf=0, ByEProf=0):
        #paratmers
        self.N  = N
        self.dr = dr
        self.r00 = r00  
        self.Bz = Bz
        self.By = By
        self.r0 = r0
        self.l = l
        self.DenProf = DenProf
        self.ByEProf = ByEProf
        # relfection matrix
        self.Refl = np.eye(9)
        self.Refl[0,0] = -1
        self.Refl[1,1] = -1 
        self.Refl[3,3] = -1 
        self.Refl[4,4] = -1 
        self.Refl[6,6] = -1 
        self.Refl[7,7] = -1 
        
    def r_int(self, i):
        return self.r00+(i)*self.dr
    
    def r_half(self, i):
        return self.r00+(i+0.5)*self.dr
    
    def Hminus(self,m, kz, i):
        ri = self.r_int(i)
        Hm = np.zeros((9, 9), dtype=complex)
        Hm[3,7] =  kz/2;
        Hm[3,8] = -m/2/ri;
        Hm[4,6] = -kz/2;
        Hm[4,8] =  1j/self.dr;
        Hm[5,6] =  m/2/ri;
        Hm[5,7] = -1j/self.dr+1j/2/ri;
        return Hm
        
    def Hplus(self,m,kz,i):
        rh = self.r_half(i)
        Hp = np.zeros((9, 9), dtype=complex)
        Hp[6,4] = -kz/2;
        Hp[6,5] =  m/2/rh;
        Hp[7,3] =  kz/2;
        Hp[7,5] =  1j/self.dr;
        Hp[8,3] = -m/2/rh;
        Hp[8,4] = -1j/self.dr-1j/2/rh;
        return Hp
        
    def HH(self,m,kz,i,omegap,by):
        ri = self.r_int(i)
        rh = self.r_half(i);
        H0 = np.zeros((9, 9), dtype=complex)

        H0[0,1] = -1j*self.Bz
        H0[0,2] =  1j*by
        H0[0,3] =  1j*omegap
        
        H0[1,0] =  1j*self.Bz
        H0[1,4] =  1j*omegap
        
        H0[2,0] = -1j*by
        H0[2,5] =  1j*omegap
        
        H0[3,0] = -1j*omegap
        H0[3,7] =  kz/2
        H0[3,8] = -m/2/ri
        
        H0[4,1] = -1j*omegap
        H0[4,6] = -kz/2
        H0[4,8] = -1j/self.dr
        
        H0[5,2] = -1j*omegap
        H0[5,6] =  m/2/ri
        H0[5,7] =  1j/self.dr+1j/2/ri
        
        H0[6,4] = -kz/2
        H0[6,5] =  m/2/rh
        
        H0[7,3] =  kz/2
        H0[7,5] = -1j/self.dr
        
        H0[8,3] = -m/2/rh
        H0[8,4] =  1j/self.dr-1j/2/rh
        return H0
        
    def OmegaP(self,i):
        ri = self.r_int(i)
        if self.DenProf==0:
            return 0.5*(math.tanh( (self.r0-abs(ri))/self.l) +1)
        return 1.
    
    def ByExt(self,i):
        ri = self.r_int(i)
        if self.ByEProf==0:
            return 0.5*(math.tanh( (abs(ri)-self.r0)/self.l) +1)*self.By
        return 1.
        
    def BigMatrix(self,m, kz):
        BM = np.zeros((self.N*9, self.N*9), dtype=complex)
        for i in range(0, self.N):
            im = (i-1)*9
            ii = (i  )*9
            ip = (i+1)*9
            Hm = self.Hminus(m,kz,i)
            Hp = self.Hplus( m,kz,i)
            HH = self.HH(    m,kz,i,self.OmegaP(i), self.ByExt(i))
            if i==0:
                HH = HH + np.matmul(HH, self.Refl)
            if i>0:
                BM[ii:ii+9,im:im+9] = Hm.copy()
            if i<self.N-1:
                BM[ii:ii+9,ip:ip+9] = Hp.copy()
            BM[ii:ii+9,ii:ii+9] = HH.copy()
        return BM
                
    def EigTopo(self,m,kz, low=0.01, high=1.2):
        eig = (linalg.eigvals(self.BigMatrix(m, kz), None, False,False)).real
        mask = (low<=eig)&(eig<=high)
        return eig[mask]
            


# flat geometry 
class TopoFlat:
    def __init__(self,N,dx,L0,Bz,By,x0,l, DenProf=0, ByEProf=0, BC=0):
        #paratmers
        self.N  = N
        self.dx = dx
        self.L0 = L0   # min of x domain 
        self.Bz = Bz
        self.By = By
        self.x0 = x0
        self.l = l
        self.DenProf = DenProf
        self.ByEProf = ByEProf
        self.BC = BC;  #0=periodic 1=open

    def x_int(self, i):
        return self.L0+i*self.dx

    def Hminus(self, ky, kz):
        Hm = np.zeros((9, 9), dtype=complex)
        Hm[3,7] =  kz/2
        Hm[3,8] = -ky/2
        Hm[4,6] = -kz/2
        Hm[4,8] =  1j/self.dx
        Hm[5,6] =  ky/2;
        Hm[5,7] = -1j/self.dx
        return Hm
        
    def Hplus(self, ky, kz):
        Hp = np.zeros((9, 9), dtype=complex)
        Hp[6,4] = -kz/2;
        Hp[6,5] =  ky/2;
        Hp[7,3] =  kz/2;
        Hp[7,5] =  1j/self.dx
        Hp[8,3] = -ky/2;
        Hp[8,4] = -1j/self.dx
        return Hp
        
    def HH(self, ky, kz, omegap, by):
        H0 = np.zeros((9, 9), dtype=complex)

        H0[0,1] = -1j*self.Bz
        H0[0,2] =  1j*by
        H0[0,3] =  1j*omegap
        
        H0[1,0] =  1j*self.Bz
        H0[1,4] =  1j*omegap
        
        H0[2,0] = -1j*by
        H0[2,5] =  1j*omegap
        
        H0[3,0] = -1j*omegap
        H0[3,7] =  kz/2
        H0[3,8] = -ky/2
        
        H0[4,1] = -1j*omegap
        H0[4,6] = -kz/2
        H0[4,8] = -1j/self.dx
        
        H0[5,2] = -1j*omegap
        H0[5,6] =  ky/2
        H0[5,7] =  1j/self.dx
        
        H0[6,4] = -kz/2
        H0[6,5] =  ky/2
        
        H0[7,3] =  kz/2
        H0[7,5] = -1j/self.dx
        
        H0[8,3] = -ky/2
        H0[8,4] =  1j/self.dx
        return H0
        
    def OmegaP(self, i):
        x = self.x_int(i)
        x = abs(x)
        #ramp
        if self.DenProf==0:
            if x<self.x0:
                return 1
            elif x>self.x0+self.l:
                return 0
            else:
                return math.sqrt(1-(x-self.x0)/self.l)

        # sharp
        if self.DenProf==1:
            if x<self.x0:
                return 1
            else:
                return 0
            
        return 1.
    
    def ByExt(self, i):
        x = self.x_int(i)
        x = abs(x)
        #tanh function
        if self.ByEProf==0:
            return 0.5*(math.tanh((x-self.x0)/self.l) +1)*self.By
        return 1.
        
    def BigMatrix(self, ky, kz):

        BM = np.zeros((self.N*9, self.N*9), dtype=complex)

        Hm = self.Hminus(ky,kz)
        Hp = self.Hplus( ky,kz)

        for i in range(0, self.N):
            im = (i-1)*9
            ii = (i  )*9
            ip = (i+1)*9
            HH = self.HH(ky,kz,self.OmegaP(i), self.ByExt(i))
            
            # make the BC
            if i==0 and self.BC==0:
                im = (self.N-1)*9
            if i==self.N-1 and self.BC==0:
                ip = 0

            if i>0 or self.BC==0:
                BM[ii:ii+9,im:im+9] = Hm.copy()
            if i<self.N-1 or self.BC==0:
                BM[ii:ii+9,ip:ip+9] = Hp.copy()

            BM[ii:ii+9,ii:ii+9] = HH.copy()


        return BM
                
    def EigTopo(self, ky, kz, low=0.01, high=1.2):
        eig = (linalg.eigvals(self.BigMatrix(ky, kz), None, False,False)).real
        mask = (low<=eig)&(eig<=high)
        return eig[mask]
            
        
    
    
        
    
    
    