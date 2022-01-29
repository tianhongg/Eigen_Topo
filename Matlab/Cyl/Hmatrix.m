%
% the eigen system is
% (H-)(psi-) +(H-\omega)(psi) +(H+)(psi+)=0;

classdef Hmatrix
    methods (Static)
        %-----------------------------------
        function r = r_int(i)
            global dr;
            global r00;
            r=r00+(i-1)*dr;
        end
        
        function r = r_half(i)
            global dr;
            global r00;
            r=r00+(i-1/2)*dr;
        end
        %-----------------------------------
        function Hm = Hminus(m,kz,i)
            global dr;
            Hm = zeros(9,9);
            ri  = Hmatrix.r_int(i);
            %-------------------
            Hm(4,8) =  kz/2;
            Hm(4,9) = -m/2/ri;
            Hm(5,7) = -kz/2;
            Hm(5,9) =  1i/dr;
            Hm(6,7) =  m/2/ri;
            Hm(6,8) = -1i/dr+1i/2/ri;
        end
        
        function Hp = Hplus(m,kz,i)
            global dr;
            Hp = zeros(9,9);
            rh  = Hmatrix.r_half(i);
            %-------------------
            Hp(7,5) = -kz/2;
            Hp(7,6) =  m/2/rh;
            Hp(8,4) =  kz/2;
            Hp(8,6) =  1i/dr;
            Hp(9,4) = -m/2/rh;
            Hp(9,5) = -1i/dr-1i/2/rh;
        end
        
        function H0 = HH(m,kz,op,oy,i)
            %op = omega_p
            global dr;
            global Omega_z;
            ri  = Hmatrix.r_int( i);
            rh  = Hmatrix.r_half(i);
            H0 = zeros(9,9);
            %-----------------------------------------------------
            H0(1,2) = -1i*Omega_z; H0(1,3) = 1i*oy; H0(1,4) = 1i*op;
            H0(2,1) =  1i*Omega_z;                  H0(2,5) = 1i*op;
            H0(3,1) = -1i*oy;                       H0(3,6) = 1i*op;
            
            H0(4,1) =  -1i*op;  H0(4,8) =  kz/2;    H0(4,9) = -m/2/ri;
            H0(5,2) =  -1i*op;  H0(5,7) = -kz/2;    H0(5,9) = -1i/dr;
            H0(6,3) =  -1i*op;  H0(6,7) =  m/2/ri;  H0(6,8) =  1i/dr+1i/2/ri;
            
            H0(7,5) = -kz/2;    H0(7,6) =  m/2/rh;
            H0(8,4) =  kz/2;    H0(8,6) = -1i/dr;
            H0(9,4) = -m/2/rh;  H0(9,5) =  1i/dr-1i/2/rh;
        end
        
        
        %plasma density
        function op = omegap(i)
            % N is the number of points
            global N;
            global r0; global l;
            op = 0;
            if(i<=N&&i>0)
                r = Hmatrix.r_int(i);
                op = 1/2*(tanh((r0-abs(r))/l) + 1);
                op = sqrt(op);
            end
        end
        
        % magnetic bias y
        function oy = Omegay(i)
            % N is the number of points
            global N;
            global r0; global l;
            global Omega_t;
            oy = 0;
            if(i<=N&&i>0)
                r = Hmatrix.r_int(i);
                oy = 1/2*(tanh((abs(r)-r0)/l) + 1)*Omega_t;
            end
        end
        
        
       
%--------------------------------------------------------------------------
        %the big eigen matrix
        function BM = BigEigenMatrix(m,kz)
            global N;
            global ref;
            global BC;
            BM = zeros(N*9,N*9);
            for i = 1:N
                im = (i-2)*9+1;
                ii = (i-1)*9+1;
                ip = (i  )*9+1;
                
                Hm = Hmatrix.Hminus(m,kz,i);
                HH = Hmatrix.HH(m,kz,Hmatrix.omegap(i),Hmatrix.Omegay(i),i);
                Hp = Hmatrix.Hplus(m,kz,i);
                
                % BC
                if(i==1); HH = (HH+HH*ref); end
                
                if(i==N && BC==1)
                   HH=HH+2*Hp;
                   Hm=Hm-Hp;
                end
                
                % assign
                if(i>1); BM(ii:ii+8,im:im+8) = Hm; end
                         BM(ii:ii+8,ii:ii+8) = HH;
                if(i<N); BM(ii:ii+8,ip:ip+8) = Hp; end 
                

            end
        end
        
        
%--------------------------------------------------------------------------
        function d = Det_M(m,kz,omega)
            d = det(Hmatrix.BigMatrix(m,kz,omega));
        end
        
        function ev = EigTopo(m,kz)
            ev = real(eig(Hmatrix.BigEigenMatrix(m,kz)));
            ev=ev(ev<=1.2&ev>=0.01);
        end
        %-----------------------------------
    end
end

