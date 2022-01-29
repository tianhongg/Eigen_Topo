%
% the eigen system is
% (H-)(psi-) +(H-\omega)(psi) +(H+)(psi+)=0;

classdef Hmatrix
    methods (Static)
        %-----------------------------------
        function Hm = Hminus(ky,kz)
            global dx;
            Hm = zeros(9,9);
            %-------------------
            Hm(4,8) =  kz/2;
            Hm(4,9) = -ky/2;
            Hm(5,7) = -kz/2;
            Hm(5,9) =  1i/dx;
            Hm(6,7) =  ky/2;
            Hm(6,8) = -1i/dx;
        end
        

        function Hp = Hplus(ky,kz)
            global dx;
            Hp = zeros(9,9);
            %-------------------
            Hp(7,5) = -kz/2;
            Hp(7,6) =  ky/2;
            Hp(8,4) =  kz/2;
            Hp(8,6) =  1i/dx;
            Hp(9,4) = -ky/2;
            Hp(9,5) = -1i/dx;
        end
        
        function H0 = HH(ky,kz,op,oy)
            %op = omega_p
            global dx;
            global Omega_z;
            H0 = zeros(9,9);
            %-------------------
            H0(1,2) = -1i*Omega_z; H0(1,3) = 1i*oy; H0(1,4) = 1i*op;
            H0(2,1) =  1i*Omega_z;                  H0(2,5) = 1i*op;
            H0(3,1) = -1i*oy;                       H0(3,6) = 1i*op;
            
            H0(4,1) =  -1i*op;  H0(4,8) =  kz/2;  H0(4,9) = -ky/2;
            H0(5,2) =  -1i*op;  H0(5,7) = -kz/2;  H0(5,9) = -1i/dx;
            H0(6,3) =  -1i*op;  H0(6,7) =  ky/2;  H0(6,8) =  1i/dx;
            
            H0(7,5) = -kz/2;  H0(7,6) =  ky/2;
            H0(8,4) =  kz/2;  H0(8,6) = -1i/dx;
            H0(9,4) = -ky/2;  H0(9,5) =  1i/dx;
        end
        
        
        %plasma density
        function op = omegap(i)
            % N is the number of points
            global dx; global L; global N;
            global x0; global l;
            op = 0;
            if(i<=N&&i>0)
                x = -L/2 + (i-1)*dx;
                op = 1/2*(tanh((x0-abs(x))/l) + 1);
                op = sqrt(op);
            end
        end
        
        % magnetic bias y
        function oy = Omegay(i)
            % N is the number of points
            global dx; global L; global N;
            global x0; global l;
            global Omega_y;
            oy = 0;
            if(i<=N&&i>0)
                x = -L/2 + (i-1)*dx;
                oy = 1/2*(tanh((abs(x)-x0)/l) + 1)*Omega_y;
            end
        end
        
        
%--------------------------------------------------------------------------
        %the big eigen matrix
        function BM = BigMatrix(ky,kz,omega)
            %periodic BC will be used.
            global N;
            BM = zeros(N*9,N*9);
            Hm = Hmatrix.Hminus(ky,kz);
            Hp = Hmatrix.Hplus(ky,kz);
            for i = 1:N
                im = (i-2)*9+1;
                ii = (i-1)*9+1;
                ip = (i  )*9+1;
                if(i==1); im = (N-1)*9+1; end  %BC
                if(i==N); ip = 1;         end  %BC
                BM(ii:ii+8,im:im+8) = Hm;
                BM(ii:ii+8,ii:ii+8) = Hmatrix.HH(ky,kz,Hmatrix.omegap(i),Hmatrix.Omegay(i))-eye(9)*omega;
                BM(ii:ii+8,ip:ip+8) = Hp;
            end
        end
       
%--------------------------------------------------------------------------
        %the big eigen matrix
        function BM = BigEigenMatrix(ky,kz)
            %periodic BC will be used.
            global N;
            BM = zeros(N*9,N*9);
            Hm = Hmatrix.Hminus(ky,kz);
            Hp = Hmatrix.Hplus(ky,kz);
            for i = 1:N
                im = (i-2)*9+1;
                ii = (i-1)*9+1;
                ip = (i  )*9+1;
                if(i==1); im = (N-1)*9+1; end  %BC
                if(i==N); ip = 1;         end  %BC
                BM(ii:ii+8,im:im+8) = Hm;
                BM(ii:ii+8,ii:ii+8) = Hmatrix.HH(ky,kz,Hmatrix.omegap(i),Hmatrix.Omegay(i));
                BM(ii:ii+8,ip:ip+8) = Hp;
            end
        end
        
        
%--------------------------------------------------------------------------
        function d = Det_M(ky,kz,omega)
            d = det(Hmatrix.BigMatrix(ky,kz,omega));
        end
        
        
        %-----------------------------------
    end
end

