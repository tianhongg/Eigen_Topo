% Eigen Mode Solver
% Main Driver
% Tianhong Wang
% tw474@cornell.edu
% 01-02-2022


clear all
close all


%------------
%E_0 to E_N
%V_0 to V_N;
%B_1/2 to B_N+1/2


global L; %Domain size;
global r00; % Axis offset
global r0;%ramp locaiton;
global l; %ramp size;  % n = 1/2*(tanh((r0-abs(r))/l) + 1);
global N; %Total Grids = N+1;
global Omega_z; %Magnetic Bias
global Omega_t; %Magnetic Bias theta
global dr;
global ref;  %reflection matrix
global BC;  % 0 = PEC??? 1 = a simple open BC??

%reflection matrix
ref = eye(9);
ref(1,1)=-1;ref(2,2)=-1;
ref(4,4)=-1;ref(5,5)=-1;
ref(7,7)=-1;ref(8,8)=-1;

%------------------
L  = 10*2*pi;
r0 = 5*2*pi;
l  = 1*2*pi;
N  = 200;
Omega_z = 0.98589;
Omega_t = 0.2*Omega_z;
BC=1;

% make the domain
N = ceil(N/2)*2; %round to even;
dr = L/N;
r00=dr/2;  % don't change
N=N+1;  %Total Grids = N+1;

% m scan
Nm  = 20;

kz=-0.75;
Band=zeros(2*Nm+1,300);

f = waitbar(0,'waiting.....');
tic;
for j = 1: 2*Nm+1
    m  = -Nm-1+j;
    e = Hmatrix.EigTopo(m,kz);
    Band(j,1:length(e))=e;
    t2=toc()/j;
    waitbar(j/2/Nm,f,sprintf('[ %d ]  [TimeLapse %5.2fm ]  [Remain %5.2fm ]',j,toc()/60,(2*Nm+1-j)*t2/60))
end
close(f)

figure
plot([-Nm:Nm],Band,...
    'Linestyle','none','Marker','o','MarkerSize',5,'MarkerEdgeColor',[17, 138, 178]/255,...
    'MarkerFaceColor',[17, 138, 178]/255)
title(['k_z = ',num2str(kz)],'fontweight','normal')
set(gcf,'WindowStyle','normal');
set(gcf,'Position',[10 5 18 16]*30);
set(gca,'Position',[0.2 0.25 0.65 0.60]);
set(gca,'linewidth',2);
set(gca,'BoxStyle','full','Box','on')
set(gca,'fontsize',24);
set(gca,'TickDir','out')
xlabel('m')
ylabel('\omega')
set(gcf,'color','w')
set(gca,'color','w')
xlim([-Nm,Nm]);
ylim([0,1.2])