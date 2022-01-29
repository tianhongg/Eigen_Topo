% Eigen Mode Solver
% Main Driver
% Tianhong Wang
% tw474@cornell.edu
% 01-02-2022


clear all
% close all


%------------
%E_0 to E_N
%V_0 to V_N;
%B_1/2 to B_N+1/2

global L; %Domain size;
global x0; %ramp locaiton;
global l; %ramp size;
global N; %Total Grids = N+1;
global Omega_z; %Magnetic Bias
global Omega_y; %Magnetic Bias
global dx;

L  = 20*2*pi; % Full size[-L/2 to L/2], not one side   
x0 = 5*2*pi;
l  = 1*2*pi;  % n = 1/2*(tanh((x0-abs(x))/l) + 1);
N  = 200; %Total Grids = N+1;
Omega_z = 1.0;
Omega_y = 0.2*Omega_z;

% make the domain
N = ceil(N/2)*2; %round to even;
dx = L/N;
N=N+1;


%scan
Nky = 100;
Ky_max =  5;
Ky_min = -5;
dky =  (Ky_max-Ky_min)/Nky;

Band=zeros(Nky+1,300);

%kz
kz = 0.8;

for j = 1: Nky+1
    tic;
    ky = Ky_min+dky*(j-1);
    e = eig(Hmatrix.BigEigenMatrix(ky,kz));
    e=e(abs(e)<1.2&abs(e)>=0);
    Band(j,1:length(e))=e;
    toc;
    disp(j)
end


figure
plot([Ky_min:dky:Ky_max],Band(:,:),...
    'Linestyle','none','Marker','o','MarkerSize',5,'MarkerEdgeColor',[17, 138, 178]/255,...
    'MarkerFaceColor',[17, 138, 178]/255)
title(['k_z = ',num2str(kz)],'fontweight','normal')
set(gcf,'WindowStyle','normal');
set(gcf,'Position',[10 5 18 16]*30);
set(gca,'Position',[0.2 0.25 0.65 0.65]);
set(gca,'linewidth',2);
set(gca,'BoxStyle','full','Box','on')
set(gca,'fontsize',24);
set(gca,'TickDir','out')
xlabel('k_y')
ylabel('\omega')
set(gcf,'color','w')
set(gca,'color','w')
xlim([-1,1]);
ylim([0,1.2])


