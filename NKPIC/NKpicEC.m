%%%%%%%%%%%%%%%%%%%%%%%
% NK PIC
% S. Markidis and G. Lapenta
% Markidis, Stefano, and Giovanni Lapenta. "The energy conserving particle-in-cell method." 
% Journal of Computational Physics 230.18 (2011): 7037-7052.
% http://www.sciencedirect.com/science/article/pii/S0021999111003445
% September 2010
%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

global L; global dx; global NG;
global DT;

global N;
global WP; global QM; global Q; global rho_back;
global x0; global v0;
global E0;

graphics=0


% parameters
L=2*pi*10;
DT=.25*10;
NT=200;
NTOUT=25;
NG=64;
N=10000;
WP=1;
QM=-1;
V0=0.2;
VT=0.01;
% perturbation 
XP1=1*0; 
V1=.1;
mode=1;
% charge and gris parameter
Q=WP^2/(QM*N/L);
rho_back=-Q*N/L;
dx=L/NG;



% 2 Stream instability
x0=linspace(0,L-L/N,N)';
v0=VT*randn(N,1); % maxwellian


pm=[1:N]';
pm=1-2*mod(pm,2);
v0=v0+pm.*V0;

% Perturbation
v0=v0+V1*sin(2*pi*x0/L*mode);
x0=x0+XP1*(L/N)*sin(2*pi*x0/L*mode);
out=(x0<0); x0(out)=x0(out)+L;
out=(x0>=L); x0(out)=x0(out)-L;

E0 = zeros(NG,1);


% absolute and relative tollerance
tol = [1E-8, 1E-8];
% prepare xkrylov vector
histEnergy = [];
histEnergyP = [];
histEnergyK = [];
histMomentum = [];
histSolverIteration = [];
spettro=[];
time_start = clock();
for it=1:NT
   it;
   xkrylov = [v0; E0]; 
   [sol, it_hist, ierr] = nsolgm(xkrylov,'residueEC',tol);
   v_average = sol(1:N);
   v0 = 2*v_average - v0;
   x0 = x0 + v_average*DT;
   out=(x0<0); x0(out)=x0(out)+L;
   out=(x0>=L);x0(out)=x0(out)-L;
   % the new electric field
   E0 = sol((N+1):(N + NG));
   Ek = 0.5*abs(Q)*sum(v0.^2);
   Ep = 0.5*sum(E0.^2)*dx;
   Etot = Ek + Ep;
   histEnergy = [histEnergy Etot];
   histEnergyP = [histEnergyP Ep];
   histEnergyK = [histEnergyK Ek];
   histMomentum = [histMomentum sum(v0)];
   histSolverIteration = [histSolverIteration it_hist(end,2)];
if(mod(it,round(NT/NTOUT))==0&graphics)
subplot(2,3,1:3)
plot(x0,v0,'.')
axis tight
title(['\omega_{pe}t = ' num2str(it*DT) '           \Delta x/\lambda_{De} = ' num2str(dx/VT)])
subplot(2,3,4)
plot(1:it,histEnergy,1:it,histEnergyP,1:it,histEnergyK)
subplot(2,3,5)
plot(histEnergy-histEnergy(1))
subplot(2,3,6)
plot(1:it,histMomentum/N./mean(abs(v0))) 
pause(.1)   
spettro= [spettro fftshift((fft(E0)))];
end

end

% timing
time_end = clock();
time_elapsed = etime(time_end,time_start)
figure
pcolor(abs(spettro))
shading interp

figure(200)
subplot(2,3,1:3)
plot(x0,v0,'.')
axis tight
title(['\omega_{pe}t = ' num2str(it*DT) '           \Delta x/\lambda_{De} = ' num2str(dx/VT)])
subplot(2,3,4)
plot(1:it,histEnergy,1:it,histEnergyP,1:it,histEnergyK)
subplot(2,3,5)
plot(histEnergy-histEnergy(1))
subplot(2,3,6)
plot(1:it,histMomentum/N./mean(abs(v0))) 

save ECpic.txt histEnergy;
save PhaseSpaceEC.txt x0 v0;
save SolverIterations.txt histSolverIteration;
save TimeEC.txt time_elapsed;
