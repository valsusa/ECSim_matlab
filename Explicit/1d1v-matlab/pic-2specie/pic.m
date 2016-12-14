%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explicit Matlab implementing the basic PIC method used in the Master Course:
% Introduction to Plasma Dynamics (B-KUL-G0P71B)
% https://arxiv.org/abs/1602.06326
% https://perswww.kuleuven.be/~u0052182/
% First implementation, September, 2010
% License:  GNU LESSER GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
% Copyright: KU Leuven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

% All units are SI
%den=1e10;
%TeV=1e5;
%T=1.16e4*TeV;
den=3e6;
T=1e6;
e=1.6e-19;
me=9.1e-31;
mi=me*100;
kB=1.38e-23;
eps0=8.85e-12;

WP=sqrt(den*e^2/me/eps0)
QMe=-e/me;
QMi=e/mi;
VTe=sqrt(kB*T/me)
VTi=sqrt(kB*T/mi)
V0=3*VTe;
Deb=VTe/WP

L=50*Deb

NT=500;
NG=160;
N=20000;
bc='open-right';

dx=L/NG;
DT=.1/WP;
xc=dx/2:dx:L-dx/2;


XP1=1;
V1=0.0;
mode=1;

Qe=-den*e*L/N;
Qi=den*e*L/N;

%rho_back=-Qe*N/L;
rho_back=0;

% 2 Stream instability
xpe=linspace(0,L-L/N,N)';
vpe=VTe*randn(N,1);
%pm=[1:N]';pm=1-2*mod(pm,2);
%vpe=vpe+pm.*V0;

xpi=linspace(0,L-L/N,N)';
vpi=VTi*randn(N,1);


% Perturbation
vpe=vpe+V1*VTe*sin(2*pi*xpe/L*mode);
xpe=xpe+XP1*(L/N)*sin(2*pi*xpe/L*mode);

un=ones(NG-1,1)*eps0;
Poisson=spdiags([un -2*un un],[-1 0 1],NG-1,NG-1);
switch(lower(bc))
    case 'open-right'
      Poisson(1,1)=Poisson(1,1)+Poisson(1,2);
end   
figure(1)
plot(xpe,vpe,'b.',xpi,vpi,'r.')
pause

PP=[];
for it=1:NT
   
   % aggiornamento xp

   xpe=xpe+vpe*DT;  
   xpi=xpi+vpi*DT;     
   [xpe,vpe]=impose_bc(xpe,vpe,L,bc);
   [xpi,vpi]=impose_bc(xpi,vpi,L,bc);
    		 
   % proiezione p->g

   
   g1=floor(xpe/dx-.5)+1;g=[g1;g1+1];
   fraz1=1-abs(xpe/dx-g1+.5);fraz=[fraz1;1-fraz1];	
   g=projection_bc(g,NG,bc); 
   
   N=max(size(xpe));
   p=1:N;p=[p p];
   mate=sparse(p,g,fraz,N,NG);
   rho=full((Qe/dx)*sum(mate))'+rho_back;
   
   g1=floor(xpi/dx-.5)+1;g=[g1;g1+1];
   fraz1=1-abs(xpi/dx-g1+.5);fraz=[fraz1;1-fraz1];	
   g=projection_bc(g,NG,bc); 
   
   N=max(size(xpi));
   p=1:N;p=[p p];
   mati=sparse(p,g,fraz,N,NG);
   rho=rho+full((Qi/dx)*sum(mati))';


   % calcolo del campo

   Phi=Poisson\(-rho(1:NG-1)*dx^2);Phi=[Phi;0];
    switch(lower(bc))
      case 'periodic'
         Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx);
      case 'open-right'
         Eg=([Phi(1); Phi(1:NG-1)]-[Phi(2:NG);0])/(2*dx);
	  otherwise
	     Eg=([0; Phi(1:NG-1)]-[Phi(2:NG);0])/(2*dx);
	  end
   
   % proiezione q->p e aggiornamento velocita'

   vpe=vpe+mate*QMe*Eg*DT;
   vpi=vpi+mati*QMi*Eg*DT;

PP=[PP,mean(Phi)];   

if(mod(it,10)==0)
figure(1)
plot(xpe,vpe,'b.',xpi,vpi,'r.')

figure(2)
plot(xc,Phi)
title('potential')
if(it==10) 
    pause
end
end
end

figure(3) 
plot(PP)
