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
kB=1.38e-23;
eps0=8.85e-12;

WP=sqrt(den*e^2/me/eps0)
QM=-e/me;
VT=sqrt(kB*T/me)
V0=3*VT
Deb=VT/WP

L=30*Deb

NT=300;
NG=320;
N=100000;
bc='periodic';

dx=L/NG;
DT=.1/WP;
xn=0:dx:L;


XP1=0.0;
V1=.3;
mode=1;

Q=-den*e*L/N;
rho_back=-Q*N/L;

% Plasma waves
xp=linspace(0,L-L/N,N)';
vp=VT*randn(N,1);

% % 2 Stream instability
% xp=linspace(0,L-L/N,N)';
% vp=VT*randn(N,1);
% pm=[1:N]';pm=1-2*mod(pm,2);
% vp=vp+pm.*V0;

% Perturbation
vp=vp+V1*VT*sin(2*pi*xp/L*mode);
xp=xp+XP1*(L/NG)*sin(2*pi*xp/L*mode);

un=ones(NG-1,1)*eps0;
Poisson=spdiags([un -2*un un],[-1 0 1],NG-1,NG-1);
figure
plot(xp,vp,'.')

Ep=[];
for it=1:NT
   
   % aggiornamento xp

   xp=xp+vp*DT;  
   switch(lower(bc))
      case 'periodic'
         out=(xp<0); xp(out)=xp(out)+L;
         out=(xp>=L);xp(out)=xp(out)-L;
	  case 'open-right'
	     out=(xp<0); xp(out)=-xp(out);vp(out)=-vp(out);
         in=(xp<=L); xp=xp(in);vp=vp(in);   
      otherwise
	     in=(xp>=0); xp=xp(in);vp=vp(in);
         in=(xp<=L); xp=xp(in);vp=vp(in);
   end
    		 
   % proiezione p->g

   g1=floor(xp/dx-.5)+1;g=[g1;g1+1];
   fraz1=1-abs(xp/dx-g1+.5);fraz=[fraz1;1-fraz1];	
   switch(lower(bc))
      case 'periodic'
	     out=(g<1);g(out)=g(out)+NG;
	     out=(g>NG);g(out)=g(out)-NG;
      otherwise
         out=(g<1);g(out)=1;
	     out=(g>NG);g(out)=NG;
   end
   
   N=max(size(xp));
   p=1:N;p=[p p];
   mat=sparse(p,g,fraz,N,NG);
   rho=full((Q/dx)*sum(mat))'+rho_back;

   % calcolo del campo

   Phi=Poisson\(-rho(1:NG-1)*dx^2);Phi=[Phi;0];
    switch(lower(bc))
      case 'periodic'
         Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx);
	  otherwise
	     Eg=([0; Phi(1:NG-1)]-[Phi(2:NG);0])/(2*dx);
	  end
   
   % proiezione q->p e aggiornamento velocita'

   vp=vp+mat*QM*Eg*DT;
   Ep=[Ep;norm(Phi)];
end

figure
plot(xp,vp,'.')

figure
semilogy(Ep)
title('potential')