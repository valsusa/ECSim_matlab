%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab implementing the ECsim method descirbed in:
% Lapenta G, Exactly Energy Conserving Implicit Moment Particle in Cell Formulation
% https://arxiv.org/abs/1602.06326
% Written by G. Lapenta
% First implementation, December, 2015
% License:  GNU LESSER GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
% Copyright: KU Leuven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

graphics=1

IMM=0; ninner=3;
explicit = 0 %it does not work
groot=1;

Nx=128;
L=2*pi*groot;

NT=500;
NTOUT=25;
dt=.25*groot*10;
dt=.125*groot;
qom=-1;

th=0.5;

Np=10000;

dx=L/Nx;

xv=linspace(0,L,Nx+1)';
xc=.5*(xv(2:end)+xv(1:end-1)); 
Vx=dx*ones(Nx,1);
Ex0=zeros(Nx,1);
Ey0=zeros(Nx,1);
Ez0=zeros(Nx,1);
mode=3;

Bx=ones(Nx,1)/10;
By=zeros(Nx,1);
%By=cos(xv(1:Nx));
Bz=zeros(Nx,1);
%Bz=sin(xv(1:Nx));

Byc=zeros(Nx,1);
%Byc=zeros(Nx,1);
%Bzc=sin(xc(1:Nx));
Bzc=zeros(Nx,1);

%Ev0=+.03*sin(2*pi*mode/L*xv(1:Nx));

curlBv_y=zeros(Nx,1);
curlBv_z=zeros(Nx,1);
curlEc_y=zeros(Nx,1);
curlEc_z=zeros(Nx,1);


%Sampling in x,y cooridnates
VTx=.01;
VTy=.01;
VTz=.01;
V0=.2;V0=.8/4;
V1=.1*VTy*0;

%xp=L*rand(Np,1);
xp=linspace(0,L-L/Np,Np)';
up=VTx*randn(Np,1);
vp=VTy*randn(Np,1);
wp=VTz*randn(Np,1);
frac1=zeros(Np,1);
sigma_x=.01;

pm=[1:Np]';
pm=1-2*mod(pm,2);
%pm=mod(pm,2);
vp=vp+pm.*V0;

%up=up(randperm(Np));
%vp=vp(randperm(Np));
%wp=wp(randperm(Np));


%xp=xp+XP1*(L/Np)*sin(2*pi*xp/L*mode);
vp=vp+V1*sin(2*pi*xp/L*mode);


%computing overlap
ix=1+floor(xp/dx); % cell of the particle, first cell is cell 1, first node is 1 last node Nx+1
frac1 = 1-(xp/dx-ix+1);
ix2=mod(ix,Nx)+1;

M=zeros(Nx,Nx);

for ip=1:Np
    M(ix(ip),ix(ip))=M(ix(ip),ix(ip))+frac1(ip).^2/2; %divided by 2 to reflect the matrix to symmetrize
    M(ix2(ip),ix(ip))=M(ix2(ip),ix(ip))+frac1(ip).*(1-frac1(ip));
    M(ix2(ip),ix2(ip))=M(ix2(ip),ix2(ip))+(1-frac1(ip)).^2/2;
end    
M=M+M';

rhotarget=exp(-(xv(1:Nx)/L).^2/sigma_x^2);
rhotarget=-1;
rhotarget=rhotarget.*Vx;
rhotildeV=M\rhotarget;
for ip=1:Np
qp(ip)=rhotildeV(ix(ip))*frac1(ip)+rhotildeV(ix2(ip))*(1-frac1(ip));
end



rho=zeros(Nx,1);
for ip=1:Np
rho(ix(ip))=rho(ix(ip))+frac1(ip)*qp(ip);
rho(ix2(ip))=rho(ix2(ip))+(1-frac1(ip))*qp(ip);   
end  
rho=rho./Vx;
rhotarget=rhotarget./Vx;
figure(4)
plot(xv(1:Nx),rho,xv(1:Nx),rhotarget,'o')

% assign initial velocity

uuu=[]

histEnergy = [];
histEnergyB = [];
histEnergyE = [];
histEnergyK = [];
histMomentum = [];
histFmode = [];

spettroEx=[];
spettroEy=[];
spettroEz=[];
spettroBy=[];
spettroBz=[];
time_matrix = 0;
time_start = clock();
for it=1:NT

ix=1+floor(xp/dx); % cell of the particle, first cell is cell 1, first node is 1 last node Nx+1
frac1 = 1-(xp/dx-ix+1);
ix2=mod(ix,Nx)+1;

By(2:Nx)=.5*(Byc(2:Nx)+Byc(1:Nx-1)); By(1)=.5*(Byc(Nx)+Byc(1));
Bz(2:Nx)=.5*(Bzc(2:Nx)+Bzc(1:Nx-1)); Bz(1)=.5*(Bzc(Nx)+Bzc(1));


for ip=1:Np
    Bxp(ip)=Bx(ix(ip))*frac1(ip)+Bx(ix2(ip))*(1-frac1(ip));
    Byp(ip)=By(ix(ip))*frac1(ip)+By(ix2(ip))*(1-frac1(ip));
    Bzp(ip)=Bz(ix(ip))*frac1(ip)+Bz(ix2(ip))*(1-frac1(ip));
    alphap(ip,1:3,1:3)=alpha(qom*dt/2,Bxp(ip),Byp(ip),Bzp(ip));
end

if(explicit) 
    Jx0=zeros(Nx,1);
    Jx0(ix)=Jx0(ix)+frac1.*qp'.*up;
    Jx0(ix2)=Jx0(ix2)+(1-frac1).*qp'.*up; 
    Jy0=zeros(Nx,1);
    Jy0(ix)=Jy0(ix)+frac1.*qp'.*vp;
    Jy0(ix2)=Jy0(ix2)+(1-frac1).*qp'.*vp; 
    Jz0=zeros(Nx,1);
    Jz0(ix)=Jz0(ix)+frac1.*qp'.*wp;
    Jz0(ix2)=Jz0(ix2)+(1-frac1).*qp'.*wp; 
    Jx0=Jx0./dx;
    Jy0=Jy0./dx;
    Jz0=Jz0./dx;
else    
Jx0=zeros(Nx,1);
for ip=1:Np
    uphat=up(ip)*alphap(ip,1,1)+vp(ip)*alphap(ip,1,2)+wp(ip)*alphap(ip,1,3);
    Jx0(ix(ip))=Jx0(ix(ip))+frac1(ip)*qp(ip)*uphat;
    Jx0(ix2(ip))=Jx0(ix2(ip))+(1-frac1(ip))*qp(ip)*uphat;   
end
Jx0=Jx0./dx;

Jy0=zeros(Nx,1);
for ip=1:Np
    uphat=up(ip)*alphap(ip,2,1)+vp(ip)*alphap(ip,2,2)+wp(ip)*alphap(ip,2,3);
    Jy0(ix(ip))=Jy0(ix(ip))+frac1(ip)*qp(ip)*uphat;
    Jy0(ix2(ip))=Jy0(ix2(ip))+(1-frac1(ip))*qp(ip)*uphat;   
end
Jy0=Jy0./dx;

Jz0=zeros(Nx,1);
for ip=1:Np
    uphat=up(ip)*alphap(ip,3,1)+vp(ip)*alphap(ip,3,2)+wp(ip)*alphap(ip,3,3);
    Jz0(ix(ip))=Jz0(ix(ip))+frac1(ip)*qp(ip)*uphat;
    Jz0(ix2(ip))=Jz0(ix2(ip))+(1-frac1(ip))*qp(ip)*uphat;   
end
Jz0=Jz0./dx;
end 

% Jx0=smooth(Jx0);
% Jy0=smooth(Jy0);
% Jz0=smooth(Jz0);


%Bxc=.5*(Bx(2:Nx)+Bx(1:Nx-1)); Bxc(Nx)=.5*(Bx(Nx)+Bx(1));
%Bx is zero or assigned in 1D (with variable x)

curlBv_y(2:Nx)=-(Bzc(2:Nx)-Bzc(1:Nx-1))/dx;curlBv_y(1)=-(Bzc(1)-Bzc(Nx))/dx;
curlBv_z(2:Nx)=(Byc(2:Nx)-Byc(1:Nx-1))/dx;curlBv_z(1)=(Byc(1)-Byc(Nx))/dx;

curlEc_y=-(Ez0(2:Nx)-Ez0(1:Nx-1))/dx;curlEc_y(Nx)=-(Ez0(1)-Ez0(Nx))/dx;
curlEc_z=(Ey0(2:Nx)-Ey0(1:Nx-1))/dx;curlEc_z(Nx)=(Ey0(1)-Ey0(Nx))/dx;

curlBv_y=curlBv_y*0;
curlBv_z=curlBv_z*0;
curlEc_y=curlEc_y*0;
curlEc_z=curlEc_z*0;

time_matrix_in  = clock();
%[sum(qp'.*up) sum(J0*dx)]
M=zeros(Nx,Nx,3,3);

if(explicit)
    %nothing
else    
for i=1:3
    for j=1:3
for ip=1:Np
   % alphap(ip,1:3,1:3)=alpha(qom*dt/2,Bx,By,Bz);
   %M(ix(ip),ix(ip),i,j)=(M(ix(ip),ix(ip),i,j))+.5*frac1(ip).^2*qp(ip)*(alphap(ip,i,j)); %divided by 2 to reflect the matrix to symmetrize
   
    M(ix(ip),ix(ip),i,j)=(M(ix(ip),ix(ip),i,j))+.5*frac1(ip).^2*qp(ip)*(alphap(ip,i,j)); %divided by 2 to reflect the matrix to symmetrize
    M(ix2(ip),ix(ip),i,j)=(M(ix2(ip),ix(ip),i,j))+frac1(ip).*(1-frac1(ip))*qp(ip).*(alphap(ip,i,j));
    M(ix2(ip),ix2(ip),i,j)=(M(ix2(ip),ix2(ip),i,j))+.5*(1-frac1(ip)).^2*qp(ip).*(alphap(ip,i,j));
end   
M(:,:,i,j)=M(:,:,i,j)+M(:,:,i,j)';
if(IMM) 
    M(:,:,i,j)=diag(sum(M(:,:,i,j))); %multiplied for 0 is explicit
end    
    end
end
end

time_matrix_out = clock();
time_matrix=time_matrix + etime(time_matrix_out,time_matrix_in);
    
% 
% M11=zeros(Nx,Nx);
% for ip=1:Np
%    % alphap(ip,1:3,1:3)=alpha(qom*dt/2,Bx,By,Bz);
%        M11(ix(ip),ix(ip))=(M11(ix(ip),ix(ip)))+.5*frac1(ip).^2*qp(ip); %divided by 2 to reflect the matrix to symmetrize
%     M11(ix2(ip),ix(ip))=(M11(ix2(ip),ix(ip)))+frac1(ip).*(1-frac1(ip))*qp(ip);
%     M11(ix2(ip),ix2(ip))=(M11(ix2(ip),ix2(ip)))+.5*(1-frac1(ip)).^2*qp(ip);
% end
un=ones(Nx,1);
AmpereX=spdiags([un],[0],Nx,Nx)+qom*dt^2*th/2*(M(:,:,1,1))/dx;

AmpereY=spdiags([un],[0],Nx,Nx)+qom*dt^2*th/2*(M(:,:,2,2))/dx;

AmpereZ=spdiags([un],[0],Nx,Nx)+qom*dt^2*th/2*(M(:,:,3,3))/dx;

Derv=spdiags([-un un],[-1 0],Nx,Nx); Derv(1,Nx)=-1; Derv=Derv/dx;

Derc=spdiags([-un un],[0 1],Nx,Nx); Derc(Nx,1)=1; Derc=Derc/dx;

Faraday=spdiags([un],[0],Nx,Nx);

bKrylov=[Ex0-Jx0*dt*th; Ey0-(Jy0-0*curlBv_y)*dt*th; Ez0-(Jz0-0*curlBv_z)*dt*th; Byc; Bzc];

% Ampere=[AmpereX qom*dt^2*th/2*(M(:,:,1,2))/dx qom*dt^2*th/2*(M(:,:,1,3))/dx;
%     qom*dt^2*th/2*(M(:,:,2,1))/dx AmpereY  qom*dt^2*th/2*(M(:,:,2,3))/dx;
%     qom*dt^2*th/2*(M(:,:,3,1))/dx qom*dt^2*th/2*(M(:,:,3,2))/dx AmpereZ];


Maxwell=[AmpereX qom*dt^2*th/2*(M(:,:,1,2))/dx qom*dt^2*th/2*(M(:,:,1,3))/dx  sparse(zeros(Nx,Nx)) sparse(zeros(Nx,Nx));
    qom*dt^2*th/2*(M(:,:,2,1))/dx AmpereY  qom*dt^2*th/2*(M(:,:,2,3))/dx sparse(zeros(Nx,Nx))   +Derv*dt*th;
    qom*dt^2*th/2*(M(:,:,3,1))/dx qom*dt^2*th/2*(M(:,:,3,2))/dx AmpereZ  -Derv*dt*th  sparse(zeros(Nx,Nx))
    sparse(zeros(Nx,Nx)) sparse(zeros(Nx,Nx))   -Derc*dt*th Faraday sparse(zeros(Nx,Nx)); 
    sparse(zeros(Nx,Nx)) +Derc*dt*th sparse(zeros(Nx,Nx))   sparse(zeros(Nx,Nx)) Faraday];

xKrylov=Maxwell\bKrylov;

Ex12=xKrylov(1:Nx);
Ey12=xKrylov(Nx+1:2*Nx);
Ez12=xKrylov(2*Nx+1:3*Nx);
Byc=(xKrylov(3*Nx+1:4*Nx)-Byc*(1-th))/th;
Bzc=(xKrylov(4*Nx+1:5*Nx)-Bzc*(1-th))/th;


%
% GMRES Below works abotu just as well as the direct solver. Ampere is very
% diagonally dominant and GMRES works very well because of the condition
% number. However the overall energy conservation is not as perfect becasue
% of the tolerance in the iterative Krylov solver.
%
%Ev12=gmres(Ampere,Ev0-J0*dt/2,10,1e-7,100);


%Ex12=smooth(Ex12);
%Ey12=smooth(Ey12);
%Ez12=smooth(Ez12);

Ex0=(Ex12-Ex0*(1-th))/th;
Ey0=(Ey12-Ey0*(1-th))/th;
Ez0=(Ez12-Ez0*(1-th))/th;


if(IMM) 
    
    upold=up;
    vpold=vp;
    wpold=wp; 
    for inner=1:ninner
for ip=1:Np

up(ip)=upold(ip)*(2*alphap(ip,1,1)-1)+vpold(ip)*(2*alphap(ip,1,2))+ wpold(ip)*(2*alphap(ip,1,3))+ ...
dt*(Ex12(ix(ip))*frac1(ip)+Ex12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,1,1) + ...
dt*(Ey12(ix(ip))*frac1(ip)+Ey12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,1,2) + ...
dt*(Ez12(ix(ip))*frac1(ip)+Ez12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,1,3);

vp(ip)=upold(ip)*(2*alphap(ip,2,1))+vpold(ip)*(2*alphap(ip,2,2)-1)+ wpold(ip)*(2*alphap(ip,2,3))+...
dt*(Ex12(ix(ip))*frac1(ip)+Ex12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,2,1) + ...
dt*(Ey12(ix(ip))*frac1(ip)+Ey12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,2,2) + ...
dt*(Ez12(ix(ip))*frac1(ip)+Ez12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,2,3);

wp(ip)=upold(ip)*(2*alphap(ip,3,1))+vpold(ip)*(2*alphap(ip,3,2))+ wpold(ip)*(2*alphap(ip,3,3)-1)+...
dt*(Ex12(ix(ip))*frac1(ip)+Ex12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,3,1) + ...
dt*(Ey12(ix(ip))*frac1(ip)+Ey12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,3,2) + ...
dt*(Ez12(ix(ip))*frac1(ip)+Ez12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,3,3);


xpbar=xp(ip)+(up(ip) +upold(ip))/4*dt;
xpbar=mod(xpbar,L);

ix(ip)=1+floor(xpbar/dx); % cell of the particle, first cell is cell 1, first node is 1 last node Nx+1
frac1(ip) = 1-(xpbar/dx-ix(ip)+1);
ix2(ip)=mod(ix(ip),Nx)+1;

    Bxp(ip)=Bx(ix(ip))*frac1(ip)+Bx(ix2(ip))*(1-frac1(ip));
    Byp(ip)=By(ix(ip))*frac1(ip)+By(ix2(ip))*(1-frac1(ip));
    Bzp(ip)=Bz(ix(ip))*frac1(ip)+Bz(ix2(ip))*(1-frac1(ip));
    alphap(ip,1:3,1:3)=alpha(qom*dt/2,Bxp(ip),Byp(ip),Bzp(ip));
end
    end

xp=xp+(up +upold)/2*dt;

else
for ip=1:Np
    upold=up(ip);
    vpold=vp(ip);
    wpold=wp(ip);
up(ip)=upold*(2*alphap(ip,1,1)-1)+vpold*(2*alphap(ip,1,2))+ wp(ip)*(2*alphap(ip,1,3))+ ...
dt*(Ex12(ix(ip))*frac1(ip)+Ex12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,1,1) + ...
dt*(Ey12(ix(ip))*frac1(ip)+Ey12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,1,2) + ...
dt*(Ez12(ix(ip))*frac1(ip)+Ez12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,1,3);

vp(ip)=upold*(2*alphap(ip,2,1))+vpold*(2*alphap(ip,2,2)-1)+ wp(ip)*(2*alphap(ip,2,3))+...
dt*(Ex12(ix(ip))*frac1(ip)+Ex12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,2,1) + ...
dt*(Ey12(ix(ip))*frac1(ip)+Ey12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,2,2) + ...
dt*(Ez12(ix(ip))*frac1(ip)+Ez12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,2,3);

wp(ip)=upold*(2*alphap(ip,3,1))+vpold*(2*alphap(ip,3,2))+ wp(ip)*(2*alphap(ip,3,3)-1)+...
dt*(Ex12(ix(ip))*frac1(ip)+Ex12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,3,1) + ...
dt*(Ey12(ix(ip))*frac1(ip)+Ey12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,3,2) + ...
dt*(Ez12(ix(ip))*frac1(ip)+Ez12(ix2(ip))*(1-frac1(ip)))*qom*alphap(ip,3,3);


%xp(ip)=xp(ip)+(up(ip) +upold)/2*dt;
xp(ip)=xp(ip)+up(ip)*dt;
end

end
   Ek = 0.5*sum(qp'.*up.^2/qom+qp'.*vp.^2/qom+qp'.*wp.^2/qom);
   Ee = 0.5*sum(Ex0.^2+Ey0.^2+Ez0.^2)*dx;
   Eb = 0.5*sum(Byc.^2+Bzc.^2)*dx;
   Etot = Ek + Ee + Eb;
   histEnergy = [histEnergy Etot];
   histEnergyE = [histEnergyE Ee];
   histEnergyB = [histEnergyB Eb];
   histEnergyK = [histEnergyK Ek];
   histMomentum = [histMomentum qp*up];

   spettroEx= [spettroEx Ex0];
   spettroEy= [spettroEy Ey0];
   spettroEz= [spettroEz Ez0];
   spettroBy= [spettroBy Byc];
   spettroBz= [spettroBz Bzc];

   uuu=[uuu up];
   f=fft(Bzc);
   histFmode = [histFmode abs(f(4))];

xp=mod(xp,L);

if(mod(it,round(NT/NTOUT))==0&graphics)
    figure(1)
subplot(3,3,1)
 ii=vp<0;xxm=xp(ii);uum=up(ii);wwm=wp(ii);
 ii=vp>0;xxp=xp(ii);uup=up(ii);wwp=wp(ii);
    
    plot(xxm,uum,'r.',xxp,uup,'b.')
  
%plot(xp,up,'.')
xlabel('x')
ylabel('v_x')
xlim([0 L])
title(['\omega_{pe}t = ' num2str(it*dt) '   CFL = ' num2str(mean(abs(up))*dt/dx) '    Dx/\lambda_{De}= ' num2str(dx/VTx) ])
subplot(3,3,2)
plot(xp,vp,'.')
xlabel('x')
ylabel('v_y')
xlim([0 L])
subplot(3,3,3)
plot(up,vp,'.')
xlabel('v_x')
ylabel('v_y')
subplot(3,3,4)
semilogy(dt*(1:it),histEnergyB,dt*(1:it),histEnergyE)
subplot(3,3,5)
plot(xc,Byc,xc,Bzc)
xlim([0 L])
subplot(3,3,6)
plot(xv(1:Nx),Ex0,xv(1:Nx),Ey0,xv(1:Nx),Ez0)
xlim([0 L])
subplot(3,3,7)
plot(1:it,histEnergy-histEnergy(1))
subplot(3,3,8)
plot(1:it,histEnergyB,1:it,histEnergyE,1:it,histEnergyK,1:it,histEnergy)
subplot(3,3,9)
plot(1:it,histMomentum/Np./mean(abs(qp'.*up)))
axis tight
  
figure(2)
subplot(2,2,1)
hist(up,100)
subplot(2,2,2)
hist(vp,100)
subplot(2,2,3)
hist(wp,100)
pause(.1) 
end

%rho=zeros(Nx,1);
%for ip=1:Np
%rho(ix(ip))=rho(ix(ip))+frac1(ip)*qp(ip);
%rho(ix2(ip))=rho(ix2(ip))+(1-frac1(ip))*qp(ip);   
%end  
%rho=rho./Vx;
end
% timing
time_end = clock();
time_elapsed = etime(time_end,time_start)

t=dt*(1:it);
gamma=.35;
figure(100)
subplot(2,3,1)
plot(xp,up,'.')
xlim([0 L])
title(['\omega_{pe}t = ' num2str(it*dt) '   CFL = ' num2str(mean(abs(up))*dt/dx) '    Dx/\lambda_{De}= ' num2str(dx/VTx) ])
subplot(2,3,2)
plot(xp,vp,'.')
xlim([0 L])
subplot(2,3,3)
semilogy(t,histEnergyE,t(1:end/2))
subplot(2,3,4)
plot(1:it,histEnergy-histEnergy(1))
subplot(2,3,5)
plot(1:it,histEnergyE,1:it,histEnergyB,1:it,histEnergyK,1:it,histEnergy)
subplot(2,3,6)
plot(1:it,histMomentum/Np./mean(abs(qp'.*up)))
axis tight
  
