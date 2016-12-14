%%%%%%%%%%%%%%%%%%%%%%%
% Function residue EC needed by NK PIC
% S. Markidis and G. Lapenta
% Markidis, Stefano, and Giovanni Lapenta. "The energy conserving particle-in-cell method." 
% Journal of Computational Physics 230.18 (2011): 7037-7052.
% http://www.sciencedirect.com/science/article/pii/S0021999111003445
% September 2010
%%%%%%%%%%%%%%%%%%%%%%%

function res = residueEC(xkrylov)
global L; global dx; global NG;
global DT;

global N;
global WP; global QM; global Q; global rho_back;
global x0; global v0;
global E0;



% calculate the x at n+1/2 time level
x_average = x0 + xkrylov(1:N)*DT/2;
out=(x_average<0); x_average(out)=x_average(out)+L;
out=(x_average>=L);x_average(out)=x_average(out)-L;
% calculate the E on particle
p=1:N;p=[p p];
g1=floor(x_average/dx-.5)+1;     
g=[g1;g1+1];
fraz1=1-abs(x_average(1:N)/dx-g1+.5); 
fraz=[(fraz1);1-fraz1];	
out=(g<1);g(out)=g(out) + NG;
out=(g>NG);g(out)=g(out)- NG;
mat=sparse(p,g,fraz,N,NG);


res = zeros(N + NG,1);
% residual for the average velocity
res(1:N,1) = xkrylov(1:N) - v0 - 0.25*mat*QM*(E0 + xkrylov((N+1):(N+NG)))*DT;
% calculate the J
fraz=[(fraz1).*xkrylov(1:N);(1-fraz1).*xkrylov(1:N)];	
mat=sparse(p,g,fraz,N,NG);
J = full((Q/dx)*sum(mat))';
% residual for the electric field
res((N+1):(N+NG)) = xkrylov((N+1):(N+NG)) - E0 + J*DT;


   
   
  


