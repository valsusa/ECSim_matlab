%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explicit Matlab implementing the basic PIC method used in the Master Course:
% Introduction to Plasma Dynamics (B-KUL-G0P71B)
% https://arxiv.org/abs/1602.06326
% https://perswww.kuleuven.be/~u0052182/
% First implementation, September, 2010
% License:  GNU LESSER GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
% Copyright: KU Leuven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xp,vp]=impose_bc(xp,vp,L,bc)

   switch(lower(bc))
      case 'periodic'
         out=(xp<0); xp(out)=xp(out)+L;
         out=(xp>L); xp(out)=xp(out)-L;
	  case 'open-right'
	     out=(xp<0); xp(out)=-xp(out);vp(out)=-vp(out);
         in=(xp<=L); xp=xp(in);vp=vp(in);   
      otherwise
	     in=(xp>=0); xp=xp(in);vp=vp(in);
         in=(xp<=L); xp=xp(in);vp=vp(in);
   end