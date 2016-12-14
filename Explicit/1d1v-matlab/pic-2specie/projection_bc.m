%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explicit Matlab implementing the basic PIC method used in the Master Course:
% Introduction to Plasma Dynamics (B-KUL-G0P71B)
% https://arxiv.org/abs/1602.06326
% https://perswww.kuleuven.be/~u0052182/
% First implementation, September, 2010
% License:  GNU LESSER GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
% Copyright: KU Leuven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function g=projection_bc(g,NG,bc)
   switch(lower(bc))
      case 'periodic'
	     out=(g<1);g(out)=g(out)+NG;
	     out=(g>NG);g(out)=g(out)-NG;
      otherwise
         out=(g<1);g(out)=1;
	     out=(g>NG);g(out)=NG;
   end
