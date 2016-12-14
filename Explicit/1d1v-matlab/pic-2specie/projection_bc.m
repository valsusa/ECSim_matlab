function g=projection_bc(g,NG,bc)
   switch(lower(bc))
      case 'periodic'
	     out=(g<1);g(out)=g(out)+NG;
	     out=(g>NG);g(out)=g(out)-NG;
      otherwise
         out=(g<1);g(out)=1;
	     out=(g>NG);g(out)=NG;
   end
