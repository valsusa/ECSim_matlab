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