function alpha = alpha(beta,Bx,By,Bz )
%B=[Bx By Bz];
%I=diag([1 1 1]);
%IcB=[cross(I(1,:),B); cross(I(2,:),B); cross(I(3,:),B)];

%alpha=(I-beta*IcB+beta^2*B'*B)/(1+beta^2*B*B');

sx=Bx*beta;sy=By*beta;sz=Bz*beta;
alpha=[1+sx*sx  sz+sx*sy   -sy+sx*sz;
      -sz+sx*sy  1+sy*sy   sx+sy*sz;
      sy+sx*sz   -sx+sy*sz    1+sz*sz]/(1+sx*sx+sy*sy+sz*sz);
end

