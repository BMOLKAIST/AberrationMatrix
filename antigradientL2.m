function y = antigradientL2(x,f)
if nargin == 1
     f = fdct2_prep(size(x(:,:,1)),x(1) );
end
x(end,:,1)=0;
x(:,end,2)=0;
y=invPfunc(x(:,:,1)+x(:,:,2)-circshift(x(:,:,1), [1,0])-circshift(x(:,:,2), [0,1]),f);
end

function phi = invPfunc(rho, dct_plan)
% 2D intput (rho) --> 2D output (U,V)
% phi == unwrapped phase
% psi == wrapped phase


[yy,xx] = size(rho);
denominator = 2*(cos( (0:yy-1).'/yy*pi ) + cos( (0:xx-1)  /xx*pi ) - 2); 
denominator(1,1) = 1; % to prevent divdided by 0

rho = fdct2(rho,dct_plan);
% rho = dct(rho);
% rho= permute(rho,[2,1]);
% rho = dct(rho);
% rho= ipermute(rho,[2,1]);

%%% deconv
phi = rho./ denominator;

%%% fidct2 by Herve
phi = ifdct2(phi,dct_plan);
% phi = idct(phi);
% phi= permute(phi,[2,1]);
% phi = idct(phi);
% phi= ipermute(phi,[2,1]);
end