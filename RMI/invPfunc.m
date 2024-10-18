function phi = invPfunc(rho, dct_plan)
% 2D intput (rho) --> 2D output (U,V)
% [1] Ghiglia, Dennis C., and Louis A. Romero. "Minimum Lp-norm two-dimensional phase unwrapping." JOSA A 13.10 (1996): 1999-2013.
% phi == unwrapped phase
% psi == wrapped phase


[yy,xx] = size(rho);
denominator = 2*(cos( (0:yy-1).'/yy*pi ) + cos( (0:xx-1)  /xx*pi ) - 2); 
denominator(1,1) = 1; % to prevent divdided by 0

%%% fdct2 by Herve
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