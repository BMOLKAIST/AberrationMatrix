function [U,V] = Wfunc(phi,f,g,p,epsilon0)
% 2D intput (phi, f, g) + scalar (p) --> 2D output (U,V)
% [1] Ghiglia, Dennis C., and Louis A. Romero. "Minimum Lp-norm two-dimensional phase unwrapping." JOSA A 13.10 (1996): 1999-2013.
% phi == unwrapped phase
% psi == wrapped phase

if nargin < 5
    epsilon0 = 0.01; % default value based on the Ref.[1]
end
    

assert( p <= 2, 'put p <= 2')
if p == 2    
    U = ones( size(phi) );
    V = ones( size(phi) );
    return
end    

[yy,xx] = size(phi);   
U = epsilon0 ./ ( abs( [diff(phi,1,1); zeros(1,xx)] - f ).^(2-p) + epsilon0 ); % Eq. (41) of Ref. [1]
V = epsilon0 ./ ( abs( [diff(phi,1,2), zeros(yy,1)] - g ).^(2-p) + epsilon0 ); % Eq. (42) of Ref. [1]









