function cout = cfunc(f,g,U,V)
% 2D input (f, g, U, V), --> 2D output (cout)

% [1] Ghiglia, Dennis C., and Louis A. Romero. "Minimum Lp-norm two-dimensional phase unwrapping." JOSA A 13.10 (1996): 1999-2013.
% phi == unwrapped phase
% psi == wrapped phase

shift10 = @(mat) circshift(mat, [1,0]);
shift01 = @(mat) circshift(mat, [0,1]);

U = f .* U;
V = g .* V;

cout = U + V - shift10(U) - shift01(V); % Eq.(37) of Ref. [1]
% cout = reshape(cout,[],1); % make column vector;



