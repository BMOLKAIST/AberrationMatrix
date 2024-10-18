function Qphi = Qfunc(phi,U,V)
% 2D input (phi, U, V), --> 2D output (Qphi)

% [1] Ghiglia, Dennis C., and Louis A. Romero. "Minimum Lp-norm two-dimensional phase unwrapping." JOSA A 13.10 (1996): 1999-2013.
% phi == unwrapped phase
% psi == wrapped phase

shift10 = @(mat) circshift(mat, [1,0]);
shift01 = @(mat) circshift(mat, [0,1]);

%% reshape
[yy,xx] = size(phi);
% phi = reshape(phi,yy,xx); % make column vector;


%% diff gen.
diff1_phi = [diff(phi,1,1); zeros(1,xx)]; % psi_{i+1,j} - psi_{i,j}
% diffm1_psi = shift10(diffp1_psi); % psi_{i,j}   - psi_{i-1,j}

diff2_phi = [diff(phi,1,2), zeros(yy,1)]; % psi_{i,j+1} - psi_{i,j}
% diffm2_psi = shift01(diffp2_psi; % psi_{i,j}   - psi_{i,j-1}

%%
U = diff1_phi .* U;
V = diff2_phi .* V;

Qphi = U + V - shift10(U) - shift01(V); % Eq.(36) of Ref. [1]
% Qphi = reshape(Qphi,[],1); % make column vector;






