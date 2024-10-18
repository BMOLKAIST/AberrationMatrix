function [phi_uw_iter,U_iter,V_iter,residue0] = unwrap2_Lp(phi_w,p,iterMax,tol)
% 2D input (phi, U, V), --> 2D output (Qphi)

% [1] Ghiglia, Dennis C., and Louis A. Romero. "Minimum Lp-norm two-dimensional phase unwrapping." JOSA A 13.10 (1996): 1999-2013.
% phi == unwrapped phase
% psi == wrapped phase

%%% ALgorithm Lp-norm (Section 3.B)
%%% step 1
if nargin <4
    tol = 1e-3;
end

if nargin <3
    iterMax = 1000;
end

if nargin < 2
    p = 0; % p-norm <= 2
end


%%% step 2
iterN = 0;
phi_uw_iter  = 0; % unwrapped phase
dct_plan = fdct2_prep( size(phi_w), phi_w(1) ); % Herve fct2 code
invPfunc_singleInput = @(rho) invPfunc(rho, dct_plan);

while true
    %%% step 3
    testPhase = wrap(phi_w - phi_uw_iter);
    [fTest,gTest,residueTest] = diffFunc(testPhase);
    residueN = sum( residueTest ~= 0,'all' );  %%% step 4
    %     fprintf('iter: %05d, residue: %05d\n',iterN,residueN);
    if iterN == 0 % Eq. (2 - 4)
        f0 = fTest;
        g0 = gTest;
        residue0 = residueTest;
        %         residueN0 = residueN;
    end

    if (~residueN) || p == 2 % if no residue exsits
        %%% 5: Algorithm LS (Least-Squares)
        c_iter = cfunc(fTest,gTest,1,1); % Eq. (37) with U = V = 1 or Eq. (A1).
        testPhase_uw = invPfunc_singleInput(c_iter);
        globalphase    = mean( wrap(testPhase - testPhase_uw) ,'all');
        phi_uw_iter = phi_uw_iter + testPhase_uw + globalphase; % Eq. (44)

        if nargout > 2
            [U_iter,V_iter] = Wfunc(phi_uw_iter,f0,g0,p); % Eq. (41-42)
        end

        break;

    else
        %%% 6
        [U_iter,V_iter] = Wfunc(phi_uw_iter,f0,g0,p); % Eq. (41-42)
        %%% 7
        c_iter = cfunc(f0,g0,U_iter,V_iter); % Eq. (37)
        %%% 8: Algorithm WLS (Weighted Least-Squares)
        Q_iter = @(phi) Qfunc(phi,U_iter,V_iter); % Eq.(36)

        phi_uw_iter = pcq_for_Lp_unwrap(Q_iter,c_iter,phi_uw_iter,invPfunc_singleInput,iterMax,tol);
        %%% 9
        iterN = iterN + 1;

        %%% 10
        if iterN > iterMax
            break;
        end
    end
end
end

%%% functions used
function f = wrap(f)

f = mod(f+pi,2*pi)-pi;
end

function [f,g,residue] = diffFunc(phi_w)
[yy,xx] = size(phi_w);
f = [diff(phi_w,1,1); zeros(1,xx)]; % zero at the end
g = [diff(phi_w,1,2), zeros(yy,1)]; % zero at the end

f = wrap(f);
g = wrap(g);
residue = g + circshift(f,[0,-1]) -circshift(g,[-1,0]) - f; % Eq. (4) of Ref. [1]

residue(:,end) = 0;
residue(end,:) = 0;
residue( abs(residue) <  pi ) = 0; % precision consideration
end

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

end

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
end

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

end

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
end

function [x_iter,rnorm] = pcq_for_Lp_unwrap(Q,c,x0,invPfunc,iterMax,tol)
% Q = function hangle 2D --> 2D
% c = weigthed laplacian of wrapped phase, 2D.
% x0 = initial guess, 2D.

innerProd = @(u,v) sum(u.*v,'all');
rnorm     = zeros(iterMax,1);

%%% step 1
iter   = 0;
x_iter = x0;        % 2D, solution
r_iter = c - Q(x0); % 2D, residual
% rnorm0     = norm(r_iter(:));

while true
    %%% step 2
    z_iter = invPfunc(r_iter); % 2D, P^-1 * r, preconditoined matrix applied

    %%% step 3
    iter = iter + 1;

    %%% step 4 , 5: get next search direction, p_iter
    if iter == 1
        p_iter = z_iter;
        rz     = innerProd(r_iter, z_iter);
    else
        rz_before = rz;
        rz        = innerProd(r_iter, z_iter);
        beta = rz/rz_before;
        p_iter = z_iter + beta * p_iter;
    end

    %%% step 6: get next solution (x_iter) and residual (r_iter)
    Qp  = Q(p_iter);
    pQp = innerProd(p_iter, Qp);
    alpha = rz/pQp;

    x_iter = x_iter + alpha*p_iter;
    r_iter = r_iter - alpha*Qp;

    rnorm(iter) = norm(r_iter(:));

    %%% vis
    %     subplot(121),imagesc(x_iter),axis image
    %     subplot(122),plot(log10(rnorm))
    %     drawnow;

    if (rnorm(iter) < tol) || (iter > iterMax)
        break;
    end

end

end




