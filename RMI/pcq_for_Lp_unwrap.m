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





