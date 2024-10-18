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