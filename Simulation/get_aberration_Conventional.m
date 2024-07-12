function [pout, pin] = get_aberration_Conventional(field,uin)


Eout = fft2(field);


uin_x = reshape(uin(1,:),1,1,[]);
uin_y = reshape(uin(2,:),1,1,[]);


S_before = 0;
for kk = 1:size(Eout,3)
    S_before = S_before + circshift(Eout(:,:,kk),-[uin_y(kk) uin_x(kk)]);
end
S_before = abs(ifft2(S_before)).^2;
S_prev = S_before;
fprintf("Intensity after %d iterations: %.5e\n",0, mean(S_before,'all'))


pout = ones(size(Eout,1),size(Eout,2));
pin = ones(1,1,size(Eout,3));

% pin_prev = pin;
% pout_prev = pout;

for big_iter = 1:20

    %----------% %This part is for accelerating convergence. It can be removed.
%     pins = pin.*pin.*conj(pin_prev);
%     pouts = pout.*pout.*conj(pout_prev);
% 
%     pin_prev = pin;
%     pout_prev = pout;
% 
%     pin = pins;
%     pout = pouts;
    %----------%
    
    % Eq. 19
    D = conj(pout).*Eout;
    for iter = 1:10
        Dv = zeros(size(Eout,1),size(Eout,2));

        % D*v
        for kk = 1:size(Eout,3)
            Dv = Dv + circshift(D(:,:,kk).*conj(pin(kk)),-[uin_y(kk) uin_x(kk)]);
        end

        % D'*D*v
        for kk = 1:size(Eout,3)
            pin(kk) = sum(conj(D(:,:,kk)).*circshift(Dv,[uin_y(kk) uin_x(kk)]),'all');
        end
        pin = exp(-1i*angle(pin));
    end

    % Eq. 20
    D = conj(pin).*Eout;
    for iter = 1:10      
        Dv = zeros(size(Eout,1),size(Eout,2));

        % (D'v)'
        for kk = 1:size(Eout,3)
            Dv = Dv + circshift(D(:,:,kk).*conj(pout),-[uin_y(kk) uin_x(kk)]);
        end

        % (D*D'*v)'
        pout = 0;
        for kk = 1:size(Eout,3)
            pout = pout + conj(D(:,:,kk)).*circshift(Dv,[uin_y(kk) uin_x(kk)]);
        end
        pout = exp(-1i*angle(pout));
    end


    S_after = 0;
    for kk = 1:size(Eout,3)
        S_after = S_after + circshift(conj(pin(kk)).*Eout(:,:,kk).*conj(pout),-[uin_y(kk) uin_x(kk)]);
    end
    S_after = abs(ifft2(S_after)).^2;

    fprintf("Intensity after %d iterations: %.5e\n",big_iter, mean(S_after,'all'))

    
    if mean(S_prev,'all')/mean(S_after,'all') > 0.999
        break
    end

    S_prev = S_after;

end

I = abs(fft2(pout));
I = I.*(I>0.5*max(I,[],'all'));
x = ifftshift((1:size(I,1))-ceil((size(I,1)+1)/2));
y = x.';

xc = round(sum(x.*I,'all')/sum(I,"all"));
yc = round(sum(y.*I,'all')/sum(I,"all"));

pout = ifft2(circshift(fft2(pout),-[yc xc]));


end
