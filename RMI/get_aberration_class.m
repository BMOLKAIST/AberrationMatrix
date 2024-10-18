function [pout_tot, pin_tot, field] = get_aberration_class(field)




res = 0;
for ii = 1:15
    for jj = 1:15
        res = res + circshift(field(:,:,ii,jj),[40*(jj-15), 40*(ii-15)]);
    end
end

Iclass_prev = mean(abs(res).^2,'all');
fprintf("Iclass before correction: %.3e\n",mean(abs(res).^2,'all'))
% figure;imagesc(abs(sum(ifft2(res),3)).^2);colorbar

field = gpuArray(single(field));
pout = gpuArray(single(ones(1024,1024)));
pin = gpuArray(single(ones(1,1,15,15)));

pout_tot = pout;
pin_tot = pin;

for big_iter = 1:5
for iter = 1:10
    tmp = field.*pin;
    pin = 0;
    for ii = 1:15
        for jj = 1:15
            pin = pin + circshift(tmp(:,:,ii,jj),[40*(jj-15), 40*(ii-15)]);
        end
    end
    tmp = pin;
    pin = gpuArray(single(zeros(1,1,15,15)));
    for ii = 1:15
        for jj = 1:15
            pin(1,1,ii,jj) = sum(conj(field(:,:,ii,jj)).*circshift(tmp,-[40*(jj-15), 40*(ii-15)]),'all');
        end
    end
    pin = exp(1i*angle(pin));
end
pin_tot = pin_tot.*pin;
field=field.*pin;
res = 0;
for ii = 1:15
    for jj = 1:15
        res = res + circshift(field(:,:,ii,jj),[40*(jj-15), 40*(ii-15)]);
    end
end
fprintf("Iclass after %d iterations: %.3e\n",big_iter, mean(abs(res).^2,'all'))
for iter = 1:10
    tmp = field.*pout;
    pout = 0;
    for ii = 1:15
        for jj = 1:15
            pout = pout + circshift(tmp(:,:,ii,jj),[40*(jj-15), 40*(ii-15)]);
        end
    end
    tmp = pout;
    pout = 0;

    for ii = 1:15
        for jj = 1:15
            pout = pout + conj(field(:,:,ii,jj)).*circshift(tmp,-[40*(jj-15), 40*(ii-15)]);
        end
    end
    pout = exp(1i*angle(pout));
end
pout_tot = pout_tot.*pout;
field = field.*pout;


res = 0;
for ii = 1:15
    for jj = 1:15
        res = res + circshift(field(:,:,ii,jj),[40*(jj-15), 40*(ii-15)]);
    end
end
fprintf("Iclass after %d iterations: %.3e\n",big_iter, mean(abs(res).^2,'all'))
Iclass = mean(abs(res).^2,'all');

if Iclass_prev/Iclass > 0.999
    break
end

Iclass_prev = Iclass;

end


res = 0;
for ii = 1:15
    for jj = 1:15
        res = res + circshift(field(:,:,ii,jj),[40*(jj-15), 40*(ii-15)]);
    end
end
% fprintf("Iclass after correction: %.3e\n",mean(abs(res).^2,'all'))
% figure;imagesc(abs(sum(ifft2(res),3).^2));colorbar

pout_tot = gather(angle(conj(pout_tot)));
pin_tot = gather(angle(conj(pin_tot)));
field = gather(field);



end