function [reflected_field,transmitted_field,internal_field_forward,internal_field_backward] = BPM(rs,ts,coordinates,constant,Efin,Ef,Eb,verbose)


rs = gpuArray(single(rs));
ts = gpuArray(single(ts));
kernal = gpuArray(single(exp(1i*2*pi*coordinates.Uz(constant)*coordinates.dz)));
prop = @(E) ifft2(fft2(E).*kernal);

Ef = zeros(size(rs),'single','gpuArray');
Ef(:,:,1) = gpuArray(single(Efin));
for k = 2:size(rs,3)
    Ef(:,:,k) = prop(ts(:,:,k-1).*Ef(:,:,k-1));
end



Eb = zeros(size(Ef),'single','gpuArray');    
for k = size(Ef,3):-1:2
    Eb(:,:,k-1) = prop(ts(:,:,k).*Eb(:,:,k)+rs(:,:,k).*Ef(:,:,k));
end

Ebin = Eb(:,:,end);



internal_field_forward = gather(Ef);
internal_field_backward = gather(Eb);

reflected_field = gather(rs(:,:,1).*Efin + ts(:,:,1).*Eb(:,:,1));
transmitted_field = gather(rs(:,:,end).*Ebin + ts(:,:,end).*Ef(:,:,end));


end

