function [reflected_field,transmitted_field,internal_field_forward,internal_field_backward] = BPM3(rs,ts,coordinates,constant,Efin,Ef,Eb,verbose)


rs = gpuArray(single(rs));
ts = gpuArray(single(ts));
kernal = gpuArray(single(exp(1i*2*pi*coordinates.Uz(constant)*coordinates.dz)));

Ef = zeros(size(rs,1),size(rs,2),size(Efin,3),size(rs,3),'single','gpuArray');
Ef(:,:,:,1) = gpuArray(single(Efin));
for k = 2:size(rs,3)
    Ef(:,:,:,k) = ifft2(fft2(ts(:,:,k-1).*Ef(:,:,:,k-1)).*kernal);
end



Eb = zeros(size(Ef),'single','gpuArray');    
for k = size(rs,3):-1:2
    Eb(:,:,:,k-1) = ifft2(fft2(ts(:,:,k).*Eb(:,:,:,k)+rs(:,:,k).*Ef(:,:,:,k)).*kernal);
end

Ebin = Eb(:,:,:,end);



internal_field_forward = Ef;
internal_field_backward = Eb;

reflected_field = rs(:,:,1).*Efin + ts(:,:,1).*Eb(:,:,:,1);
transmitted_field = rs(:,:,end).*Ebin + ts(:,:,end).*Ef(:,:,:,end);


end

