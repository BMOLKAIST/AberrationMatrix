function field_t = simulate3(uin,zm,r,t,wavelength,spectrum,constant,coordinates,wavefront)
constant = copy(constant);
field_t = zeros(coordinates.Nx,coordinates.Nx,size(uin,2));

for k = 1:length(wavelength)

    fprintf("%d / %d \n",k,length(wavelength))
    constant.wavelength = wavelength(k);
    Ein = coordinates.gaussian(10,-zm,uin,constant); % 10: size of the window_in
    if ~isempty(wavefront)
        Ein = ifft2(fft2(Ein).*wavefront);
    end
    f1 = BPM3(r,t,coordinates,constant,Ein,[],[],false);
    kernal = gpuArray(single(exp(-1i*2*pi*coordinates.Uz(constant)*zm)));
    field_t = field_t + sqrt(spectrum(k))*gather(ifft2(fft2(f1).*kernal));
end

end

