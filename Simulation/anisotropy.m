function [Tc, cos_avg, std_ang] = anisotropy(r0,t0,coordinates,constant)

if (size(r0,3)>1)|(size(t0,3)>1)
    return;
end

Tc = abs(sum(t0,'all'))^2/numel(t0)^2; %collimated transmission: proportion of the unscattered light



r0 = abs(fft2(r0)).^2;
t0 = abs(fft2(t0)).^2;


% t0(1) = 0; %exclude the unscattered light when calculating the anisotropy

cos_angle = coordinates.prop_mask(constant).*coordinates.Uz(constant)/(constant.mediumRI/constant.wavelength);
X = [-cos_angle(coordinates.prop_mask(constant)); cos_angle(coordinates.prop_mask(constant))];
Y = [r0(coordinates.prop_mask(constant)); t0(coordinates.prop_mask(constant))];


cos_avg = sum(X.*Y,'all')/sum(Y,'all');
std_ang = sqrt(sum(asin(sqrt(1-X.^2)).^2.*Y,'all')/sum(Y,'all'));

end

