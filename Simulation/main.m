close all
clear
restoredefaultpath;
%%
constant = Constant();
constant.mediumRI = 1.4;
constant.wavelength = 0.85;
constant.NA = 1.2;

coordinates = Coordinates();

coordinates.dx = constant.wavelength/constant.mediumRI/3;
coordinates.dz = 10;

coordinates.Nx = ceil(200/coordinates.dx);
coordinates.Nz = 15;


coordinates.update_parameters();

%% Load data of the optical system
load("system.mat")
%% Estimate the scattering parameters
L = coordinates.dz * (1:coordinates.Nz);
Tc = zeros(1,length(L));
cos_avg = zeros(1,length(L));


for k = 1:coordinates.Nz
    [reflected,transmitted] = BPM(r(:,:,1:k),t(:,:,1:k),coordinates,constant,ones(coordinates.Nx),[],[],true);
    [Tc(k),cos_avg(k)] =anisotropy(reflected,transmitted,coordinates,constant);
end



subplot(1,2,1)
plot(fit(L.',-log(cos_avg.'),"poly1"),L.',-log(cos_avg.'))
subplot(1,2,2)
plot(fit(L.',-log(Tc.'),"poly1"),L.',-log(Tc.'))


lt = 1/fit(L.',-log(cos_avg.'),"poly1").p1;
ls = 1/fit(L.',-log(Tc.'),"poly1").p1;
g = 1 - ls/lt;

fprintf("ls : %.4g    g : %.4g \n",ls,g)

%% Embed the target image to the imaging plane
zm = 100; % depth of the imaging plane and the reference mirror
target = load("target.mat").target;
pos = round(coordinates.Nx/2 - size(target,1)/2);
pos = pos : pos+size(target,1) - 1;

r(pos,pos,round(zm/coordinates.dz)+1) = r(pos,pos,round(zm/coordinates.dz)+1) .* target;

%% Simulate the time-gated reflection image
wavelength = 1./linspace(1/0.83,1/0.87,40);
spectrum = hann(length(wavelength));
axial_res = 2*log(2)/pi*(constant.wavelength^2/((wavelength(end)-wavelength(1))/2))/constant.mediumRI;
% freq = 3e8*1e6./wavelength;
% z= 0:0.5:10;
% plot(z,abs(nufft(spectrum,freq,2*z*constant.mediumRI/3e8/1e6)));


% Incident angles for conventional AO method.
angle1 = 2*pi/9*(0:8);
angle2 = 2*pi/18*(0:17)+pi/18;
angle3 = 2*pi/27*(0:26);
uin1 = [0.3*cos(angle1) 0.6*cos(angle2) 0.9*cos(angle3); 0.3*sin(angle1) 0.6*sin(angle2) 0.9*sin(angle3)]/constant.wavelength;
uin1 = round(uin1/coordinates.dux)*coordinates.dux;

% Incident angles for our method.
angle1 = pi/3*(0:5);
angle2 = pi/6*(0:11)+pi/12;
uin2 = [0.4*cos(angle1) 0.8*cos(angle2); 0.4*sin(angle1) 0.8*sin(angle2)]/constant.wavelength;
uin2 = round(uin2/coordinates.dux);
step = 1;
uin2 = [uin2 uin2+[step;0] uin2+[0;step]]*coordinates.dux;

% generate a random surface profile
surface_profile = load("surface_profile.mat").surface_profile;

% apply the surface profiles
t1  = t;
t1(:,:,1) = t(:,:,1).*exp(1i*5*surface_profile);

t2  = t;
t2(:,:,1) = t(:,:,1).*exp(1i*10*surface_profile);

% light simulation based on the BPM.
% replace t2 by t or t1 to adjust the surface profile.
ts = t2;
field_conv = ifft2(fft2(simulate3(uin1,zm,r,ts,wavelength,spectrum,constant,coordinates,[])).*coordinates.NAmask(constant));
field_our = ifft2(fft2(simulate3(uin2,zm,r,ts,wavelength,spectrum,constant,coordinates,[])).*coordinates.NAmask(constant));

field_x = field_our(:,:,size(field_our,3)/3+1:2*size(field_our,3)/3);
field_y = field_our(:,:,2*size(field_our,3)/3+1:end);
field_0 = field_our(:,:,1:size(field_our,3)/3);

%% adaptive optics
p = 3;5; %size of the window

aberration = 1;

% our method
num_iter_windowing = 8;
for k = 1:num_iter_windowing
    field_0_c = ifft2(fft2(field_0).*conj(aberration)) .* abs(coordinates.gaussian(p,0,[0 0],constant));
    field_x_c = ifft2(fft2(field_x).*conj(aberration)) .* abs(coordinates.gaussian(p,0,[0 0],constant));
    field_y_c = ifft2(fft2(field_y).*conj(aberration)) .* abs(coordinates.gaussian(p,0,[0 0],constant));
    aberration = aberration .* get_aberration(field_0_c,field_x_c,field_y_c,step,coordinates.NAmask(constant));

    I = abs(fft2(aberration));
    I = I.*(I>0.5*max(I,[],'all'));
    x = ifftshift((1:size(I,1))-ceil((size(I,1)+1)/2));
    y = x.';
    
    xc = round(sum(x.*I,'all')/sum(I,"all"));
    yc = round(sum(y.*I,'all')/sum(I,"all"));
    
    aberration = ifft2(circshift(fft2(aberration),-[yc xc]));

    if k == 1
        aberration1 = aberration;
    end
end


% conventional
field_c = field_conv .* abs(coordinates.gaussian(p,0,[0 0],constant));
aberration_conventional = get_aberration_conventional(field_c,round(uin1/coordinates.dux));

% ideal
point_source = zeros(coordinates.Nx,coordinates.Nx);
point_source(1) = 1;
[~, transmitted] = BPM(r(:,:,1:round(zm/coordinates.dz)+1),flip(ts(:,:,1:round(zm/coordinates.dz)+1),3),coordinates,constant,fftshift(point_source),[],[],false);
kernal = exp(-1i*2*pi*real(coordinates.Uz(constant))*zm);
aberration_ideal = exp(1i*angle(fft2(ifftshift(transmitted)).*kernal.*coordinates.NAmask(constant)));

figure;
subplot(2,3,1)
I = ind2rgb(1+round((angle(zeros(size(aberration1)))+pi)/(2*pi)*255),turbo(256)).*coordinates.NAmask(constant) + ~coordinates.NAmask(constant);
imshow(fftshift(fftshift(I,2),1));axis image
title("uncorrected")

subplot(2,3,2)
I = ind2rgb(1+round((angle(aberration1)+pi)/(2*pi)*255),turbo(256)).*coordinates.NAmask(constant) + ~coordinates.NAmask(constant);
imshow(fftshift(fftshift(I,2),1));axis image
title("our method")

subplot(2,3,3)
I = ind2rgb(1+round((angle(aberration)+pi)/(2*pi)*255),turbo(256)).*coordinates.NAmask(constant) + ~coordinates.NAmask(constant);
imshow(fftshift(fftshift(I,2),1));axis image
title("our method with iterative windowing")

subplot(2,3,4)
I = ind2rgb(1+round((angle(aberration_conventional)+pi)/(2*pi)*255),turbo(256)).*coordinates.NAmask(constant) + ~coordinates.NAmask(constant);
imshow(fftshift(fftshift(I,2),1));axis image
title("conventional")

subplot(2,3,5)
I = ind2rgb(1+round((angle(aberration_ideal)+pi)/(2*pi)*255),turbo(256)).*coordinates.NAmask(constant) + ~coordinates.NAmask(constant);
imshow(fftshift(fftshift(I,2),1));axis image
title("ideal")


% calculate the aberration of the incoming light using time reversal symmetry
aberration1_in = circshift(flip(flip(coordinates.NAmask(constant).*aberration1,1),2),[1 1]);
aberration_in = circshift(flip(flip(coordinates.NAmask(constant).*aberration,1),2),[1 1]);
aberration_in_conventional = circshift(flip(flip(coordinates.NAmask(constant).*aberration_conventional,1),2),[1 1]);
aberration_in_ideal = circshift(flip(flip(coordinates.NAmask(constant).*aberration_ideal,1),2),[1 1]);

% wavefront shaping using the retrieved aberration.
kernal = exp(-1i*2*pi*real(coordinates.Uz(constant))*zm);
wavefront_uncorrected = fftshift(ifft2(coordinates.NAmask(constant).*kernal));
wavefront_corrected1 = ifft2(conj(aberration1_in).*fft2(wavefront_uncorrected));
wavefront_corrected = ifft2(conj(aberration_in).*fft2(wavefront_uncorrected));
wavefront_conventional = ifft2(conj(aberration_in_conventional).*fft2(wavefront_uncorrected));
wavefront_ideal = ifft2(conj(aberration_in_ideal).*fft2(wavefront_uncorrected));

% simulate the focus formed inside the sample 
constant.wavelength = 0.85;
[~,~,internal_field,~] = BPM(r,ts,coordinates,constant,wavefront_uncorrected,[],[],false);
[~,~,internal_field_corrected1,~] = BPM(r,ts,coordinates,constant,wavefront_corrected1,[],[],false);
[~,~,internal_field_corrected,~] = BPM(r,ts,coordinates,constant,wavefront_corrected,[],[],false);
[~,~,internal_field_conventional,~] = BPM(r,ts,coordinates,constant,wavefront_conventional,[],[],false);
[~,~,internal_field_ideal,~] = BPM(r,ts,coordinates,constant,wavefront_ideal,[],[],false);

focal_plane = @(f) ifft2(fft2(f(:,:,ceil(zm/coordinates.dz)+1)).*exp(1i*2*pi*real(coordinates.Uz(constant))*(-mod(coordinates.dz-zm,coordinates.dz))));


upsamp = @(x,p) ifft2(ifftshift(padarray(fftshift(fft2(x(pos,pos))),[p p])));

focus_uncorrected = abs(upsamp(focal_plane(internal_field),900)).^2;
focus_corrected1 = abs(upsamp(focal_plane(internal_field_corrected1),900)).^2;
focus_corrected = abs(upsamp(focal_plane(internal_field_corrected),900)).^2;
focus_conventional = abs(upsamp(focal_plane(internal_field_conventional),900)).^2;
focus_ideal = abs(upsamp(focal_plane(internal_field_ideal),900)).^2;

figure;
subplot(2,3,1)

imagesc(focus_uncorrected./max(focus_uncorrected(:)));colormap gray;axis image
title("uncorrected Strehl:"+num2str(max(focus_uncorrected(:))/max(focus_ideal(:))))
subplot(2,3,2)

imagesc(focus_corrected1./max(focus_corrected1(:)));colormap gray;axis image
title("our method Strehl:"+num2str(max(focus_corrected1(:))/max(focus_ideal(:))))
subplot(2,3,3)

imagesc(focus_corrected./max(focus_corrected(:)));colormap gray;axis image
title("our method with iterative windowing Strehl:"+num2str(max(focus_corrected(:))/max(focus_ideal(:))))
subplot(2,3,4)

imagesc(focus_conventional./max(focus_conventional(:)));colormap gray;axis image
title("conventional Strehl:"+num2str(max(focus_conventional(:))/max(focus_ideal(:))))

subplot(2,3,5)
imagesc(focus_ideal./max(focus_ideal(:)));colormap gray;axis image
title("ideal Strehl:"+num2str(max(focus_ideal(:))/max(focus_ideal(:))))



%% Imaging

I_un = sum(abs(ifft2(fft2(cat(3,field_conv,field_0,field_x,field_y)).* 1)).^2,3);
I_our1 = sum(abs(ifft2(fft2(cat(3,field_conv,field_0,field_x,field_y)).* conj(aberration1))).^2,3);
I_our = sum(abs(ifft2(fft2(cat(3,field_conv,field_0,field_x,field_y)).* conj(aberration))).^2,3);
I_conv = sum(abs(ifft2(fft2(cat(3,field_conv,field_0,field_x,field_y)).* conj(aberration_conventional))).^2,3);
I_ideal= sum(abs(ifft2(fft2(cat(3,field_conv,field_0,field_x,field_y)).* conj(aberration_ideal))).^2,3);

figure;
subplot(2,3,1)

imagesc(I_un./max(I_un(:)));colormap parula;axis image
xlim([pos(1) pos(end)]);ylim([pos(1) pos(end)]);clim([0 1])
title("uncorrected")
subplot(2,3,2)

imagesc(I_our1./max(I_our1(:)));colormap parula;axis image
xlim([pos(1) pos(end)]);ylim([pos(1) pos(end)]);clim([0 1])
title("our method")

subplot(2,3,3)

imagesc(I_our./max(I_our(:)));colormap parula;axis image
xlim([pos(1) pos(end)]);ylim([pos(1) pos(end)]);clim([0 1])
title("our method with iterative windowing")
subplot(2,3,4)

imagesc(I_conv./max(I_conv(:)));colormap parula;axis image
xlim([pos(1) pos(end)]);ylim([pos(1) pos(end)]);clim([0 1])
title("conventional")

subplot(2,3,5)
imagesc(I_ideal./max(I_ideal(:)));colormap parula;axis image
xlim([pos(1) pos(end)]);ylim([pos(1) pos(end)]);clim([0 1])
title("ideal")
