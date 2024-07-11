%% REMARK
% We provide CPU version code for those who don't have a CUDA compatible GPU.
% Therefore, it takes "much longer" than the GPU version. 
% You can simply use the GPU computing toolbox to run this code on a GPU.
close all
clear
restoredefaultpath;
%% Load data (10 um)
num_pixel_z = 200;
field_max = 1.7;
load("dataset_10um_aberrated.mat");
Ein_pf=load("dataset_10um_polymer_free.mat").Ein_pf;
Eout_pf=load("dataset_10um_polymer_free.mat").Eout_pf;
u_in_pf=load("dataset_10um_polymer_free.mat").u_in_pf;
disp("Loading done.")
%% Load data (50 um)
num_pixel_z = 600;
field_max = 2.2;
load("dataset_50um_aberrated.mat");
Ein_pf=load("dataset_50um_polymer_free.mat").Ein_pf;
Eout_pf=load("dataset_50um_polymer_free.mat").Eout_pf;
u_in_pf=load("dataset_50um_polymer_free.mat").u_in_pf;
disp("Loading done.")
%% Load data (100 um)
load("dataset_100um_aberrated.mat");
field_max = 3;
Ein_pf=load("dataset_100um_polymer_free.mat").Ein_pf;
Eout_pf=load("dataset_100um_polymer_free.mat").Eout_pf;
u_in_pf=load("dataset_100um_polymer_free.mat").u_in_pf;
disp("Loading done.")
%% Aberration matrix
Eout0=Eout(:,:,ind0);
Eout1=Eout(:,:,ind1); % tilted by du1 from Eout0
Eout2=Eout(:,:,ind2); % tilted by du2 from Eout0


if mask_radius == inf
    mask = ones(size(NA));
else
    mask = imresize(NA,mask_radius);
    mask = padarray(mask, floor([size(NA,1)-size(mask,1) size(NA,2)-size(mask,2)]/2));
    if mod(size(NA,1)-size(mask,1),2) == 1
        mask = padarray(mask,[1 0],"pre");
    end
    if mod(size(NA,2)-size(mask,2),2) == 1
        mask = padarray(mask,[0 1],"pre");
    end
end

aberration = get_aberration_AberrationMatrix(Eout0,Eout1,Eout2,u_in(ind0,:),du1,du2,NA,threshold,mask);
imagesc(angle(aberration));axis image;clim([-pi pi]);colormap turbo
%% 2 AO method (CLASS or Distortion matrix)
%Note: Conventional methods are very unstable for transmission imaging of thick samples. 
% Changing the starting solution or iteration number of the inner loop gives different results.

% Confirm that this code works for 2D samples (S) using Eout_thin.
% S = randn(size(Ein(:,:,1)))+1i*randn(size(Ein(:,:,1)));
% Eout_thin = ifft2(fft2(Ein.*S).*ifftshift(aberration));
% [aberration_conventional_test, S_before, S_after] = get_aberration_conventional(Eout_thin,u_in,NA);
% imagesc(angle(aberration_conventional_test));

[aberration_conventional, S_before, S_after] = get_aberration_conventional(Eout,u_in,NA);

figure;
subplot(1,2,1)
imagesc(S_before);axis image;colorbar;clim([0 2000]);title("S(r) befor correction")
subplot(1,2,2)
imagesc(S_after);axis image;colorbar;clim([0 600000]);title("S(r) after correction")

%% Deconvolution 
Eout_AberrationMatrix = ifft2(fft2(Eout).*ifftshift(conj(aberration)));
Eout_conventional = ifft2(fft2(Eout).*ifftshift(conj(aberration_conventional)));

%% Fields
load("colormap2D.mat")
figure;
ind = 1; %normal illumination

subplot(2,2,1)
offset = exp(1i*angle(sum(Eout_pf(:,:,ind)./Ein_pf(:,:,ind),"all")))./exp(1i*angle(sum(Eout(:,:,ind)./Ein(:,:,ind),"all")));
imshow(complex2rgb(offset.*exp(1i*2*pi*0.4)*(Eout(:,:,ind)./Ein(:,:,ind)),field_max,colormap2D))
title("Polymer added")

subplot(2,2,2)
offset = exp(1i*angle(sum(Eout_pf(:,:,ind)./Ein_pf(:,:,ind),"all")))./exp(1i*angle(sum(Eout_AberrationMatrix(:,:,ind)./Ein(:,:,ind),"all")));
imshow(complex2rgb(offset.*exp(1i*2*pi*0.4)*(Eout_AberrationMatrix(:,:,ind)./Ein(:,:,ind)),field_max,colormap2D))
title("Corrected")

subplot(2,2,3)
imshow(complex2rgb(exp(1i*2*pi*0.4)*(Eout_pf(:,:,ind)./Ein_pf(:,:,ind)),field_max,colormap2D))
title("Without Polymer")

subplot(2,2,4)
offset = exp(1i*angle(sum(Eout_pf(:,:,ind)./Ein_pf(:,:,ind),"all")))./exp(1i*angle(sum(Eout_conventional(:,:,ind)./Ein(:,:,ind),"all")));
imshow(complex2rgb(offset.*exp(1i*2*pi*0.4)*(Eout_conventional(:,:,1)./Ein(:,:,1)),2.5,colormap2D))
title("Conventional")


%% Phase Correlation
z = dx.*reshape(-70:70,1,1,[]);
ux=(1:size(Eout,1))-ceil((size(Eout,1)+1)/2);
ux = ux * (1/dx/size(Eout,1));
uy = ux';
uz = real(sqrt((mediumRI/wavelength)^2-(ux).^2-(uy).^2));

fftshift2 = @(x) fftshift(fftshift(x,1),2);
ifftshift2 = @(x) ifftshift(ifftshift(x,1),2);

%Eq. (30), using the fact that IFT_z delta(k_3D - k_medium) (z) = exp(1i*sqrt(k_medium^2-kx^2-ky^2).*z)
corr_2D_corrected = exp(1i*angle(fft2(Eout_AberrationMatrix).*conj(fft2(Eout_pf))));
corr_2D_corrected = NA.*fftshift(exp(1i*(angle(sum(corr_2D_corrected,3)))));
corr_3D_corrected = abs(fftshift2(ifft2(ifftshift2(corr_2D_corrected.*exp(1i*2*pi*uz.*z))))).^2./abs(max(ifft2(NA),[],'all')).^2;
[peak_corr, pos]=max(corr_3D_corrected(:));
[y0, x0, z0] = ind2sub(size(corr_3D_corrected),pos);

%Eq. (30)
corr_2D_aberrated = exp(1i*angle(fft2(Eout).*conj(fft2(Eout_pf))));
corr_2D_aberrated = NA.*fftshift(exp(1i*(angle(sum(corr_2D_aberrated,3)))));
corr_3D_aberrated = abs(fftshift2(ifft2(ifftshift2(corr_2D_aberrated.*exp(1i*2*pi*uz.*z))))).^2./abs(max(ifft2(NA),[],'all')).^2;
[peak_corr_aberrated, pos]=max(corr_3D_aberrated(:));
[y0_aberrated, x0_aberrated, z0_aberrated] = ind2sub(size(corr_3D_aberrated),pos);

figure;
subplot(1,2,1)
imshow(cat(3,corr_3D_corrected(y0+(-30:30),x0+(-30:30),z0)/peak_corr,    corr_3D_aberrated(y0+1+(-30:30),x0+2+(-30:30),z0_aberrated)/peak_corr_aberrated,zeros(size(corr_3D_corrected(y0+(-30:30),x0+(-30:30),z0))) ))
subplot(1,2,2)
imshow(cat(3,squeeze(corr_3D_corrected(y0,x0+(-30:30),z0+(-30:30))).'/peak_corr,  squeeze(corr_3D_aberrated(y0+1,x0+2+(-30:30),z0_aberrated+(-30:30))).'/peak_corr_aberrated , zeros(size(squeeze(corr_3D_corrected(y0,x0+(-30:30),z0+(-30:30))).'))))

%% Angular spectrum
figure;
imagesc(log(abs(fftshift(fft2(Eout(:,:,1))))));clim(max(log(abs(fftshift(fft2(Eout(:,:,1))))),[],'all')+[-6 0]);colormap gray;colorbar;axis image

%% Tomograms
z_batch = 200; %Reduce this if you are running out of GPU memory.
use_GPU = true;

%polymer added
RI_polymer = get_tomogram_RYTOV(Eout(:,:,ind0),Ein(:,:,ind0),u_in(ind0,:),NA,mediumRI,wavelength,dx,num_pixel_z,z_batch,use_GPU);

%corrected using our method
RI_aberration_matrix = get_tomogram_RYTOV(Eout_AberrationMatrix(:,:,ind0),Ein(:,:,ind0),u_in(ind0,:),NA,mediumRI,wavelength,dx,num_pixel_z,z_batch,use_GPU);

%no polymer
RI_polymer_free = get_tomogram_RYTOV(Eout_pf(:,:,ind0),Ein_pf(:,:,ind0),u_in_pf(ind0,:),NA,mediumRI,wavelength,dx,num_pixel_z,z_batch,use_GPU);

%corrected using the conventional method
RI_conventional = get_tomogram_RYTOV(Eout_conventional(:,:,ind0),Ein(:,:,ind0),u_in(ind0,:),NA,mediumRI,wavelength,dx,num_pixel_z,z_batch,use_GPU);

figure;orthosliceViewer(horzcat(RI_polymer,RI_aberration_matrix,RI_polymer_free,RI_conventional)); title('Tomograms: polymer added/corrected/without polymer/conventional');clim([1.47 1.55]);colormap inferno
