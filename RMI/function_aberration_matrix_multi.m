function [result, phi_in, phi_out] = function_aberration_matrix_multi(distortion_ur_refocused_t122,pos_x,pos_y)
%% prepare the reflection matrix. 
% the number 40 corresponds to the scanning interval of the tilt angle (wave vector).

R = permute(distortion_ur_refocused_t122,[2 3 1]);
R = reshape(R,1024,1024,15,15);
for ii = 1:15
    for jj = 1:15
        R(:,:,ii,jj) = circshift(fft2(R(:,:,ii,jj)),[-40*(jj-15), -40*(ii-15)]);
    end
end


%% Generate the window.

wx = 1024*pos_x;
wy = 1024*pos_y;

window_out = zeros(1024,1024);
x = -60:60;
y = x.';
f = exp(-(x.^2+y.^2)/(2*20));
window_out(1:121,1:121) = f;
window_out = abs(ifft2(window_out));
window_out = fftshift(window_out/window_out(1));
window_out = circshift(window_out,round([wy,wx])-ceil((1024+1)/2));
%
%imagesc(window_out)

R_windowed = permute(distortion_ur_refocused_t122,[2 3 1]);
R_windowed = reshape(R_windowed,1024,1024,15,15);
for ii = 1:15
    for jj = 1:15
        R_windowed(:,:,ii,jj) = circshift(fft2(R_windowed(:,:,ii,jj).*window_out),[-40*(jj-15), -40*(ii-15)]);
    end
end

%% Calculate the x-direction aberration matrix

X = R_windowed;
A = circshift(X(:,:,1:14,:),[0 -40]).*conj(X(:,:,2:15,:));

M = reshape(A,1024*1024,14*15);
M = M'*M;
x = ones(14*15,1);
for iter = 1:200
    x_prev = x;
    x = M * x;
    x = exp(1i*angle(x));
end

grad_in_x=conj(reshape(x,1,1,14,15));
grad_out_x = exp(1i*angle(sum(A.*conj(grad_in_x),[3,4])));

grad_in_x = padarray(squeeze((grad_in_x)).',[0,1],1,'pre');

%% Calculate the y-direction aberration matrix

A = circshift(X(:,:,:,1:14),[-40 0]).*conj(X(:,:,:,2:15));

M = reshape(A,1024*1024,14*15);
M = M'*M;
x = ones(14*15,1);
for iter = 1:200
    x_prev = x;
    x = M * x;
    x = exp(1i*angle(x));
end

grad_in_y=conj(reshape(x,1,1,15,14));
grad_out_y = exp(1i*angle(sum(A.*conj(grad_in_y),[3,4])));


grad_in_y = padarray(squeeze((grad_in_y)).',[1,0],1,'pre');


mask = sqrt(((1:1024)-260-5).^2  +  ((1:1024).'-260-5).^2) < 260;
mask_in = angle(grad_in_x) ~=0 & angle(grad_in_y) ~=0 ;

%% Phase unwrapping is neccesary due to the large tilt angle
[grad_out_x_u,~,~,~]  = unwrap2_Lp(angle(grad_out_x(1:526,1:526).*mask(1:526,1:526)), 0);%y
[grad_out_y_u,~,~,~]  = unwrap2_Lp(angle(grad_out_y(1:526,1:526).*mask(1:526,1:526)), 0);%x

grad_out_x(:) = 0;
grad_out_y(:) = 0;
grad_out_x(1:526,1:526) = grad_out_x_u;
grad_out_y(1:526,1:526) = grad_out_y_u;

% subplot(1,2,1)
% imagesc(squeeze((grad_out_x)));axis image;
% subplot(1,2,2)
% imagesc(squeeze((grad_out_y)));axis image;

%% Phase unwrapping is neccesary due to the large tilt angle
[grad_in_x,~,~,~]  = unwrap2_Lp(angle(grad_in_x.*mask_in), 0);%y
[grad_in_y,~,~,~]  = unwrap2_Lp(angle(grad_in_y.*mask_in), 0);%x

% subplot(1,2,1)
% imagesc(squeeze((grad_in_x)));axis image;
% subplot(1,2,2)
% imagesc(squeeze((grad_in_y)));axis image;

%% Inverse gradient

g_out=[];
g_out(:,:,1)=squeeze(grad_out_y)/40;
g_out(:,:,2)=squeeze(grad_out_x)/40;

g_in=[];
g_in(:,:,1)=-squeeze(grad_in_y);
g_in(:,:,2)=-squeeze(grad_in_x);

shift = mean(g_out(220:240,220:240,:),[1,2]);

g_out = g_out - shift;
g_in = g_in - 40*shift;

phi_out =  antigradient2((g_out), mask*1,0,10);
figure(2)
imagesc(wrapToPi(phi_out))



phi_in =  antigradient2((g_in), mask_in*1,0,10);
figure(3)
imagesc(wrapToPi(phi_in))

%% Generates the corrected image

R_corrected = zeros(1024,1024,15,15);
for ii = 1:15
    for jj = 1:15
        R_corrected(:,:,ii,jj) = circshift((mask.*exp(-1i*phi_out).*R(:,:,ii,jj)),[40*(jj-15), 40*(ii-15)]);
        R_corrected(:,:,ii,jj) = mask_in(jj,ii).*exp(-1i*phi_in(jj,ii))* R_corrected(:,:,ii,jj);
    end
end
R_corrected = reshape(R_corrected,1024,1024,[]);
R_corrected = ifft2(R_corrected);
result = abs(sum(R_corrected,3));
end
