function aberration_out = get_aberration_AberrationMatrix(Eout0,Eout1,Eout2,u_in,du1,du2,NA,threshold,mask)

uin_x = round(u_in(:,1));
uin_y = round(u_in(:,2));

fftshift2 = @(x) fftshift(fftshift(x,1),2);
ifftshift2 = @(x) ifftshift(ifftshift(x,1),2);

window=tukeywin(size(Eout0,1),1);
window=window*window';

Eout0 = Eout0 .* window;
Eout1 = Eout1 .* window;
Eout2 = Eout2 .* window;
%% Contruct Aberration matix

% Align (Fig. b(i))
x=1:size(Eout0,2);
y=(1:size(Eout0,1))';

tilt1=exp(1i*2*pi*(du1(1).*x/length(x)+du1(2).*y/length(y)));
tilt2=exp(1i*2*pi*(du2(1).*x/length(x)+du2(2).*y/length(y)));

Eout0=fftshift2(fft2(ifftshift2(Eout0)));
Eout1=fftshift2(fft2(ifftshift2(Eout1./tilt1)));
Eout2=fftshift2(fft2(ifftshift2(Eout2./tilt2)));

% Aberration matrix (Fig. b(ii))
A_du1 = Eout1 .* conj(Eout0) ; 
A_du2 = Eout2 .* conj(Eout0) ; 

%% Reweight aberration matrix


% thresholding Eq. (11)
A_du1 = (abs(A_du1) > threshold) .* exp(1i*angle(A_du1));  
A_du2 = (abs(A_du2) > threshold) .* exp(1i*angle(A_du2));  

% masking (see Methods section)
for ii = 1:size(Eout0,3)
    mask_diag = circshift(mask,[uin_y(ii),uin_x(ii)]); % mask centered on the diagonal elements.

    %Remove wrapped parts
    if uin_y(ii) > 0
        mask_diag(1:uin_y(ii),:) = 0;
    else
        mask_diag((end+uin_y(ii)+1):end,:) = 0;
    end
        
    if uin_x(ii) > 0
        mask_diag(:,1:uin_x(ii)) = 0;
    else
        mask_diag(:,(end+uin_x(ii)+1):end) = 0;
    end
    
    A_du1(:,:,ii) = A_du1(:,:,ii) .* mask_diag .* NA;
    A_du2(:,:,ii) = A_du2(:,:,ii) .* mask_diag .* NA;
end

%% Factorization

% Eq. (12)
M_du1 = squeeze(sum(conj(A_du1).* reshape(A_du1, size(A_du1,1),size(A_du1,2),1,size(A_du1,3)),[1,2]));
M_du2 = squeeze(sum(conj(A_du2).* reshape(A_du2, size(A_du2,1),size(A_du2,2),1,size(A_du2,3)),[1,2]));

% projection power method (Eq. (14))
v1=ones(size(M_du1(:,1)));
v2=ones(size(M_du2(:,1)));

for k = 1:100
    v1=M_du1*v1;
    v2=M_du2*v2;
    v1=v1./abs(v1);
    v2=v2./abs(v2);
end

% Eq. (15)
dphi_out_du1 = sum(A_du1.*reshape(v1,1,1,[]),3);
dphi_out_du2 = sum(A_du2.*reshape(v2,1,1,[]),3);

% calculate the argument (L0 unwrap)
dphi_out_du1 = NA.*angle(dphi_out_du1.*NA./exp(1i*angle(sum(dphi_out_du1.*NA,'all'))));
dphi_out_du2 = NA.*angle(dphi_out_du2.*NA./exp(1i*angle(sum(dphi_out_du2.*NA,'all'))));
% if use_unwrap
% end
%% Calculate the phase gradient (Fig. 2b (iii))

dphi_out_dx = ((dphi_out_du1.*du2(2))-(dphi_out_du2.*du1(2)))./(du1(1).*du2(2)-du2(1).*du1(2));
dphi_out_dy = ((dphi_out_du2.*du1(1))-(dphi_out_du1.*du2(1)))./(du1(1).*du2(2)-du2(1).*du1(2));

%% Integrate the gradeint

g_real=[];
g_real(:,:,1)=dphi_out_dy;
g_real(:,:,2)=dphi_out_dx;

%ramp correction
g_real =  g_real - mean(g_real(round(size(g_real,1)/2-25):round(size(g_real,1)/2+25),round(size(g_real,2)/2-25):round(size(g_real,2)/2+25),:),[1,2]);

phi_out =  antigradient2(g_real, NA*1,0,10);
phi_out(NA)=phi_out(NA)-mean(phi_out(NA),'all');
aberration_out = exp(1i*phi_out).*NA;






