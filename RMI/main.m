%% load data
load("distortion_ur_refocused_t122.mat")

%% without correction
result_uncorrected = uncorrected_image(distortion_ur_refocused_t122);


%% correction using our method

result=zeros(1024,1024,15,15);
phi_in=zeros(15,15,15,15);
phi_out=zeros(1024,1024,15,15);


N = 15;
for ii = 1:15
    for jj = 1:15
        [result(:,:,jj,ii),phi_in(:,:,jj,ii),phi_out(:,:,jj,ii)] = function_aberration_matrix_multi(distortion_ur_refocused_t122,ii/N-0.5/N,jj/N-0.5/N);
    end
    
end

%% correction using class
result_class=zeros(1024,1024,15,15);
phi_in_class=zeros(15,15,15,15);
phi_out_class=zeros(1024,1024,15,15);

N = 15;
for ii = 1:15
    for jj = 1:15
        [result_class(:,:,jj,ii),phi_in_class(:,:,jj,ii),phi_out_class(:,:,jj,ii)] = function_class_multi(distortion_ur_refocused_t122,ii/N-0.5/N,jj/N-0.5/N);
    end
    
end



%% stitching
Ilog = 10*log10(stitch(result).^2); 
Ilog = Ilog - max(Ilog(:));



Ilog_class = 10*log10(stitch(result_class).^2);
Ilog_class = Ilog_class - max(Ilog_class(:));

%% figure for images

figure(4);subplot(1,3,1);imagesc(10*log10(result_uncorrected.^2)-max(10*log10(result_uncorrected(:).^2)));clim([-60,0]);colormap gray;axis image;title("uncorrected")
figure(4);subplot(1,3,2);imagesc(Ilog);clim([-60,0]);colormap gray;axis image;title("our method")
figure(4);subplot(1,3,3);imagesc(Ilog_class);clim([-60, 0]);colormap gray;axis image;title("class")
%% figure for aberrations
figure(5);

phi_in_merged = zeros(15,15,3,15,15);
mask_in = sum(phi_in,[3 4]) ~= 0;

for ii = 1:15
    for jj = 1:15
        phi_in_merged(:,:,:,jj,ii) = ~mask_in+mask_in.*ind2rgb(1+round(((wrapToPi(phi_in(:,:,jj,ii)))+pi)/(2*pi)*255),turbo(256));
    end
end

phi_in_merged = reshape(permute(phi_in_merged,[1 4 2 5 3]),225,225,3);

subplot(2,2,1)
imshow(phi_in_merged);title("incoming, our method")

phi_out_merged = zeros(100,100,3,15,15);
mask = sum(phi_out(1:526,1:526),[3 4]) ~= 0;


for ii = 1:15
    for jj = 1:15
        phi_out_merged(:,:,:,jj,ii) = imresize(~mask+mask.*ind2rgb(1+round(((wrapToPi(phi_out(1:526,1:526,jj,ii)))+pi)/(2*pi)*255),turbo(256)),[100 100]);
    end
end

phi_out_merged = reshape(permute(phi_out_merged,[1 4 2 5 3]),1500,1500,3);

subplot(2,2,2)
imshow(phi_out_merged);title("outgoing, our method")


phi_in_merged_class = zeros(15,15,3,15,15);

for ii = 1:15
    for jj = 1:15
        phi_in_merged_class(:,:,:,jj,ii) = ~mask_in+mask_in.*ind2rgb(1+round(((wrapToPi(phi_in_class(:,:,jj,ii).'))+pi)/(2*pi)*255),turbo(256));
    end
end

phi_in_merged_class = reshape(permute(phi_in_merged_class,[1 4 2 5 3]),225,225,3);

subplot(2,2,3)
imshow(phi_in_merged_class);title("incoming, CLASS")

phi_out_merged_class = zeros(100,100,3,15,15);

for ii = 1:15
    for jj = 1:15
        phi_out_merged_class(:,:,:,jj,ii) = imresize(~mask+mask.*ind2rgb(1+round(((wrapToPi(phi_out_class(1:526,1:526,jj,ii)))+pi)/(2*pi)*255),turbo(256)),[100 100]);
    end
end

phi_out_merged_class = reshape(permute(phi_out_merged_class,[1 4 2 5 3]),1500,1500,3);

subplot(2,2,4)
imshow(phi_out_merged_class);title("outgoing, CLASS")