function phi = get_aberration(field_0,field_x,field_y,step,mask)
fftshift2 = @(x) fftshift(fftshift(x,1),2);
mask = fftshift2(mask);


gradx = fftshift2(circshift(fft2(field_x),[0 -step]).*conj(fft2(field_0)));
grady = fftshift2(circshift(fft2(field_y),[-step 0]).*conj(fft2(field_0)));

gradx = gradx.*mask;
grady = grady.*mask;

AAx = squeeze(sum(reshape(gradx,size(field_0,1),size(field_0,1),[],1).*reshape(conj(gradx),size(field_0,1),size(field_0,1),1,[]),[1 2]));
AAy = squeeze(sum(reshape(grady,size(field_0,1),size(field_0,1),[],1).*reshape(conj(grady),size(field_0,1),size(field_0,1),1,[]),[1 2]));

v1=ones(size(AAx(:,1)));
v2=ones(size(AAy(:,1)));


for k = 1:100
    v1=AAx*v1;
    v2=AAy*v2;
    v1=v1./abs(v1);
    v2=v2./abs(v2);
end

gradx = sum(gradx./reshape(v1,1,1,[]),3);
grady = sum(grady./reshape(v2,1,1,[]),3);

% subplot(1,3,1)
% imagesc(angle(gradx))
% 
% subplot(1,3,2)
% imagesc(angle(grady))

g_real=[];
g_real(:,:,1)=grady;
g_real(:,:,2)=gradx;

g_real = angle(g_real./exp(1i*angle(sum(g_real.*mask,[1 2]))))/step;

phi =  ifftshift(mask.*exp(1i*antigradient2(g_real, mask*1,0,10)));



% subplot(1,3,3)
% imagesc(fftshift(angle(phi)))

end

