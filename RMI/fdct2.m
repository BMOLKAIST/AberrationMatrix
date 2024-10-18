function C = fdct2(C,plan)
%ORIGINAL ALGORITHM FROM : JOHN MAKHOUL / A Fast Cosine Transform in One and Two Dimensions
%CODE BY : Herve Hugonnet
if ~plan{4}
   C=single(C);
end
sz=size(C);
if length(sz)==2
   sz=[sz(1),sz(2),1]; 
end
if size(C,3)>1
    C=reshape(C,sz(1)*sz(2),sz(3));
    C(plan{3},:)=C(:,:);
    C=reshape(C,sz(1),sz(2),sz(3));
else
    C(plan{3})=C(:);
end
reps=1;

if ~isreal(C)
    [C1,C2]=fft2_re_im(C);
    C=C1;
    reps=2;
else
    C=fft2(C);
end

for ii=1:reps
    if ii==2
        C1=C;
        C=C2;
    end
    %C = real(sqrt(1/(16*size(IMG,1)*size(IMG,2)))*(MULT1.*(MULT2.*trans+conj(MULT2).*circshift(flip(trans,2),[0,1]))));
    C =( plan{1}.*plan{2}.*(C));
    D=circshift(flip(imag(C),2),[0,1]);
    if sz(3)>1
        D(:,1,:)=0;%D(:,2);
    else
        D(:,1)=0;%D(:,2);
    end
    C= real(C)-D;
    C=real(sqrt(1/(16*sz(1)*sz(2))))*C;
    if sz(3)>1
        C(:,1,:)=sqrt(2)*C(:,1,:);
        C(1,:,:)=sqrt(1/2)*C(1,:,:);
    else
        C(:,1)=sqrt(2)*C(:,1);
        C(1,:)=sqrt(1/2)*C(1,:);
    end
end
if reps==2
    C=C+1i*C1;
end
end

function [re_fft_IMG,img_fft_IMG] = fft2_re_im(IMG)
%give the fft of the real and imaginary part separatly with only one fft call
fft_IMG=fft2(IMG);
fft_IMG_flip=fft_IMG;
fft_IMG_flip(2:end,2:end,:)=flip(flip(fft_IMG(2:end,2:end,:),1),2);
fft_IMG_flip(1,2:end,:)=flip(fft_IMG(1,2:end,:),2);
fft_IMG_flip(2:end,1,:)=flip(fft_IMG(2:end,1,:),1);

re_fft_IMG= (1/2)*(fft_IMG+conj(fft_IMG_flip));
img_fft_IMG=(-1i/2)*(fft_IMG-conj(fft_IMG_flip));
end