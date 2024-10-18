function I = stitch(imgs)
%code written by Herve


s_fft2 = @(in) fftshift(fft2(ifftshift(in)));
s_ifft2 = @(in) fftshift(ifft2(ifftshift(in)));

result_single = zeros(1024,1024,'single','gpuArray');
temp_single = zeros(1024,1024,'single','gpuArray');
temp_overlap = zeros(1024,1024,'single','gpuArray');
overlap_single = zeros(1024,1024,'single','gpuArray');
imgs=gpuArray(single(imgs));


co=make_coo(size(result_single));
cox=gpuArray(single(co{1}));
coy=gpuArray(single(co{2}));

N = 15;

list_ii=1:N;
list_jj=1:N;

for ii = list_ii
    for jj = list_jj      
        xi = max(ceil((ii-2)*1024/N),1); xf = min(ceil((ii+1)*1024/N),1024);
        yi = max(ceil((jj-2)*1024/N),1); yf = min(ceil((jj+1)*1024/N),1024);
        
        
        window=(gpuArray(single(tukeywin(length(yi:yf)))).*gpuArray(single(tukeywin(length(xi:xf))))');
        
        
        % fine tune position 
        temp_single(:)=0;
        temp_overlap(:)=0;
        temp_single(yi:yf,xi:xf)=window.*imgs(yi:yf,xi:xf,jj,ii);
        temp_overlap(yi:yf,xi:xf)=window;
        
        temp_single=s_fft2(temp_single);
        temp_overlap=s_fft2(temp_overlap);
        
        coor=temp_single.*conj(result_single); 
        coor=exp(1i.*angle(coor));
        coor=coor.*(sqrt(cox.^2+coy.^2)>0.001);
        coor=s_ifft2(coor);
        coor=coor.*(sqrt(cox.^2+coy.^2)<0.01);
        [~,m]=max(coor(:));
        [m1,m2]=ind2sub(size(coor),m);
        m1=m1-floor(size(coor,1)/2)-1;
        m2=m2-floor(size(coor,1)/2)-1;

        if m1>100
            error('stop')
        end
        result_single = result_single+temp_single.*exp(2i*pi*(m1.*cox+m2.*coy));
        overlap_single = overlap_single + temp_overlap.*exp(2i*pi*(m1.*cox+m2.*coy));

    end
    
end

result_single=gather(real(s_ifft2(result_single)));
overlap_single=gather(real(s_ifft2(overlap_single)));

I = sum(result_single,[3 4])./(sum(overlap_single,[3 4])+0.1);
I(I<0)=0;

end

function [co] = make_coo(sz)
co={};
for k=1:length(sz)
    co_crr=((1:sz(k))-floor(sz(k)/2)-1)./sz(k); 
    shape=ones(1,length(sz));
    shape(k)=length(co_crr);
    co_crr=reshape(co_crr,shape);
    co{k}=co_crr;
end
end