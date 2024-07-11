function RI = get_tomogram_RYTOV(Eout,Ein,u_in,NA,mediumRI,wavelength,dx,num_pixel_z,z_batch,use_GPU)
            

            fftshift2 = @(x) fftshift(fftshift(x,1),2);
            ifftshift2 = @(x) ifftshift(ifftshift(x,1),2);

            field_normalized = Eout./Ein;

            num_angle=size(field_normalized,3);

            f_dx=u_in(:,1);
            f_dy=u_in(:,2);
                       
            k0_x=(1/dx/size(Eout,1)).*f_dx;
            k0_y=(1/dx/size(Eout,1)).*f_dy;
            k0_z=real(sqrt((mediumRI/wavelength)^2-(k0_x).^2-(k0_y).^2));
            z=(1:num_pixel_z)-ceil((num_pixel_z+1)/2);
            z=dx*reshape(z,1,1,[]);

            kx=(1:size(Eout,1))-ceil((size(Eout,1)+1)/2);
            kx = kx * (1/dx/size(Eout,1));
            ky = kx';
            
            Kzi=single(zeros(size(Eout,1),size(Eout,1),num_angle,'single'));
            NAmask=single(zeros(size(Eout,1),size(Eout,1),num_angle,'single'));
            Weight=single(zeros(size(Eout,1),size(Eout,1),1,num_angle,'single'));
            Count=single(zeros(size(Eout,1),size(Eout,1),num_pixel_z,'single'));
            

            xx=zeros(size(Eout,1),size(Eout,1));
            yy=xx;
            xx=xx+(1:size(Eout,1));
            yy=yy+(1:size(Eout,1))';


            for kk= 1 :num_angle
                kz = real(sqrt((mediumRI/wavelength)^2-(kx+k0_x(kk)).^2-(ky+k0_y(kk)).^2));
                Kz = kz-k0_z(kk);
                Kzi(:,:,kk)=Kz;
                NAmask(:,:,kk) = circshift(NA,[-round(k0_y(kk)/(1/dx/size(Eout,1))) -round(k0_x(kk)/(1/dx/size(Eout,1)))]);
                zind = round(Kz/(1/dx/num_pixel_z))+ceil((num_pixel_z+1)/2);
                zind(zind<1)=1;zind(zind>num_pixel_z)=num_pixel_z;
                zind = zind(NAmask(:,:,kk)==1);
                xind = xx(NAmask(:,:,kk)==1);
                yind = yy(NAmask(:,:,kk)==1);
                ind = sub2ind(size(Count),yind,xind,zind);
                Count(ind) = Count(ind) + 1;

            end

            Count(Count~=0)=1./Count(Count~=0);

            for kk= 1:num_angle
                zind = round(Kzi(:,:,kk)/(1/dx/num_pixel_z))+ceil((num_pixel_z+1)/2);
                zind(zind<1)=1;zind(zind>num_pixel_z)=num_pixel_z;
                zind = zind(NAmask(:,:,kk)==1);
                xind = xx(NAmask(:,:,kk)==1);
                yind = yy(NAmask(:,:,kk)==1);
                ind = sub2ind(size(Count),yind,xind,zind);
                tmp=zeros(size(Eout,1),size(Eout,1));
                tmp(NAmask(:,:,kk)==1)=Count(ind);
                Weight(:,:,1,kk)=tmp;
            end
            clear Count
            
            gradV = (single(zeros(size(Eout,1),size(Eout,1),num_pixel_z,3)));
            if use_GPU
                z = gpuArray(z);
            end
            for kk= 1 :num_angle
                fprintf("%d / %d\n",kk,num_angle)

                kz = real(sqrt((mediumRI/wavelength)^2-(kx+k0_x(kk)).^2-(ky+k0_y(kk)).^2));
                TF = 2*2*pi*kz/1i; 

                if use_GPU
                    TF = gpuArray(TF.*NAmask(:,:,kk));
                    zz=gpuArray(1:z_batch:num_pixel_z);
                    weight = gpuArray(Weight(:,:,1,kk));
                    Field = gpuArray(single(field_normalized(:,:,kk)));
                    kzi = gpuArray(Kzi(:,:,kk));
                else
                    TF = TF.*NAmask(:,:,kk);
                    zz=1:z_batch:num_pixel_z;
                    weight = Weight(:,:,1,kk);
                    Field = single(field_normalized(:,:,kk));
                    kzi = Kzi(:,:,kk);
                end

                Field = fftshift(fft2(ifftshift(Field)));
                Field = Field .* NAmask(:,:,kk);
     

                for kkk = 1:length(zz)
                    if kkk == length(zz)
                        z_batch_ind=[zz(end):num_pixel_z, 1];
                    else
                        z_batch_ind=zz(kkk):zz(kkk+1);
                    end
                    %fields propagated to the 3D FOV
                    Fieldz = fftshift2(ifft2(ifftshift2(Field.*exp(1i*2*pi*kzi.*z(z_batch_ind)))));

                    % Eq. (29)
                    Fieldz = gradient3(log(1e-4+abs(Fieldz))) + 1i*wrapToPi(gradient3(angle(Fieldz)));
                    
                    % Scattering potential (V) is combined in k-space.
                    % Weights are multiplied to deconvolve the Green's funtion (or Ewald sphere).

                    Fieldz = fftshift2(fft2(ifftshift2(Fieldz)));
                    if kkk == length(zz)
                        Fieldz(:,:,end-1,3)=0;
                    end
                    z_batch_ind = z_batch_ind(1:end-1);
                    Fieldz = Fieldz(:,:,1:end-1,:);
                    gradV(:,:,z_batch_ind,:) = gradV(:,:,z_batch_ind,:) + gather(weight.*Fieldz.*TF);
                end
            end
            clear Fieldz Weight Kzi kzi NAmask

            gradV = fftshift2(ifft2(ifftshift2((gradV))));
            

            gradV = gradV * (1/dx/num_pixel_z);
            gradV = real(gradV);

            % Inverse gradient
            V = zeros(size(gradV(:,:,:,1)));
            dct_plan = fdct2_prep( size(gradV(:,:,1,1)), gradV(1)  );
            RR=floor(0.3*size(Eout,1)):floor(0.7*size(Eout,1));
            for k = 1:size(gradV,3)
                V(:,:,k)=antigradientL2(double((squeeze(gradV(:,:,k,1:2)))),dct_plan);
                if k>1
                    V(:,:,k)=V(:,:,k)-mean(V(RR,RR,k)-V(RR,RR,k-1)-gradV(RR,RR,k-1,3),'all');
                end
            end
            POI=floor(size(V,3)/2);   
            V = V - mean(V(:,:,POI),'all');

            % Scattering potential to refractive index
            RI=mediumRI*sqrt(1+V/(2*pi*mediumRI/wavelength)^2);
end