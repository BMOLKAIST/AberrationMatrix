classdef Coordinates < matlab.mixin.Copyable
    %COORDINATES 이 클래스의 요약 설명 위치
    %   자세한 설명 위치
    
    properties
        Nx = [];
        Nz = [];
        x = [];
        y = [];
        z = [];

        ux = [];
        uy = [];
        uz = [];
        

        dx = [];
        dz = [];

        Lx = [];
        Lz = [];
        
        dux = [];
        duz = [];

    end

    methods
        
        function update_parameters(obj)

            obj.Lx = obj.Nx*obj.dx;
            obj.dux = 1/obj.Lx;
            obj.x = obj.dx*reshape((1:obj.Nx)-ceil((obj.Nx+1)/2),1,[]);
            obj.y = obj.dx*reshape((1:obj.Nx)-ceil((obj.Nx+1)/2),[],1);
            obj.ux = obj.dux*reshape(ifftshift((1:obj.Nx)-ceil((obj.Nx+1)/2)),1,[]);
            obj.uy = obj.dux*reshape(ifftshift((1:obj.Nx)-ceil((obj.Nx+1)/2)),[],1);

            if ~isempty(obj.Nz)
                obj.Lz = obj.Nz*obj.dz;
                obj.duz = 1/obj.Lz;
                obj.z = obj.dz*reshape((0:obj.Nz-1),1,1,[]);
                obj.uz = obj.duz*reshape(ifftshift((1:obj.Nz)-ceil((obj.Nz+1)/2)),1,1,[]);
            end
        end

        function mask = prop_mask(obj,constant)
            u0 = constant.mediumRI/constant.wavelength;
            mask = u0.^2 > obj.ux.^2 + obj.uy.^2;
        end



        function q = Uz(obj,constant)
            u0 = constant.mediumRI/constant.wavelength;
            q = real(sqrt(u0.^2-obj.ux.^2-obj.uy.^2));
        end


        function mask = NAmask(obj,constant)
            mask = (obj.ux.^2+obj.uy.^2) < (constant.NA/constant.wavelength)^2;
        end

        function f = gaussian(obj,sigma,z,uin,constant)
            uin = reshape(uin,2,1,[]);
            f = exp(-(obj.x.^2+obj.y.^2)/(2*sigma^2));
            f = f.*exp(1i*2*pi*(uin(1,1,:).*obj.x+uin(2,1,:).*obj.y));
            f = ifft2(fft2(f).*exp(1i*2*pi*real(obj.Uz(constant))*z));

        end
    end
    
end

