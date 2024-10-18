function C = ifdct2(IMG,plan)
%ALGORITHM INSPIRED FROM : JOHN MAKHOUL / A Fast Cosine Transform in One and Two Dimensions
%error('not finished need to remove the need for 2 fft for the complex and time it')
C=real(IMG);
if ~plan{4}
   C=single(C);
end
C=real(1/sqrt(1/(16*size(IMG,1)*size(IMG,2))))*C;
if size(IMG,3)>1
    C(:,1,:)=sqrt(2)*C(:,1,:);
    C(1,:,:)=sqrt(1/2)*C(1,:,:);
else
    C(:,1)=sqrt(2)*C(:,1);
    C(1,:)=sqrt(1/2)*C(1,:);
end
%C = real(sqrt(1/(16*size(IMG,1)*size(IMG,2)))*(MULT1.*(MULT2.*trans+conj(MULT2).*circshift(flip(trans,2),[0,1]))));
%C=(1/real(sqrt(1/(16*size(IMG,1)*size(IMG,2)))))*C;
D=flip(circshift(C,[0,-1]),2);
if size(IMG,3)>1
    D(:,1,:)=0;
else
    D(:,1)=0;
end
C= C - 1i*D ;%(real(C)-circshift(flip(imag(C),2),[0,1]));
C =( (1./plan{1}).*(1./plan{2}).*(C));
C=real(ifft2(C));
if size(IMG,3)>1
    C=reshape(C,size(IMG,1)*size(IMG,2),size(IMG,3));
    C(:,:)=C(plan{3},:);
    C=reshape(C,size(IMG,1),size(IMG,2),size(IMG,3));
else
    C(:)=C(plan{3});
end
end