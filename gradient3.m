function g = gradient3(x)
%GRADIENTS 이 함수의 요약 설명 위치
%   자세한 설명 위치

g=circshift(x,[-1, 0, 0])-x;
g(:,:,:,2)=circshift(x,[0, -1, 0])-x;
g(:,:,:,3)=circshift(x,[0, 0, -1])-x;
% g(end,:,:,1)=g(end-1,:,:,1);
% g(:,end,:,1)=g(:,end-1,:,1);
% g(:,:,end,1)=g(:,:,end-1,1);
g(end,:,:,1)=0;
g(:,end,:,2)=0;
g(:,:,end,3)=0;
end

