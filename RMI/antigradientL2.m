function y = antigradientL2(x,f)
%ANTIGRADIENTL2 이 함수의 요약 설명 위치
%   자세한 설명 위치
if nargin == 1
     f = fdct2_prep(size(x(:,:,1)),x(1) );
end
x(end,:,1)=0;
x(:,end,2)=0;
y=invPfunc(x(:,:,1)+x(:,:,2)-circshift(x(:,:,1), [1,0])-circshift(x(:,:,2), [0,1]),f);
end

