function RGB = complex2rgb(C,Ampmax,colormap2D)
Cth = min(abs(C),0.99*Ampmax);
C = Cth.*exp(1i*angle(C));
delta=2*Ampmax/(colormap2D.N-1);
X=round((real(C)+Ampmax)/delta)+1;
Y=round((-imag(C)+Ampmax)/delta)+1;
X(X<1)=1; X(X>colormap2D.N)=colormap2D.N;
Y(Y<1)=1; Y(Y>colormap2D.N)=colormap2D.N;
Ind = Y + colormap2D.N * (X - 1);
RGB = colormap2D.data(Ind);
RGB(:,:,2) = colormap2D.data(Ind+colormap2D.N^2);
RGB(:,:,3) = colormap2D.data(Ind+2*colormap2D.N^2);
end

