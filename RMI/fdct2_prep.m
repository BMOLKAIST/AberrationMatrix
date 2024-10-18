function plan = fdct2_prep(sizing,classRepresent)
%%% Herve code, slightly modified by KRLee 20220308
    
sizing=squeeze([sizing(1),sizing(2)]);

vec1 = cast(1:sizing(1),'like',classRepresent);
vec2 = cast(1:sizing(2),'like',classRepresent);
[k2,k1]=meshgrid(vec2, vec1);

POS1=zeros(sizing,'like',classRepresent);
POS1(1:2:sizing(1),:,:)=k1(1:length(1:2:sizing(1)),:,:);
POS1(2:2:sizing(1),:,:)=flip(k1(end-length(2:2:sizing(1))+1:end,:,:),1);

%POS1=k1;

POS2=zeros(sizing,'like',classRepresent);
POS2(:,1:2:sizing(2),:)=k2(:,1:length(1:2:sizing(2)),:);
POS2(:,2:2:sizing(2),:)=flip(k2(:,end-length(2:2:sizing(2))+1:end,:),2);
%POS2=k2;

POS=round(sub2ind(sizing, POS1,POS2));
POS = cast(POS,'like',classRepresent);

%MULT=single(2*exp(-1i.* pi.*((k1-1)./(2*sizing(1))+(k2-1)./(2*sizing(2)))));
MULT1=2*exp(-1i.* pi.*((k1(:,1)-1)./(2*sizing(1))));
MULT2=2*exp(-1i.* pi.*((k2(1,:)-1)./(2*sizing(2))));
%MULT=MULT1.*MULT2;
plan=squeeze(cell(3,1));

plan{1}=MULT1;
plan{2}=MULT2;
plan{3}=POS;
plan{4}=false;%indicate that computation ar to be done in double
end