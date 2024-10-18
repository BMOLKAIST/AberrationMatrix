function [result, phi_in, phi_out] = function_class_multi(distortion_ur_refocused_t122,pos_x,pos_y)
%%
R = permute(distortion_ur_refocused_t122,[2 3 1]);
R = reshape(R,1024,1024,15,15);
for ii = 1:15
    for jj = 1:15
        R(:,:,ii,jj) = circshift(fft2(R(:,:,ii,jj)),[-40*(jj-15), -40*(ii-15)]);
    end
end


%%
mask = sqrt(((1:1024)-260-5).^2  +  ((1:1024).'-260-5).^2) < 260;
illumination_idx = [21	22	23	24	25	34	35	36	37	38	39	40	41	42	48	49	50	51	52	53	54	55	56	57	58	63	64	65	66	67	68	69	70	71	72	73	74	77	78	79	80	81	82	83	84	85	86	87	88	89	92	93	94	95	96	97	98	99	100	101	102	103	104	105	107	108	109	110	111	112	113	114	115	116	117	118	119	120	122	123	124	125	126	127	128	129	130	131	132	133	134	135	137	138	139	140	141	142	143	144	145	146	147	148	149	150	153	154	155	156	157	158	159	160	161	162	163	164	168	169	170	171	172	173	174	175	176	177	178	179	184	185	186	187	188	189	190	191	192	193	200	201	202	203	204	205	206	207	217	218	219	220];
mask_in=zeros(15,15);
mask_in(illumination_idx)=1;

wx = 1024*pos_x;
wy = 1024*pos_y;

window_out = zeros(1024,1024);
x = -60:60;
y = x.';
f = exp(-(x.^2+y.^2)/(2*20));
window_out(1:121,1:121) = f;
window_out = abs(ifft2(window_out));
window_out = fftshift(window_out/window_out(1));
window_out = circshift(window_out,round([wy,wx])-ceil((1024+1)/2));
%
%imagesc(window_out)

R_windowed = permute(distortion_ur_refocused_t122,[2 3 1]);
R_windowed = reshape(R_windowed,1024,1024,15,15);
for ii = 1:15
    for jj = 1:15
        R_windowed(:,:,ii,jj) = circshift(fft2(R_windowed(:,:,ii,jj).*window_out),[-40*(jj-15), -40*(ii-15)]);
    end
end
%%
[pout_tot, pin_tot] = get_aberration_class(R_windowed);
pout_tot = gpuArray(single(pout_tot));
pin_tot = gpuArray(single(pin_tot));
R = gpuArray(single(R));
mask = gpuArray(single(mask));
mask_in = gpuArray(single(mask_in));

R_corrected1 = zeros(1024,1024,15,15);
for ii = 1:15
    for jj = 1:15
        R_corrected1(:,:,ii,jj) = circshift((mask.*exp(-1i*pout_tot).*R(:,:,ii,jj)),[40*(jj-15), 40*(ii-15)]);
        R_corrected1(:,:,ii,jj) = mask_in(jj,ii).*exp(-1i*pin_tot(1,1,ii,jj))* R_corrected1(:,:,ii,jj);

    end
end


R_corrected1 = reshape(R_corrected1,1024,1024,[]);
R_corrected1 = ifft2(R_corrected1);

result = gather(abs(sum(R_corrected1,3)));
phi_in = gather(pin_tot);
phi_out = gather(pout_tot);
end