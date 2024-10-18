function result = uncorrected_image(distortion_ur_refocused_t122)

R = permute(distortion_ur_refocused_t122,[2 3 1]);
R = reshape(R,1024,1024,15,15);
for ii = 1:15
    for jj = 1:15
        R(:,:,ii,jj) = fft2(R(:,:,ii,jj));
    end
end


%%
mask = sqrt(((1:1024)-260-5).^2  +  ((1:1024).'-260-5).^2) < 260;
illumination_idx = [21	22	23	24	25	34	35	36	37	38	39	40	41	42	48	49	50	51	52	53	54	55	56	57	58	63	64	65	66	67	68	69	70	71	72	73	74	77	78	79	80	81	82	83	84	85	86	87	88	89	92	93	94	95	96	97	98	99	100	101	102	103	104	105	107	108	109	110	111	112	113	114	115	116	117	118	119	120	122	123	124	125	126	127	128	129	130	131	132	133	134	135	137	138	139	140	141	142	143	144	145	146	147	148	149	150	153	154	155	156	157	158	159	160	161	162	163	164	168	169	170	171	172	173	174	175	176	177	178	179	184	185	186	187	188	189	190	191	192	193	200	201	202	203	204	205	206	207	217	218	219	220];
mask_in=zeros(15,15);
mask_in(illumination_idx)=1;

R = gpuArray(single(R));
mask = gpuArray(single(mask));
mask_in = gpuArray(single(mask_in));


for ii = 1:15
    for jj = 1:15
        R(:,:,ii,jj) = mask_in(jj,ii).*mask.*R(:,:,ii,jj);
    end
end

R = reshape(R,1024,1024,[]);
R = ifft2(R);

result = gather(abs(sum(R,3)));

end

