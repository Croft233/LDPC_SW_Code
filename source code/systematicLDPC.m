z = 28
H = LDPC.construct802p16e(672,336);
[m,n] = size(H);
H = H(1:end-28,:);

%this is the decoder version
Hdec = H;
figure(1)
spy(Hdec)
for jj = 0:9
    for ii = 1:z
        H(ii+z+z*jj,:) = mod(H(ii+z*jj+z,:) + H(ii+z*jj,:),2);
    end
end
Henc = H;
figure(2)
spy(Henc)

save LDPCn672k364 Hdec Henc