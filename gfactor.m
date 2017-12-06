M=200;
noise=zeros(size(kspace));
for i=1:M
    noise = (2*rand(size(kspace))-1)+noise;
end
subkspace = subkspace+noise + noise*sqrt(-1);
for i=1:floor(a/R)
    B2=[];
    if(i==1)
        B2=[B2;zeros(c,a)];
        for j=1:kernel-1
            B2=[B2;transpose(squeeze(subkspace(1+(i-1)*r+(j-1)*r,:,:)))];
        end
    elseif(i>=2 && i<=floor(b/r)-2)
        for j=1:kernel
            B2=[B2;transpose(squeeze(subkspace(1+(i-1)*r+(j-2)*r,:,:)))];
        end
    elseif(i==floor(b/r)-1)
        for j=1:kernel-1
            B2=[B2;transpose(squeeze(subkspace(1+(i-1)*r+(j-2)*r,:,:)))];
        end
        B2=[B2;zeros(c,a)];
    elseif(i==floor(b/r))
        for j=1:kernel-2
            B2=[B2;transpose(squeeze(subkspace(1+(i-1)*r+(j-2)*r,:,:)))];
        end
        B2=[B2;zeros(2*c,a)];
    end
    A2=co*B2;
    for j=1:r-1
        for s = 1:c
            subkspace(2+(j-1)+(i-1)*r,:,s)= A2(s+(j-1)*c,:);
        end
    end
end

for i=1:8
    subkspace(start:finish,:,i) = kspace(start:finish,:,i);
end
noisegrappaim=zeros(size(kspace));
for i = 1:c
    noisegrappaim(:,:,i)=ifft2(subkspace(:,:,i));
end
finalnoiseim = sos(noisegrappaim);
finalnoiseim = finalnoiseim./M;
g_map = finalgrappaim./finalnoiseim;
