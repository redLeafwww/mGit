load mprage_8ch_k-space_577203641.mat;
[a,b,c]=size(kspace);
%% 
%展示原始的k-space和原始的图像
figure('Name','Original k-space');
montage(reshape(log(abs(kspace)),a,b,1,c),'DisplayRange',[]);%原始k-space
im=ifft2(fftshift(kspace,2));%原始图像
imax1=max(max(max(abs(im))));%最大的那个数
figure('Name', 'Images from phase array channels');
montage(reshape(abs(im),256,256,1,8),[0 imax1]);%原始图像
recon1=sqrt(sum(abs(im).*abs(im),3));%SOS
figure('Name','Multiple coil combination with SOS');
imshow(mat2gray(abs(recon1))),title 'SOS';
%% under-sampled
ACS=32;
kernel=4;
R=4;
kspace_us=zeros(a,b,c);
kspace_us(1:R:a,:,:)=kspace(1:R:a,:,:);
kspace_us((a-ACS)/2+1:a,:,:)=kspace((a-ACS)/2+1:a,:,:);
%figure('Name','under sampled k-space');
%montage(reshape(log(abs(kspace_us)),a,b,1,c),'DisplayRange',[]);%原始k-space
im_us=ifft2(fftshift(kspace_us,2));%kspace
imax_us=max(max(max(abs(im_us))));%最大的那个数
figure('Name', 'under sampled images from phase array channels');
montage(reshape(abs(im_us),256,256,1,8),[0 imax_us]);%under-sampled后的图像
recon_us=sqrt(sum(abs(im_us).*abs(im_us),3));%SOS
figure('Name','under sampled Multiple coil combination with SOS');
imshow(mat2gray(abs(recon_us))),title 'SOS';
%% GRAPPA
kspace_gr=kspace_us;
%每幅图4个点解出1个未知点，8幅图共32个点解1个点
%解出n
n_pre=zeros(ACS-kernel,32,c);%每一个slice上可以解出28个n，每个n有32个数字
n=zeros(32,c);%n_pre最小二乘得到
S=zeros(32,b);%32个点，共kx个点
for k=1:8%8幅图，每幅图解n_pre
    for m=(a-ACS)/2+1+2:(a+ACS)/2-2%28条已知线，解出28个n
        y=squeeze(kspace_gr(m,:,k));
        for ll=1:256
            for kk=1:8%8幅图
                %for nn=1:4%每个图4个点
                    S(4*kk-3,ll)=kspace_gr(m-2,ll,kk);
                    S(4*kk-2,ll)=kspace_gr(m-1,ll,kk);
                    S(4*kk-1,ll)=kspace_gr(m+1,ll,kk);
                    S(4*kk,ll)=kspace_gr(m+2,ll,kk);
                %end
            end
        end
        n_pre(m-114,:,k)=transpose(S)\transpose(y);
    end
end
n=squeeze(n_pre(14,:,:));
for k=1:8
    for m=4:R:113
        for ll=1:256
            for kk=1:8%8幅图
                %for nn=1:4%每个图4个点
                    S(4*kk-3,ll)=kspace_gr(m-2,ll,kk);
                    S(4*kk-2,ll)=kspace_gr(m-1,ll,kk);
                    S(4*kk-1,ll)=kspace_gr(m+1,ll,kk);
                    S(4*kk,ll)=kspace_gr(m+2,ll,kk);
                %end
            end
        end
        kspace_gr(m,:,k)=transpose(n(:,k))*S;
    end
end
im_gr=ifft2(fftshift(kspace_gr,2));%kspace
imax_gr=max(max(max(abs(im_gr))));%最大的那个数
figure('Name', 'GRAPPA images from phase array channels');
montage(reshape(abs(im_gr),256,256,1,8),[0 imax_gr]);%under-sampled后的图像
recon_gr=sqrt(sum(abs(im_gr).*abs(im_gr),3));%SOS
figure('Name','GRAPPA Multiple coil combination with SOS');
imshow(mat2gray(abs(recon_gr))),title 'SOS';
RMSE=norm(abs(recon_gr)-abs(recon1),'fro')/norm(abs(recon1),'fro');
%% g factor
MM=100
for i=1:MM
    for coil=1:c
        noise(:,:,coil)=wgn(a,b,20);
    end
    noise1=noise*0;noise1(1:2:end,:,:)=noise(1:2:end,:,:);
    knoise=GRAPPA(noise1);%此处是将前文GRAPPA段代码变成函数
    for coil=1:c
        noiseim(:,:,i)=ifft2(fftshift(knoise(:,:,coil),2));
    end
    d=(res)-sos(noiseim);
end
gmap=d;
figure,imshow(matgray(abs(gmap)));
        