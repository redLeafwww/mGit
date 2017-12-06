a_kernel=-2:2;
b_kernel=-R-1:R:R-1;
kernel=zeros(a,b,c,c);
kspace_gr=kspace_us;
for coil=1:c%对于每一个线圈的图
    for shape=1:R-1
        A=[];
        B=[];
        for i=1-a_kernel(1):a-a_kernel(end)
            for j=ACSline(1):ACSline(2)
                if mod(j-1,R)==shape
                    A0=reshape(kspace_us(i+a_kernel,j-shape+1+b_kernel,:),1,[]);
                    A=[A;A0];
                    B=[B;kspace_us(i,j,coil)];
                end
            end
        end
        x=pinv(A)*B;
        for i=1-a_kernel(1):a-a_kernel(end)
            for j=1:b
                if (((j+b_kernel(1)>0&&j<ACSline(1))||(j>ACSline(2)&&j+b_kernel(end)<b))&&(mod(j-1,R)==shape))
                    A0=reshape(kspace_gr(i+a_kernel,j-shape+1+b_kernel,:),1,[]);
                    B0=A0*x;
                    kspace_gr(i,j,coil)=B0;
                end
            end
        end
    end
end
im_gr=fftshift(ifft2(fftshift(kspace_gr,2)),2);%kspace
recon_gr=sqrt(sum(im_gr.*im_gr,3));

imax_gr=max(max(max(abs(im_gr))));%最大的那个数
figure('Name', 'GRAPPA images from phase array channels');
montage(reshape(abs(im_gr),256,256,1,8),[0 imax_gr]);%under-sampled后的图像
recon_gr=sqrt(sum(abs(im_gr).*abs(im_gr),3));%SOS
figure('Name','GRAPPA Multiple coil combination with SOS');
imshow(mat2gray(abs(recon_gr))),title 'SOS';
                
                
                