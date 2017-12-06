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

%% sensitivity map

%低分辨率
kspacenew=zeros(a,b,c);
for i=1:8
    for j=1:20
        for k=1:20
            %kspacenew
            kspacenew(a/2-10+j,b/2-10+k,i)=kspace(a/2-10+j,b/2-10+k,i);
        end
    end
end
imnew=ifft2(fftshift(kspacenew,2));%低分辨率的图像
imax2=max(max(max(abs(imnew))));%最大的那个数
%figure('Name', 'Low resolution images from phase array channels');
%montage(reshape(abs(imnew),256,256,1,8),[0 imax2]);
recon2=sqrt(sum(abs(imnew).*abs(imnew),3));%低分辨率的SOS
%figure('Name','Low resolution multiple coil combination with SOS');
%imshow(mat2gray(abs(recon2))),title 'SOS';
%求sensitivity map
for i=1:8
    map(:,:,i)=imnew(:,:,i)./recon2;
end
maxmap=max(max(max(abs(map))));
figure('Name', 'Maps from phase array channels');
montage(reshape(abs(map),256,256,1,8),[0 maxmap]);

%
%% Generate aliased images
R=2;
subsamplekspace=zeros(floor(a/R),b,c);
for i=1:floor(a/R)
    subsamplekspace(i,:,:)=kspace(R*i-1,:,:);
end
for i=1:c
    imaliased(:,:,i)=ifft2(fftshift(subsamplekspace(:,:,i)));
end
recon3=sqrt(sum(abs(imaliased).*abs(imaliased),3));
figure('Name','Aliased images');
imshow(mat2gray(abs(recon3))),title 'Aliased images SOS';

%Sense reconstruction
imsense=zeros(a,b);
gfactor=zeros(a,b);
%求M，利用公式，S、Ma(r)也就是imaliased
%S=zeros(8,1);
gfactor=zeros(256,256);
%for i=1:a/R
%    for j=1:b
%        S=squeeze(map(i:a/R:a,j,:));
%        SH=transpose(S);
%        y=squeeze(imaliased(i,j,:));
        %if (i>128)
        %    m=i-128;
        %else
        %    m=i;
        %end
        %M(i,j)=inv((S'*S))*(S')*squeeze(imaliased(m,j,:));%解码后的图像
        %M(i,j)=inv(transpose(S)*S)*transpose(S)*squeeze(imaliased(m,j,:));
%        M(i:a/R:a,j)=SH\y;
%        gfactor(i:256/R:256,j)=diag(inv(SH*S)).*diag(SH*S);%sqrt((S'*S)\(S'*S));%g factor
%    end
%end
for pe = 1:b/R,
   for fe = 1:a,
       AT = squeeze(map(pe:256/R:256, fe, :));
       A = transpose(AT);
       y = squeeze(imaliased(pe, fe, :));
       x = A\y;
       %x = pinv(A)*y;
       imsense(pe:256/R:256, fe) = x;
       CC = AT*A;
       invCC = inv(CC);
       gfactor(pe:256/R:256, fe) = [diag(invCC()).*diag(CC)];
   end
end
figure; imshow(abs(imsense),[]);
imcontrast;
figure; imshow(abs(gfactor),[]);
imcontrast;
dif=abs(imsense)./max(max(abs(imsense)))-abs(recon1)./max(max(abs(recon1)));
figure;imshow(abs(dif),[]);
rmse=sqrt(sum(sum(dif.^2))./(a*b));
%recon4=sqrt(sum(abs(M).*abs(M),3));%SENSE后的SOS
%figure('Name','SENSE SOS');
%imshow(mat2gray(abs(recon4))),title 'SENSE SOS';
%for i=1:floor(a/R)
%    for j=1:b
%        S=[];
%        for k=1:R
%            S=[S,squeeze(map(i+(k-1)*floor(a/R),j,:))];
%        end
%        X=(S'*S)\(S');
%        M=X*
%    end
%end

