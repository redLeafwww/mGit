clear;
close all;
load mprage_8ch_k-space_577203641;
fimage = ifft(fftshift(kspace, 1), [], 1);
fimage = ifft(fftshift(fimage, 2), [], 2);
ndim = size(fimage)
montage(reshape(abs(fimage), [ndim(1), ndim(2), 1, ndim(3)]));
sos = sqrt(sum(abs(fimage).*abs(fimage), 3));
figure; imshow(sos, []);

%=========Sensitivity map calculation with low resolution images==========
nlow = 20;
lkspace = zeros(ndim);
lkspace((ndim(1)-nlow)/2+1:(ndim(1)+nlow)/2, (ndim(2)-nlow)/2+1:(ndim(2)+nlow)/2, :) = kspace((ndim(1)-nlow)/2+1:(ndim(1)+nlow)/2, (ndim(2)-nlow)/2+1:(ndim(2)+nlow)/2, :);
limage = ifft(fftshift(lkspace, 1), [], 1);
limage = ifft(fftshift(limage, 2), [], 2);
figure; montage(reshape(abs(limage), [ndim(1), ndim(2), 1, ndim(3)]));
soslow = sqrt(sum(abs(limage).*abs(limage), 3));
figure; imshow(soslow, []);
%mask = soslow>(max(max(soslow))*0.05);
%figure; imshow(mask);
soslow = reshape(soslow, [ndim(1)*ndim(2), 1])*ones(1,ndim(3));
%masklow = reshape(mask, [ndim(1)*ndim(2), 1])*ones(1,ndim(3));
sensemap = reshape(limage, [ndim(1)*ndim(2), ndim(3)])./soslow;%.*masklow;
sensemap = reshape(sensemap, [ndim(1), ndim(2), ndim(3)]);
figure; montage(reshape(abs(sensemap), [ndim(1), ndim(2), 1, ndim(3)]));

%=========Undersampling data==========
R = 4;
pkspace = kspace(1:R:ndim(1), :, :);
pimage = ifft(fftshift(pkspace, 1), [], 1);
pimage = ifft(fftshift(pimage, 2), [], 2);
npdim = size(pimage)
figure;montage(reshape(abs(pimage), [npdim(1), npdim(2), 1, npdim(3)]));

%=========Unfold data==========
ufimage = zeros(ndim(1), ndim(2));
gfactor = zeros(ndim(1), ndim(2));
for pe = 1:npdim(1),
   for fe = 1:ndim(2),
       AT = squeeze(sensemap(pe:ndim(1)/R:ndim(1), fe, :));
       A = transpose(AT);
       y = squeeze(pimage(pe, fe, :));
       x = A\y;
       %x = pinv(A)*y;
       ufimage(pe:ndim(1)/R:ndim(1), fe) = x;
       CC = AT*A;
       invCC = inv(CC);
       gfactor(pe:ndim(1)/R:ndim(1), fe) = [diag(invCC()).*diag(CC)];
   end
end
figure; imshow(abs(ufimage),[]);
imcontrast;
figure; imshow(abs(gfactor),[]);
imcontrast;
