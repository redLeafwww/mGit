ufimage = zeros(256, 256);
gfactor = zeros(256, 256);
for pe = 1:128,
   for fe = 1:256,
       AT = squeeze(map(pe:256/R:256, fe, :));
       A = transpose(AT);
       y = squeeze(imaliased(pe, fe, :));
       x = A\y;
       %x = pinv(A)*y;
       ufimage(pe:256/R:256, fe) = x;
       CC = AT*A;
       invCC = inv(CC);
       gfactor(pe:256/R:256, fe) = [diag(invCC()).*diag(CC)];
   end
end
figure; imshow(abs(ufimage),[]);
imcontrast;
figure; imshow(abs(gfactor),[]);
imcontrast;