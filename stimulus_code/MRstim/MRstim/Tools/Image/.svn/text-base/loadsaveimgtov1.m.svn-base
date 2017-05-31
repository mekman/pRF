imgbase = 'newton';
img = double(imread([imgbase '.tif']));

img = flipdim(img,2);

bgindex = [0 0 255];
for n=1:3
    [imout(:,:,n), orgim(:,:,n), orgimtot(:,:,n)] = fromimagetov1(img(:,:,n),bgindex(n));
end

imwrite(uint8(imout),[imgbase 'v1.tif'],'tif');
imwrite(flipdim(uint8(orgim),2),[imgbase 'org.tif'],'tif');
imwrite(flipdim(uint8(orgimtot),2),[imgbase 'orgtot.tif'],'tif');
