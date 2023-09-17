function [labels] = generateSuperpixel(img,K)
[m,n,d]=size(img);
Y=reshape(img,m*n,d);
p = 1;%
[Y_pca] = pca(Y, p);
%
% Y_pca=reshape(Y_pca',m,n,p);
% img=mat2gray(Y_pca);
% img=imfilter(img,fspecial('average',3),'replicate','same');
% img1=im2uint8(img);
img1 = im2uint8(mat2gray(reshape(Y_pca',m,n, p)));
labels = mex_ers(double(img1),K);
labels=labels+1;
disp(size(labels));
label=reshape(labels,m,n);
end

