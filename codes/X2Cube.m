function DataCube=X2Cube(I)
%Input Ximea 1088*2048 image
%Output 256*512*14 datacube
block=4;
I=double(I);
bandNumber=16;
% black=double(imread('black.bmp'));
% I=I(:,:,1)-black(:,:,1);
[x,y]=size(I);
DataCube = im2colstep(I',[block block],[block,block]);
sz=[y/block, x/block  bandNumber];
DataCube=reshape(DataCube',sz);
DataCube=permute(DataCube,[2,1,3]);

end
