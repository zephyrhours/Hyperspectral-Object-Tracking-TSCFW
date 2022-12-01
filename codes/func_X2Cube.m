function hsi = func_X2Cube(Img,win)
%FUNC_X2Cube Convert the XIMEA or XIMEC mosaic Image to a cube dataset
%% Function Usage
% hsi = func_X2Cube(Img,win)
%
% INPUTS:
%       Img -> the XIMEA or XIMEC mosaic Image;
%       win -> the mosaic size(By default, XIMEA is 4, XIMEC is 5);
%
% OUTPUS:
%       hsi -> the hyperspectral dataset(cube dataset);
%
% Author: Zephyr Hou
% Time: 2021-09-07

% All rights reserved.
% Email: zephyrhours@gmail.com

%% Main Function
% By default, the dataset is XIMEC dataset with the mosaic size of 5.
if nargin < 2 || isempty(win)
    win = 5;  
    warning('Missing input parameter. By default, the dataset is XIMEC dataset with the mosaic size of 5!')
end

%% Method 1:(slow)
% image = double(Img);
% [rows,cols] = size(image);
% hsi = zeros(floor(rows/win),floor(cols/win),floor(win*win)); 
% for ii = 1:rows/win
%     for jj = 1:cols/win
%         patchi = image((ii-1)*win+1:ii*win,(jj-1)*win+1:jj*win);
%         hsi(ii,jj,:) = reshape(patchi',[],1);
%     end
% end

%% Method 2:(fast)
rows=floor(size(Img,1)/win);
cols=floor(size(Img,2)/win);
bands=win*win;
cutImg=double(Img(1:rows*win,1:cols*win));
patchImg=im2col(cutImg',[win,win],'distinct'); % bands x (rows x cols)
hsi=reshape(patchImg',cols,rows,bands);    % cols x rows x bands
hsi=permute(hsi,[2,1,3]);                  % rows x cols x bands

end

