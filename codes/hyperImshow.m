function [rgb] = hyperImshow( hsi, RGBbands )
%% Hyperspectral Image color display
% Author: Zephyr Hou
% Time: 2019-12-02
% Function Usage
%   Default pseudo color display, you can also define RGB band to display the
%   true color
% Input:
%    hsi¡ªthe 3D hyperspectral dataset with the size of rows x cols x bands
%    RGBbands¡ª the RGB bands to be displayed, with the format [R G B]
% Output:
%    rgb- the finual RGB map with the RGB bands(the size is rows x cols x 3)  
%% Main Function

hsi=double(hsi);
[rows, cols, bands] = size(hsi);

minVal =min(hsi(:));
maxVal=max(hsi(:));
normalizedData=hsi-minVal;

if(maxVal==minVal)
    normalizedData=zeros(size(hsi));
else
    normalizedData=normalizedData./(maxVal-minVal);
end

hsi=normalizedData;

if (nargin == 1)
    RGBbands = [bands round(bands/2) 1];
end

if (bands ==1)
    red = hsi(:,:);
    green = hsi(:,:);
    blue = hsi(:,:);
else
    red = hsi(:,:,RGBbands(1));
    green = hsi(:,:,RGBbands(2));
    blue = hsi(:,:,RGBbands(3));
end
rgb = zeros(size(hsi, 1), size(hsi, 2), 3);
rgb(:,:,1) = adapthisteq(red);      % Adaptive histogram equalization
rgb(:,:,2) = adapthisteq(green);
rgb(:,:,3) = adapthisteq(blue);
imshow(rgb); axis image;
end
