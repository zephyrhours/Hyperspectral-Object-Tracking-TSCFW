function [ feature_image ] = get_hog3d( im, fparam, gparam )
%extract fhog features using piotrs toolbox. Currently takes no parameters
%except hog-cell-size
% if ~isfield(fparam, 'nOrients')
%     fparam.nOrients = 9;
% end
[im_height, im_width, num_im_chan, num_images] = size(im);
feature_image = zeros(floor(im_height/gparam.cell_size), floor(im_width/gparam.cell_size), fparam.nDim, num_images, 'single');
for k = 1:num_images
    hog_image = hog3d(single(im(:,:,:,k)), gparam.cell_size, fparam.hOrients,fparam.vOrients);
    %the last dimension is all 0 so we can discard it
    feature_image(:,:,:,k) = hog_image;
end
end