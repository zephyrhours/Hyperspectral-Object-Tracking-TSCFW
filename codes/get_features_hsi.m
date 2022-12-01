function x = get_features_hsi(hsi, features, cell_size, cos_window)
%GET_TENSOR_FEATURES
%   Extracts tensor dense features from image.
%
%   X = GET_TENSOR_FEATURES(IM, FEATURES, CELL_SIZE)
%   Extracts features specified in struct FEATURES, from image IM. The
%   features should be densely sampled, in cells or intervals of CELL_SIZE.
%   The output has size [height in cells, width in cells, features].
%
%   To specify HOG features, set field 'hog' to true, and
%   'hog_orientations' to the number of bins.
%
%   To experiment with other features simply add them to this function
%   and include any needed parameters in the FEATURES struct. To allow
%   combinations of features, stack them with x = cat(3, x, new_feat).

%%  Feacture Extraction
nbins = 8; % bins of intensity historgram
nwindow = floor(features.window_size);
feature_type = features.feature_type;

bands = size(hsi,3);
im_features=[];

switch feature_type
	case 'gray'
        % gray features
        for i = 1:bands
             im = hsi(:,:,i);
             x = im - mean(im(:));
             x = cell_grayscale(im, cell_size);
            im_features=cat(3,im_features,x);
        end
        x = im_features;
        
	case 'shog'
        % stacked hog features
        for i = 1:bands
             im = hsi(:,:,i);
             x = double(fhog(single(im), cell_size, features.hog_orientations));
             x(:,:,end) = [];  % remove all-zeros channel ("truncation feature")
             im_features = cat(3,im_features,x);
        end
        % intensity information
        im = sum(hsi,3);
        x = double(fhog(single(im), cell_size, features.hog_orientations));
        x(:,:,end) = [];  %remove all-zeros channel ("truncation feature")
        x = cat(3,im_features,x);
    case 'sshmg'
        % sshmg features
        x = hog3d(single(hsi), cell_size,features.hog_3d_hOrients,features.hog_3d_vOrients);
end

%% process with cosine window if needed
if ~isempty(cos_window)
    x = bsxfun(@times, x, cos_window);
end
    
end
