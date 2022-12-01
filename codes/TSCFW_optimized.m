function [results] = TSCFW_optimized(params)
%% setting parameters
video_path     = params.video_path;
visualization  = params.visualization;
target_sz      = floor(params.target_sz);
s_frames       = params.s_frames;
s_frames_rgb   = params.s_frames_rgb;
output_sigma_factor = params.output_sigma_factor;  %spatial bandwidth (proportional to target)
feature_type   = params.feature_type;
num_frames     = params.num_frames;
pos            = floor(params.init_pos);  % center point position

% model parameters
sigma_d        = params.sigma_d;
sigma_e        = params.sigma_e;
lambda         = params.lambda;  
lambda_sr      = params.lambda_sr;
lambd_ca       = params.lambd_ca;

% scale parameters
nScales        = params.nScales;
scale_step     = params.scale_step;

% video name 
slash=strfind(video_path,'\');
video_name = video_path(slash(end)+1:end);
save_path = params.save_path;

% others
filter_max_area = params.filter_max_area;
init_target_sz      = target_sz;
featureRatio = params.t_global.cell_size;   

% target_area determined
if prod(target_sz) < 200   
    featureRatio = 1;
end
    
search_area_scale   = params.search_area_scale;
% when the number of cells are small, choose a smaller cell size
search_area = prod(init_target_sz / featureRatio * search_area_scale);

if search_area < params.t_global.cell_selection_thresh * filter_max_area
    featureRatio = params.t_global.cell_size;
    search_area = prod(init_target_sz / featureRatio * search_area_scale);
end

if search_area > filter_max_area
    currentScaleFactor = repmat(sqrt(search_area / filter_max_area),1,2);
else
    currentScaleFactor = [1.0, 1.0];
end

% target size at the initial scale
base_target_sz = target_sz./currentScaleFactor;
% window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        sz = floor( base_target_sz * search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        sz = floor(repmat(sqrt(prod(base_target_sz * search_area_scale)), 1, 2)); % square area, ignores the target aspect ratio
    case 'fix_padding'
        sz = floor(base_target_sz + sqrt(prod(base_target_sz * search_area_scale) + (base_target_sz(1) - base_target_sz(2))/4) - sum(base_target_sz)/2); % const padding
    otherwise
        error('Unknown "params.search_area_shape". Must be "proportional", "square" or "fix_padding"');
end

% set the size to exactly match the cell size
sz = round(sz / featureRatio) * featureRatio;
use_sz = floor(sz/featureRatio);

% construct the label function-correlation output, 2D gaussian function,
% with a peak located upon the target
output_sigma = sqrt(prod(floor(base_target_sz/featureRatio))) * output_sigma_factor;
rg           = circshift(-floor((use_sz(1)-1)/2):ceil((use_sz(1)-1)/2), [0 -floor((use_sz(1)-1)/2)]);
cg           = circshift(-floor((use_sz(2)-1)/2):ceil((use_sz(2)-1)/2), [0 -floor((use_sz(2)-1)/2)]);
[rs, cs]     = ndgrid( rg,cg);
y            = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf           = fft2(y);   %   FFT of y.

cos_window = single(hann(use_sz(1))*hann(use_sz(2))');

im = imread([video_path '/HSI/' s_frames{1}]);
im = X2Cube(im);

%% 2D space scale factors
if nScales > 0
    scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2));
    scaleFactors = scale_step .^ scale_exp;
    min_scale_factor = scale_step ^ ceil(log(max(5./ sz)) / log(scale_step));
    max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
end

total_nScales = nScales*nScales;
[scaleFactorsX,scaleFactorsY] = meshgrid(scaleFactors,scaleFactors');
scaleFactorsX = reshape(scaleFactorsX,1,total_nScales);   
scaleFactorsY = reshape(scaleFactorsY,1,total_nScales);   

%% newton interpolate
% Pre-computes the grid that is used for socre optimization
ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';
newton_iterations = params.newton_iterations;

%% feature selection
features.feature_type   = params.feature_type;
switch feature_type
	case 'gray'
		learning_rate = params.learning_rate;  %linear interpolation factor for adaptation
		cell_size = featureRatio;
	case 'shog'
		learning_rate = params.learning_rate; 
		features.hog_orientations = 9;
		cell_size = featureRatio;
    case 'sshmg'
		learning_rate = params.learning_rate; 
		features.hog_3d_hOrients = 9;
        features.hog_3d_vOrients = featureRatio;
		cell_size = featureRatio;      
end

features.window_size = sz;

%% Initialization
% initialize the projection matrix (x,y,h,w)
rect_position  = zeros(num_frames, 4);  
response_template = [];
time = 0;
loop_frame = 1;

%% TSCFW tracker
for frame = 1:numel(s_frames)
    % load image
    try
        im = imread([video_path '/HSI/' s_frames{frame}]);
    catch
        try
            im = imread([s_frames{frame}]);
        catch
            im = imread([video_path '/' s_frames{frame}]);
        end
    end
    im = X2Cube(im);
    im = func_Normalized(im,2);
    
    if iscell(s_frames_rgb)
        try
            try
                RGBim = imread([video_path '/HSI-FalseColor/' s_frames_rgb{frame}]);
            catch
                RGBim = imread([video_path '/RGB/' s_frames_rgb{frame}]);
            end
            if size(RGBim,1) ~= size(im,1) || size(RGBim,2) ~= size(im,2)
                s_frames_rgb = 0;   
            end
        catch
            s_frames_rgb = false;   
        end
    end

    tic();
    if frame > 1
         for i = 1:total_nScales
            tmp_sz = round(sz.*currentScaleFactor.*[scaleFactorsX(i),scaleFactorsY(i)]);
            patch_scale = get_subwindow(im, pos, tmp_sz);
            patch = imresize(patch_scale,sz);           
            zf = fft2(get_features_hsi(patch, features, cell_size, cos_window)); 
            response_template(:,:,i) = sum(model_wf .* zf,3);  %equation for fast detection
         end
        responsef_padded = resizeDFT2(response_template, use_sz);
        
        % response in the spatial domain
        response = ifft2(responsef_padded, 'symmetric');
        [disp_row, disp_col, szid,~,~] = resp_newton(response, responsef_padded, newton_iterations, ky, kx, use_sz);
        
        translation_vec = round([disp_row, disp_col]*cell_size.*currentScaleFactor.*[scaleFactorsX(szid),scaleFactorsY(szid)]);
        pos = pos + translation_vec;

        currentScaleFactor = currentScaleFactor.*[scaleFactorsX(szid),scaleFactorsY(szid)];
        
       if currentScaleFactor(1) < min_scale_factor
            currentScaleFactor(1) = min_scale_factor;
       elseif currentScaleFactor(1) > max_scale_factor
            currentScaleFactor(1) = max_scale_factor;
       end
       
      if currentScaleFactor(2) < min_scale_factor
            currentScaleFactor(2) = min_scale_factor;
       elseif currentScaleFactor(2) > max_scale_factor
            currentScaleFactor(2) = max_scale_factor;
      end    
    end
    
    % update template by using tensor processing
    target_sz = round(base_target_sz.*currentScaleFactor);

    % window size, taking padding into account
    switch params.search_area_shape
        case 'proportional'
            tmp_sz = floor(target_sz * search_area_scale);     % proportional area, same aspect ratio as the target
        case 'square'
            tmp_sz = floor(repmat(sqrt(prod(target_sz * search_area_scale)), 1, 2)); % square area, ignores the target aspect ratio
        case 'fix_padding'
            tmp_sz = floor(target_sz + sqrt(prod(target_sz * search_area_scale) + (target_sz(1) - target_sz(2))/4) - sum(target_sz)/2); % const padding
        otherwise
            error('Unknown "params.search_area_shape". Must be ''proportional'', ''square'' or ''fix_padding''');
    end
 
    patch_scale = get_subwindow(im, pos, tmp_sz);
    patch = imresize(patch_scale,sz); 
    xf = fft2(get_tensor_features_hsi(double(patch), features, cell_size, cos_window));  

    %% Content-Aware
     offset = [-target_sz(1) 0; 0 -target_sz(2); target_sz(1) 0; 0 target_sz(2)];
    
    kfn = zeros([size(xf) length(offset)]);
    for j=1:length(offset)
        %obtain a subwindow close to target for regression to 0
        patchn = get_subwindow(im, pos+offset(j,:),target_sz);
        patchn = imresize(patchn,sz); 
        xfn = fft2(get_features_hsi(patchn, features, cell_size, cos_window));
        kfn(:,:,:,j) = conj(xfn) .*xfn;
    end
     
    kf = conj(xf) .* xf; 
    sz_xf = size(xf);
    patch_ssw = imresize(patch, [sz_xf(1) sz_xf(2)]); 
    phi = func_SSW(patch_ssw,sigma_d,sigma_e);
    fz = conj(xf).*yf;
    fm = (1+ lambda_sr)*kf+lambda*(conj(phi) .* phi)+lambd_ca.*sum(kfn,4);
    wf = fz./fm;
    
    %first frame, train with a single image
    if frame == 1  
        model_wf = wf;
    else
        model_wf = (1 - learning_rate) * model_wf + learning_rate * wf;
    end
    
    box = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])]; % (x,y,w,h)
    
    % save position and time
	time = time + toc();
    rect_position(frame,:) = box;

	%visualization
    if visualization == 1
        rect_position_vis = box;
        im_to_show = im;
        if size(im_to_show,3) == 1
            im_to_show = repmat(im_to_show, [1 1 3]);
        end
        if frame == 1
            fig_handle = figure('Name', ['Tracking_TSCFW','     VideoName: ',video_name ]);
            if iscell(s_frames_rgb)
                imshow(RGBim,[]);
            else
                hyperImshow(im_to_show);               
            end
            hold on;
            axis off;axis image;
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(20, 30, ['# Frame : ' int2str(1) ' / ' int2str(num_frames)], 'color', [1 0 0], 'BackgroundColor', [1 1 1], 'fontsize', 16);
            hold off;
         
            set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
        else
            figure(fig_handle);
            if iscell(s_frames_rgb)
                imshow(RGBim,[]);
            else
                 hyperImshow(im_to_show);
            end
            hold on;
                      
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(20, 30, ['# Frame : ' int2str(loop_frame) ' / ' int2str(num_frames)], 'color', [1 0 0], 'BackgroundColor', [1 1 1], 'fontsize', 16);
            text(20, 60, ['FPS : ' num2str(loop_frame/time)], 'color', [1 0 0], 'BackgroundColor', [1 1 1], 'fontsize', 16);  
            hold off;
        end
        drawnow
    end
    loop_frame = loop_frame + 1; 
end

fps = loop_frame / time;
results.type = 'rect';
results.res  = rect_position;
results.fps  = fps;
results.time = time;
end

