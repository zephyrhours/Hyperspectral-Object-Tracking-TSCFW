% Compiled by Zengfu Hou 
function results = run_TSCFW(seq, video_path, lr)

params.target_sz = [seq.init_rect(1,4), seq.init_rect(1,3)];  
params.video_path = video_path;
params.learning_rate       = lr;                   

% select the feature type and expanding region 
target_area = prod(params.target_sz);
aspect_ratio = params.target_sz(1)/params.target_sz(2);

% padding selection
if aspect_ratio < 1
    params.padding = 2;                            
    params.search_area_shape   = 'proportional';     
elseif (aspect_ratio >= 1 && aspect_ratio < 1.4) && (target_area < 1000)
    params.padding = 2.5;                            
    params.search_area_shape   = 'proportional'; 
elseif (aspect_ratio >= 1 && aspect_ratio < 1.4) && (target_area >=1000)
    params.padding = 1.5;                             
    params.search_area_shape   = 'proportional';     
elseif (aspect_ratio >= 1.4 && aspect_ratio < 2)   
    params.padding = 1.5;                        
    params.search_area_shape   = 'square';       
elseif aspect_ratio >= 2.9  && (target_area > 800) 
    params.padding = 1;                     
    params.search_area_shape   = 'fix_padding';  
else
    params.padding = 1.5;                      
    params.search_area_shape   = 'proportional';         
end
    
params.feature_type        = seq.feature_type; 

params.search_area_scale   = 1+ params.padding;

%standard deviation of the desired correlation output (proportional to target)
params.output_sigma_factor = seq.output_sigma_factor;   % [1/16,0.1]	 
params.visualization       = seq.visualization;         % 1 or 0 

params.filter_max_area     = seq.filter_max_area;
params.newton_iterations   = 50;

% Threshold for reducing the cell size in low-resolution cases
params.t_global.cell_selection_thresh = 0.75^2; 

% feature cell size
switch params.feature_type
    case 'gray'
        params.t_global.cell_size  = 1;                  
    otherwise
        params.t_global.cell_size  = 4;                  
end

% size, position, frames initialization
params.init_pos       = [seq.init_rect(1,2), seq.init_rect(1,1)] + floor(params.target_sz/2);
params.s_frames       = seq.s_frames;
params.s_frames_rgb   = seq.s_frames_rgb;
params.num_frames     = seq.en_frame - seq.st_frame + 1;
params.seq_st_frame   = seq.st_frame;
params.seq_en_frame   = seq.en_frame;
params.save_path      = seq.save_path;

% model parameters
params.sigma_d = 1;
params.sigma_e = 1;
params.lambda         = seq.lambda;  %[1e-6,...,1e5];  
params.lambda_sr      = seq.lambda_sr;
params.lambd_ca       = seq.lambd_ca;

% scale parameters
params.nScales        = seq.nScales;
params.scale_step     = seq.scale_step;   

% run the main function
results = TSCFW_optimized(params);
end
