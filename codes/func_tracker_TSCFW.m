function [pd_boxes, time] = func_tracker_TSCFW(video_path, img_files, ground_truth)

%% Parameters Setting
feature_type   = 'sshmg';  % 'gray','shog','sshmg','sshmgshog';
learning_rate  = 0.0025;
nScales        = 7;
scale_step     = 1.01;   %1.005
lambda         = 1e1;  %[1e-6,...,1e5];  
lambda_sr      = 1;
lambd_ca       = 0.1;

% standard deviation of the desired correlation output (proportional to target)
output_sigma_factor = 0.1;	        % [1/16,0.1]	

% result process
visualization  = 0;   % 0 or 1
save_gif = 0;         % 0 or 1

%% Display Parameters
disp('--------------------------- TSCFW Mode ----------------------------')
disp('Parameters Setting£º')
disp(['feature_type: ',feature_type])
disp(['learning_rate: ',num2str(learning_rate)])
disp(['nScales: ',num2str(nScales)])
disp(['scale_step: ',num2str(scale_step)])
disp(['lambda: ',num2str(lambda)])
disp(['lambda_sr: ',num2str(lambda_sr)])
disp(['lambd_ca: ',num2str(lambd_ca)])
disp(['output_sigma_factor: ',num2str(output_sigma_factor)])
disp('-------------------------------------------------------------------')
%% Initialization
base_path = video_path;
slash=strfind(base_path,'\');
data_name = base_path(slash(end-1)+1:slash(end)-1);

seq.len = size(ground_truth, 1);
seq.init_rect = ground_truth(1,:);
seq.s_frames = img_files;
seq.VidName = data_name;
st_frame = 1;
en_frame = seq.len;
seq.st_frame = st_frame;
seq.en_frame = en_frame;
seq.s_frames_rgb = rgb_img_files;

seq.visualization = visualization; 
seq.save_gif = save_gif;

seq.nScales = nScales;
seq.scale_step = scale_step;  
seq.lambda = lambda;    
seq.lambda_sr = lambda_sr;
seq.lambd_ca = lambd_ca;
seq.feature_type = feature_type;
seq.output_sigma_factor = output_sigma_factor;
    
% Run - main function
results = run_TSCFW(seq, video_path, learning_rate);
pd_boxes = results.res;
time =  results.time;
    
FPS = results.fps;
disp(['FPS: ',num2str(FPS)])

end

