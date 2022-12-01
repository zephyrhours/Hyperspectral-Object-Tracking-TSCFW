%% TSCFW Tracker for Hyperspectral Object Tracking
clear;clc;close all;
disp('--------------TSCFW Tracker-------------')
disp('Author: Zengfu Hou')
disp(['Current Time: ',datestr(now)])
curr_time = num2str(fix(clock));
curr_time = strrep(curr_time,' ','');
result_name = [getenv('UserProfile'),'\Desktop\','Result_',curr_time,'.xlsx'];
%% Datasets Setting
base_path = 'E:\data\whisper\test';  % your hyperspectral videos path
save_path = strcat(cd,'\results');
if ~exist('results','dir')
   mkdir(save_path);
end
delete(strcat(save_path,'\*.mat'))
addpath(genpath(pwd));
%% Scenarios

videos={'ball';'basketball';'board';'book';'bus';'bus2';'campus';'car';'car2';'car3';'card';'coin';'coke';...
    'drive';'excavator';'face';'face2';'forest';'forest2';'fruit';'hand';'kangaroo';'paper';'pedestrain';...
    'pedestrian2';'player';'playground';'rider1';'rider2';'rubik';'student';'toy1';'toy2';'truck';'worker'};

version = 'v1';
%% Parameters Setting
visualization  = 1;         % 0 or 1
feature_type   = 'shog';  % default feature
learning_rate  = 0.0019;
nScales        = 3;
scale_step     = 1.05; 
lambda         = 1e1;  %[1e-6,...,1e5];  
lambda_sr      = 1;
lambd_ca       = 0.1;

% standard deviation of the desired correlation output (proportional to target)
output_sigma_factor = 0.1;	        % [1/16,0.1]	
filter_max_area   = 50^2;        % the size of the training/detection area in feature grid cells

%% Display Parameters
disp('--------------------------- TSCFW Mode ----------------------------')
disp('Parameters Setting£º')
disp(['deault feature_type: ',feature_type])
disp(['learning_rate: ',num2str(learning_rate)])
disp(['nScales: ',num2str(nScales)])
disp(['scale_step: ',num2str(scale_step)])
disp(['lambda: ',num2str(lambda)])
disp(['lambda_sr: ',num2str(lambda_sr)])
disp(['lambd_ca: ',num2str(lambd_ca)])
disp(['output_sigma_factor: ',num2str(output_sigma_factor)])
disp('-------------------------------------------------------------------')
%% Initialization
videoNum=numel(videos); 
total_frame =zeros(videoNum,1);
total_gt=[];
total_pd=[];
total_time =[];
total_nums =[];

%% Main 

tic;
for vid = 1:videoNum
    close all;
    disp(['The Scene is: ', videos{vid}])
    
    video_path = [base_path '\' videos{vid}];
    [seq, ground_truth] = load_video_info(video_path);
    seq.VidName = videos{vid};
    st_frame = 1;
    en_frame = seq.len;
    seq.st_frame = st_frame;
    seq.en_frame = en_frame;
   
    seq.visualization = visualization; 
    seq.save_path = save_path;
    
    seq.nScales = nScales;
    seq.scale_step = scale_step;  
    seq.lambda = lambda;    
    seq.lambda_sr = lambda_sr;
    seq.lambd_ca = lambd_ca;
    seq.feature_type = feature_type;
    seq.output_sigma_factor = output_sigma_factor;
    seq.filter_max_area  = filter_max_area;
      
    % Run - main function
    results  = run_TSCFW(seq, video_path, learning_rate);
    pd_boxes = results.res;
    FPS      = results.fps;
    time     = results.time;
    disp(['FPS: ',num2str(results.fps)])

    Results = func_metricEvaluation(ground_truth,pd_boxes,{'TSCFW'});
    
    % stack all to one for overall evaluation
    total_gt   = [total_gt; ground_truth];
    total_pd   = [total_pd; pd_boxes];
    total_nums = [total_nums; seq.len];
    total_time = [total_time; time];
    
    % move result to the specified folder
    results.gt = ground_truth;
    result=matfile(sprintf(strcat(videos{vid},'trackingTSCFW.mat')),'writable',true);
    result.results=results;
    movefile(strcat(videos{vid},'trackingTSCFW.mat'),save_path)  
  
end
overall_time = toc;

%% Evaluation
close all;

Results = func_metricEvaluation(total_gt,total_pd,{'TSCFW'});
disp(['average FPS: ',num2str(sum(total_nums)/sum(total_time))])

disp(['Finised Time: ',datestr(now)])



