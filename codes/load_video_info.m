function [seq, ground_truth] = load_video_info(video_path)
% Author: ZephyrHou
% Time: 2022-01-05

ground_truth = dlmread([video_path '/HSI/groundtruth_rect.txt']);
% ground_truth=ground_truth(11:end,:);

% 利用第一帧数据的ground truth进行初始化
seq.len = size(ground_truth, 1);
seq.init_rect = ground_truth(1,:);

img_path = [video_path '/HSI/'];
img_files = dir(fullfile(img_path, '*.png'));
img_files = {img_files.name};

% img_files = img_files(4:end);
% img_files = [img_path img_files];
% if exist([img_path num2str(1, '%04i.png')], 'file'),
%     img_files = num2str((1:seq.len)', [img_path '%04i.png']);
% elseif exist([img_path num2str(1, '%04i.jpg')], 'file'),
%     img_files = num2str((1:seq.len)', [img_path '%04i.jpg']);
% elseif exist([img_path num2str(1, '%04i.bmp')], 'file'),
%     img_files = num2str((1:seq.len)', [img_path '%04i.bmp']);
% else
%     error('No image files to load.')
% end
% img_files=img_files(11:end);

seq.s_frames = cellstr(img_files);

% load RGB image for display
if exist([video_path '/HSI-FalseColor/'],'dir')
    RGBimg_files = dir(fullfile([video_path '/HSI-FalseColor/'], '*.jpg'));
    RGBimg_files = {RGBimg_files.name};
    seq.s_frames_rgb = cellstr(RGBimg_files);
elseif exist([video_path '/RGB/'],'dir')
    RGBimg_files = dir(fullfile([video_path '/RGB/'], '*.jpg'));
    RGBimg_files = {RGBimg_files.name};
    seq.s_frames_rgb = cellstr(RGBimg_files);
else
   seq.s_frames_rgb = false;
end

end
