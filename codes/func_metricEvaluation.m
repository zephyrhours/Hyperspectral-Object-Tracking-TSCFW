function  Results = func_metricEvaluation(ground_truth,varargin)
%RESULTS = FUNC_METRICEVALUATION(ground_truth,varargin)
% Precision evaluation metrics for tracking traget, including:
%   1. MDP -> average(mean) distance precision (20 pixels as default)
%   2. CEL -> average center location error
%   3. MOP -> average overlap precision (0.5 as default)
%   4. AUC for precision plot
%   5. AUC for success plot
%   6. Precision plots of OPE
%   7. Success plots of OPE
%% Function Usage
% Results = func_metricEvaluation(ground_truth,pd_boxes,varargin)
% INPUTS:              
%  ground_truth -> ground truth, which must be Nx4-matrices where N is
%                 the number of sequence image in the tracking.The format
%                 is [X,Y,W,H], where X,Y is the left-top coordinates, W,H is 
%                 the width and height of tracking region, respectively.
%     varargin -> the bounding boxes predicted by different methods, which 
%                 aims to evaluate the accuracy of multiple methods simultaneously. 
%                 The format is:
%                 Results = func_evaluation(ground_truth,pd_boxes1,pd_boxes2,
%                 ...,pd_boxesn)
%                 where, pd_boxes is the tracker positions, which must be 
%                 Nx4-matrices where N is the number of sequence image in 
%                 the tracking. The format is [X,Y,W,H], where X,Y is the
%                 left-top coordinates, W,H is the width and height of
%                 tracking region, respectively.
%                 varargin{end}: the method name, where the format is 
%                 {'method1','method2',...,'methodn'}
% TIPS: 
% 1. The order of x and y does not matter here. Therefore, when calculating 
%    these evaluation indicators, there is no position exchange of x and y.
%
% Author: Zengfu Hou 
% Time: 2022-03-24
% Email: zephyrhou@126.com
% Website: https://zephyrhours.github.io/

% Copyright (c) 2022, Zengfu Hou
% All rights reserved.

%% ====================== Main function ===========================
%% Parameters Setting
thr_dp = 20;
MOP_threshold = 0.5;
%% Nargin Preference
nums_methods = length(varargin);
if nums_methods == 0
    error('There is no predicted location input here!');
    return
end

Results={};
ori_ground_truth = ground_truth;

if length(varargin)>=2
    if iscell(varargin{end-1}) && iscell(varargin{end})
        nums_methods = length(varargin)-2;
        names=varargin{end-1};
        location=varargin{end};
    elseif  ~iscell(varargin{end-1}) && iscell(varargin{end})
        nums_methods = length(varargin)-1;
        names=varargin{end};
        location = 'northeast';
    else
        nums_methods = length(varargin);
        names=varargin{end};
        location = 'northeast';
    end 
else
     names = false;
     location = false;
end


for i = 1:nums_methods
    positions = varargin{i};
    ground_truth = ori_ground_truth; %
    
    if size(positions,1) ~= size(ground_truth,1)
        n = min(size(positions,1), size(ground_truth,1));
        positions(n+1:end,:) = [];
        ground_truth(n+1:end,:) = [];
        disp(['Number of ground truth frames does not match number of tracked frames,','method-','num2str(i)']);
    end
    
    %% 1. Distance Precision (DP)
    % calculate the center location of bounding box
    ground_truth = [ground_truth(:,1:2) + (ground_truth(:,3:4) - 1) / 2 , ground_truth(:,3:4)]; 
    positions =[positions(:,1:2) + (positions(:,3:4) - 1) / 2 , positions(:,3:4)]; 

    distances = sqrt((positions(:,1) - ground_truth(:,1)).^2 + ...
        (positions(:,2) - ground_truth(:,2)).^2);
    distances(isnan(distances)) = [];
    % calculate distance precision
    distance_precision = nnz(distances < thr_dp) / numel(distances);

    %% 2. average Center Location Error (CLE)
    %calculate average center location error (CLE)
    average_CLE = mean(distances); 

    %% 3. Mean Overlap Precision (MOP)
    % calculate the overlap in each dimension
    overlap_height = min(positions(:,1) + positions(:,3)/2, ground_truth(:,1) + ground_truth(:,3)/2) ...
        - max(positions(:,1) - positions(:,3)/2, ground_truth(:,1) - ground_truth(:,3)/2);
    overlap_width = min(positions(:,2) + positions(:,4)/2, ground_truth(:,2) + ground_truth(:,4)/2) ...
        - max(positions(:,2) - positions(:,4)/2, ground_truth(:,2) - ground_truth(:,4)/2);

    % if no overlap, set to zero
    overlap_height(overlap_height < 0) = 0;
    overlap_width(overlap_width < 0) = 0;

    % remove NaN values (should not exist any)
    valid_ind = ~isnan(overlap_height) & ~isnan(overlap_width);

    % calculate area
    overlap_area = overlap_height(valid_ind) .* overlap_width(valid_ind);
    tracked_area = positions(valid_ind,3) .* positions(valid_ind,4);
    ground_truth_area = ground_truth(valid_ind,3) .* ground_truth(valid_ind,4);

    % calculate overlaps rate
    overlaps = overlap_area./ (tracked_area + ground_truth_area - overlap_area);

    % calculate Mean Overlap Precision (MOP)[success rate]
    MOP_precision = nnz(overlaps >= MOP_threshold) / numel(overlaps);
   
   %% ============================================================
   %% 1. Precision Plot of OPE (One-Pass Evaluation)
    % calculate and show precision plot
    max_threshold = 50;  % used for graphs in the paper
    precisions_rate = zeros(1,max_threshold);

    % calculate distances to ground truth over all frames
    distances = sqrt((positions(:,1) - ground_truth(:,1)).^2 + ...
                     (positions(:,2) - ground_truth(:,2)).^2);
    distances(isnan(distances)) = [];

    % compute precisions
    for p = 1:max_threshold
        precisions_rate(p) = nnz(distances < p) / numel(distances);
    end

    % calculate AUC value
%     LocErr=linspace(0.02,1,max_threshold);
%     LocErr=0.02:0.02:1;
%     AUC_precision=sum((LocErr(2:end)-LocErr(1:end-1)).*(precisions(2:end)+precisions(1:end-1))/2);
    AUC_precision = mean(precisions_rate);

   %% 2. Success Plot of OPE (One-Pass Evaluation)
    % calculate the overlap in each dimension
    overlap_width = min(positions(:,1) + positions(:,3)/2, ground_truth(:,1) + ground_truth(:,3)/2) ...
        - max(positions(:,1) - positions(:,3)/2, ground_truth(:,1) - ground_truth(:,3)/2);
    overlap_height = min(positions(:,2) + positions(:,4)/2, ground_truth(:,2) + ground_truth(:,4)/2) ...
        - max(positions(:,2) - positions(:,4)/2, ground_truth(:,2) - ground_truth(:,4)/2);

    % if no overlap, set to zero
    overlap_width(overlap_width < 0) = 0;
    overlap_height(overlap_height < 0) = 0;

    % remove NaN values (should not exist any)
    valid_ind = ~isnan(overlap_width) & ~isnan(overlap_height);

    % calculate area
    overlap_area = overlap_width(valid_ind) .* overlap_height(valid_ind);
    tracked_area = positions(valid_ind,3) .* positions(valid_ind,4);
    ground_truth_area = ground_truth(valid_ind,3) .* ground_truth(valid_ind,4);

    % calculate overlaps (intersection over union, IOU)
    overlaps = overlap_area ./ (tracked_area + ground_truth_area - overlap_area);

    % generate overlap rate thresholds
%     thr_overlap = linspace(0.02,1,max_threshold);
    thr_overlap=0.02:0.02:1; 
    success_rate=zeros(1,length(thr_overlap));

    for kk = 1:length(thr_overlap)
        % calculate success rate
        success_rate(kk) = nnz(overlaps >= thr_overlap(kk)) / numel(overlaps);
    end

    % calculate AUC value
%     AUC_success=sum((thr_overlap(2:end)-thr_overlap(1:end-1)).*(success_rate(2:end)+success_rate(1:end-1))/2);
    AUC_success = mean(success_rate);
       
   %% Display and Save the Results
   if iscell(names)
       disp(['-----------------Method:',names{i},'----------------------'])
   else
       disp(['------------------Method-',num2str(i),'----------------------'])
   end
    disp(['average distance precision(MDP@20p):',num2str(distance_precision)])
    disp(['average center location error(CEL):',num2str(average_CLE)])
    disp(['average overlap precision(MOP):',num2str(MOP_precision)])
    disp(['AUC value of precision plot:',num2str(AUC_precision)])
    disp(['AUC value of success plot:',num2str(AUC_success)])  
    disp('----------------------------------------------------------------')
    if iscell(names)
        results{1,1} = 'Method'; results{1,2} = names{i};
    else
        results{1,1} = 'Method'; results{1,2} = num2str(i);
    end
    results{2,1} = 'average distance precision(MDP@20p)'; results{2,2} = distance_precision;
    results{3,1} = 'average center location error(CEL)'; results{3,2} = average_CLE;
    results{4,1} = 'average overlap precision(MOP)'; results{4,2} = MOP_precision;
    results{5,1} = 'AUC value of precision plot'; results{5,2} = AUC_precision;
    results{6,1} = 'AUC value of success plot'; results{6,2} = AUC_success;
    results{7,1} = 'precisions'; results{7,2} = precisions_rate;
    results{8,1} = 'success rate'; results{8,2} = success_rate;
    Results{i} = results;
end

%% =============== Draw the Precisions/Success Plots===============
% 1. the precisions plots of one pass evaluation
figure
for i = 1:nums_methods
    results = Results{i};
    precisions_rate = results{7,2};
    AUC_precision = results{5,2};
    AUC_precision = roundn(AUC_precision,-3);
    plot(precisions_rate,'LineWidth',2);hold on;  
    if iscell(names)
        str{i}=[names{i},' [',num2str(AUC_precision),']'];
    else
        str{i}=['method',num2str(i),' [',num2str(AUC_precision),']'];
    end
end
xlabel('Location error threshold');
ylabel('Precision')
title('Precision plots of OPE')
if iscell(location)
    if size(location,2) >= 1
        legend(str,'Location',location{1,1})
    else
        legend(str) 
    end
else
    legend(str)
end
 
% 2.the success plots of one pass evaluation
figure;
for i = 1:nums_methods
    results = Results{i};
    success_rate = results{8,2};
    AUC_success = results{6,2};
    AUC_success = roundn(AUC_success,-3);
    plot(success_rate, 'LineWidth',2); hold on;
    if iscell(names)
        str2{i}=[names{i},' [',num2str(AUC_success),']'];
    else
        str2{i}=['method',num2str(i),' [',num2str(AUC_success),']'];
    end
end

XLabs=0:0.1:1; % Limited number of labels displayed 
for i =1:size(XLabs,2)
    XLable{i}=num2str(XLabs(i));
end
set(gca,'XTickLabel',XLable) 
xlabel('Overlap threshold');
ylabel('Success rate');
title('Success plots of OPE');

if iscell(location)
    if size(location,2) == 2
        legend(str2,'Location',location{1,2})
    elseif size(location,2) == 1
        legend(str2,'Location',location{1,1})
    else
        legend(str2)
    end
else
     legend(str2)
end


end
