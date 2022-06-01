%% Takes all the checkpoint data files and combines them into a 
%{1x9} nested cell arrary. Each cell has a {1x42} cell array (one for each horis displacement
% Each endtry in that cell array containing a Nx6 array
%   [[srcXYZ destXYZ]_1; ... ;[srcXYZ destXYZ]_N ]
%where N is total number of rays detected in each trial.

if ~exist('data','var')
        error('The workspace does not contain any data. Load data with ./load_data.m');
end

reshaped_data=cell(1,9);
for i = 1:length(data)
    for j=length(data{i})
        if(isfield(data{i}{j},'detected_this_trial'))
            temp_cell=data{i}{j}.detected_this_trial;
        else
            temp_cell=data{i}{j}.detected;
        end
        
        for k = 1:length(temp_cell)
            if(~isempty(temp_cell{k}))
                srcXY=temp_cell{k}(:,1:2);
                destXY=temp_cell{k}(:,4:5);
                reshaped_data{i}=[reshaped_data{i} ; srcXY destXY-1];
            end
        end
    end
end