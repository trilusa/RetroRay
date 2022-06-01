%% Takes all the checkpoint data files and combines them into a 
%{9x42}nested cell arrary. Each cell has an Nx6 array
%   [[srcXYZ destXYZ]_1; ... ;[srcXYZ destXYZ]_N ]
%where N is total number of rays detected in each trial.

if ~exist('data','var')
        error('The workspace does not contain any data. Load data with ./load_data.m');
end

%%
reshaped_data=cell(9,41);

for s = 1:9 %there are 9 cases
    for x=1:41
       for t = 1:length(data{s}) %t trial, this varies depending on how much data there is
           if(isfield(data{s}{t},'detected_this_trial'))
               reshaped_data{s,x} = [reshaped_data{s,x}; data{s}{t}.detected_this_trial{x}];
           else
               reshaped_data{s,x} = [reshaped_data{s,x}; data{s}{t}.detected{x}];
           end
       end
    end
end

%%


%%
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

%%
test=cell(9,1)
for i=1:9
    test{i}=cell(randi(9,1,1)+1,1);
end
test
for i=1:9
    for t=1:length(test{i})
    test{i}{t}=rand(randi(9,1,1)-1,6);
    end
end
