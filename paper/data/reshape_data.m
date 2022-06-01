%%
%remove the h
reshaped_data=cell(1,9);

clc


if ~exist('data','var')
        error('The workspace does not contain any data. Load data with ./load_data.m');
end

for i = 1:length(data)
    tic
    for j=length(data{i})
        if(isfield(data{i}{j},'detected_this_trial'))
            data{i}{j}.detected=data{i}{j}.detected_this_trial;
        end
        for k = 1:length(data{i}{j}.detected)
            if(~isempty(data{i}{j}.detected{k}))
                srcXY=data{i}{j}.detected{k}(:,1:2);
                destXY=data{i}{j}.detected{k}(:,4:5);
                reshaped_data{i}=[reshaped_data{i} ; srcXY destXY-1];
            end
        end
    end
    toc
end