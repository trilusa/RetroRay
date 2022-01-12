
% D_d_sim = 0.004:0.002:0.016;     % detector diameter (1mm to simulate point source)
D_d_sim = (.1:.05:4)/1000;
v_sim = 0:0.2:.8;    % horizontal displacements to simulate
h=1.5;

suffix = 1;%:25;
N=1e3;
data = zeros(length(D_d_sim),length(v_sim));
for s = suffix
data_all_temp=[];
for d=1:length(D_d_sim)
    data_temp = [];
    fn =['monte-carlo-alpha_N' num2str(N) '_PD' num2str(D_d_sim(d)*1000) 'mm_L' num2str(50) 'mm(' num2str(suffix(s)) ').mat'];
    load(fn, 'detected_this_trial');
    for i=1:length(detected_this_trial)
        data_temp = [data_temp length(detected_this_trial{i})];
    end
    data_all_temp = [data_all_temp; data_temp];
end
data = data + data_all_temp;
end
data=data/6;
thetas = atand(v_sim/h);
data_cos  = data .* cosd(thetas);
data_cos_d = data_cos ./ (2* sqrt(v_sim.^2 + h^2));
data_cos_d_norm = data_cos_d ./ data_cos_d(:,1);

