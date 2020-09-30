%% Load data after setup.m or during run.m
dataname = ['../results/',BASIC.SIMID,'/data_%.2d.mat'];
dataname = sprintf(dataname,SID);
disp(['Loading ',dataname,'...'])
load(dataname,'ZZ','NN','PP','Ts','Skip'); Ns = numel(Ts);
disp('.. Done');
[~ , its0] = min(abs(t_rst-Ts));
ZZ_0 = ZZ(:,:,its0);
NN_0 = NN(:,:,its0);
PP_0 = PP(:,:,its0);
it0  = 2; t0 = Ts(its0);