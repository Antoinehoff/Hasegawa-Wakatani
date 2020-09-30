%% Save variables
if START 
    dataname = 'data_00.mat';
else
    dataname = ['data_',sprintf('%.2d',SID),'.mat'];
end
save([FIGDIR,dataname],'ZZ','NN','PP','Ts','Skip');
disp('Variables saved')