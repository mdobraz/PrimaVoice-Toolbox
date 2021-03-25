clear all
close all
root='/hpc/banco/Primavoice_PNH/sounds/Stim_DB_2018/training_set/ALL_cut/'; % to change
cd(root);


filesall = dir('*.wav');

 % select all audio files in this folder
 
 nfiles=length(filesall);
 
 for k = 1:nfiles
     info{k}=audioinfo(filesall(k).name);
     names{k,1}=info{1,k}.Filename;
     [filepath{k,1},name{k,1},ext{k,1}] = fileparts(names{k,1});
     dur{k,1}=info{1,k}.Duration;  % get the durations
     Sf(k,1)=info{1,k}.SampleRate; % get the sample rates/freq 
 end
 
result=[name,dur]
