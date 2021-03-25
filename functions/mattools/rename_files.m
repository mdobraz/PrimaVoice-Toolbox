

clear all
close all

a ='/hpc/banco/Primavoice_PNH/sounds/Stim_DB_2018/training_set/rename/NON_VOCAL/natural/living/';  % path
cd(a)
A =dir( fullfile(a, '*.*') ); % files to rename
  A(1:2,:)=[]

for i=1:numel(A)  % tous les noms
    CurrentName = A(i).name;  % noms maintenant
   NewName = sprintf('3_NV1_%d.wav', i); % prefix a ajouter
    movefile(CurrentName,NewName) % changer nom en newname
end