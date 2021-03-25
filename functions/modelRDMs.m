%  modelRDMs is a user-editable function which specifies the models which
%  brain-region RDMs should be compared to, and which specifies which kinds of
%  analysis should be performed.
%
%  Models should be stored in the "Models" struct as a single field labeled
%  with the model's name (use underscores in stead of spaces).
%  

% function Models = modelRDMs()

clear Models

nconditions=16;
% define index vectors
Human = 1:4;
Macaque = 5:8;
Marmoset = 9:12;
NonVocal = 13:16;

voice = [Human Macaque Marmoset];
nonhuman = [Macaque Marmoset NonVocal];
nonmacaque = [Human Marmoset NonVocal];
nonmarmoset = [Human Macaque NonVocal];

%% Human
Models(1).RDM = ones(nconditions,nconditions);
Models(1).RDM(Human,Human)=0;
Models(1).RDM(nonhuman,nonhuman)=0;
Models(1).RDM(logical(eye(nconditions)))=0; % fix the zero-diagonal
Models(1).name = 'Human';
Models(1).color = rand(1,3);

%% Macaque
Models(2).RDM = ones(nconditions,nconditions);
Models(2).RDM(Macaque,Macaque)=0;
Models(2).RDM(nonmacaque,nonmacaque)=0;
Models(2).RDM(logical(eye(nconditions)))=0; % fix the zero-diagonal
Models(2).name = 'Macaque';
Models(2).color = rand(1,3);

% %% Marmoset
% Models(3).RDM = ones(nconditions,nconditions);
% Models(3).RDM(Marmoset,Marmoset)=0;
% Models(3).RDM(nonmarmoset,nonmarmoset)=0;
% Models(3).RDM(logical(eye(nconditions)))=0; % fix the zero-diagonal
% Models(3).name = 'Marmoset';
% Models(3).color = rand(1,3);

%% Voice
Models(3).RDM = ones(nconditions,nconditions);
Models(3).RDM(voice,voice)=0;
Models(3).RDM(NonVocal,NonVocal)=0;
Models(3).RDM(logical(eye(nconditions)))=0; % fix the zero-diagonal
Models(3).name = 'NonVocal';
Models(3).color = rand(1,3);

%% Human vs Macaque
Models(4).RDM = zeros(nconditions,nconditions);
% Models(4).RDM(Human,Human)=0;
% Models(4).RDM(Macaque,Macaque)=0;
Models(4).RDM(Human,Macaque)=1;
Models(4).RDM(Macaque,Human)=1;
% Models(4).RDM(logical(eye(nconditions)))=0; % fix the zero-diagonal
Models(4).name = 'HvM';
Models(4).color = rand(1,3);

% %% Random
% Models(5).RDM = squareform(pdist(rand(nconditions,nconditions)));
% Models(5).name = 'Random';
% Models(5).color = rand(1,3);


