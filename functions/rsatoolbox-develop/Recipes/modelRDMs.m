%  modelRDMs is a user-editable function which specifies the models which
%  brain-region RDMs should be compared to, and which specifies which kinds of
%  analysis should be performed.
%
%  Models should be stored in the "Models" struct as a single field labeled
%  with the model's name (use underscores in stead of spaces).
%  
%  Cai Wingfield 11-2009

function Models = modelRDMs()

nconditions=16;
% define index vectors
Human = 1:4;
Macaque = 5:8;
Marmoset = 9:12;
NonVocal = 13:16;

voice = [Human Macaque Marmoset];
nonhuman = [Macaque Marmoset NonVocal];
nonmacaque = [Human Marmoset NonVocal];

Models.human = ones(nconditions,nconditions);
Models.human(Human,Human)=0;
Models.human(nonhuman,nonhuman)=0;
Models.human(logical(eye(nconditions)))=0; % fix the zero-diagonal

Models.macaque = ones(nconditions,nconditions);
Models.macaque(Macaque,Macaque)=0;
Models.macaque(nonmacaque,nonmacaque)=0;
Models.macaque(logical(eye(nconditions)))=0; % fix the zero-diagonal


Models.voice = ones(nconditions,nconditions);
Models.voice(voice,voice)=0;
Models.voice(NonVocal,NonVocal)=0;
Models.voice(logical(eye(nconditions)))=0; % fix the zero-diagonal

Models.random = squareform(pdist(rand(nconditions,nconditions)));

