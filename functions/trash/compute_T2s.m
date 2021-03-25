% function [T2s,T2s_map] = compute_T2s(ME_file,ME_file_info,Yroi)

n_echos = 3;
fields = num2str((1:n_echos)');
fields = cellstr([repmat('echo',3,1) fields]);

ME_file = cell(n_echos,1);
ME_file_info = cell(n_echos,1);
SI_ME = cell(n_echos,1);
x = nan(1,n_echos);

for echo = 1:n_echos
    ME_file{echo} = sprintf('/hpc/banco/Primavoice_Data_and_Analysis/sub-Maga/ses-02/anat/sub-Maga_ses-02_run-02_echo-%02.0f_epi.nii',echo);
    ME_file_info{echo} = sprintf('/hpc/banco/Primavoice_Data_and_Analysis/sub-Maga/ses-02/anat/sub-Maga_ses-02_run-02_echo-%02.0f_epi.json',echo);
    SI_ME{echo} = spm_jsonread(ME_file_info{echo});
    x(echo) = SI_ME{echo}.EchoTime * 1000;
end

for echo = 1:n_echos
    P = spm_vol(ME_file{echo});
    Ytmp = spm_read_vols(P);
    if echo == 1
        Echos = mean(Ytmp,4);
    else
        Echos(:,:,:,echo) = mean(Ytmp,4);
    end
end
clear Ytmp


if ~exist('Yroi','var')
    Yroi = ones(size(Echos(:,:,:,1)));
end

% monoexponential function
y = @(b,x) b(1).*exp(-b(2).*x);


opts = optimset('MaxFunEvals',100000,'MaxIter',100000);
T2s_map = nan(size(Yroi));
% S0map = nan(size(Yroi));

fprintf('Computing T2s_map...\n')
tic
for i = 1:size(Echos,1)
    for j = 1:size(Echos,2)
        for k = 1:size(Echos,3)
            if Yroi(i,j,k) % if voxel is in the ROI mask
                yx = squeeze(Echos(i,j,k,:))'; % all echos
                [~,I] = sort(yx,'descend');
                if all(I==1:n_echos)
                    OLS = @(b) sum((y(b,x)-yx).^2);
                    B = fminsearch(OLS,rand(2,1),opts);
                    T2s_map(i,j,k) = 1/B(2);
                else
                    T2s_map(i,j,k) = nan;
                end
                if T2s_map(i,j,k) > 100
                    T2s_map(i,j,k) = nan;
                end
%                 S0map(i,j,k)=B(1);
            end
        end
    end
end
toc

T2s.mean_T2 = nanmean(T2s_map(:)) / 1000; % mean T2 of the ROI from ms to sec
T2s.mean_R = 1 / T2s.mean_T2; % in s-1
T2s.scan_time = SI_ME{1}.AcquisitionTime;



T2map_file = '/hpc/banco/Primavoice_Data_and_Analysis/sub-Maga/ses-02/anat/T2map_3TE.nii';
PT2 = P(1);
PT2.fname = T2map_file;
spm_write_vol(PT2,T2s_map);


