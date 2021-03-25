
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        %% Filename
        bold_file = SR(s).filename{r};
        [~,bold_name,ext] = fileparts(bold_file);

        out_file = fullfile(paths.resliced,[reslice_flags.prefix bold_name ext]);
       
        %% Put to double datatype
        [~,datatype_out] = system([paths.FSL_prefix 'fslval ' out_file ' datatype']);
        datatype_out = str2double(datatype_out);
        
        if datatype_out ~= 64
            fprintf('Converting session %i run %i to ''double''\n',session,run)
            system(sprintf('gzip %s',out_file));
            system(sprintf('%sfslmaths %s %s -odt double',paths.FSL_prefix,out_file,out_file)); % apply values to resampled template
            system(sprintf('gunzip %s',out_file));
        end
      
    end
end