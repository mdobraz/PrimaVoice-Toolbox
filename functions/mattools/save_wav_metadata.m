function save_wav_metadata(stims_dir,txt_file)
%%%% Add metadata from textfile to the 'Comment' field

a = tdfread(txt_file);
fn = fieldnames(a);

for i = 1:length(a.(fn{1}))
    filename = sprintf('%s/%s.wav',stims_dir,strtrim(a.(fn{1})(i,:))); % strtrim removes leading and trailing whitespaces from string array
    [y,fs] = audioread(filename);
    
    categ_str = [];
    for j = 2:length(fn)
        categ_str = [categ_str strtrim(a.(fn{j})(i,:))];
        if j < length(fn)
            categ_str = [categ_str ','];
        end
    end
    audiowrite(filename,y,fs,'Comment',categ_str)
end

