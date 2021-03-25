function rel_cog = vol_rel_cog(in_file,paths)

[~,cogchar] = system(sprintf('%sfslstats %s -C',paths.FSL_prefix,in_file));
cog = str2double(strsplit(cogchar,' '));
cog = cog(1:3);
dims = zeros(1,3);
for i = 1:3
    [~,comres] = system(sprintf('%sfslval %s dim%i',paths.FSL_prefix,in_file,i));
    dims(i) = str2double(comres);
end
rel_cog = cog ./ dims;


