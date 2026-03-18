FID = fopen('IDs.txt', 'r');
C = textscan(FID, '%s');
fclose(FID);
prefixes = C{1};

for prefixIDX = 1:length(prefixes)
    currentprefix = prefixes{prefixIDX}; 
    cd ([currentprefix]);
    addpath /home/clairea/wa26_scratch/Claire_Heart_Rate_Per_Slice/
    siemens_vd_ppu3t_matlab_script;
    cd ..;
end
