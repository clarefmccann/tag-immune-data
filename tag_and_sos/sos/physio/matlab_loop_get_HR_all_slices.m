FID = fopen('IDs.txt', 'r');
C = textscan(FID, '%s');
fclose(FID);
prefixes = C{1};

%The concatenating/interleaving part of the loop below is from Lars: https://github.com/translationalneuromodeling/tapas/issues/213

hr3 = zeros(17424,length(prefixes));

for prefixIDX = 1:length(prefixes)
    currentprefix = prefixes{prefixIDX};
    for iSlice = 1:72
    load(fullfile(currentprefix, 'physio_out', sprintf('physio_slice%03d.mat', iSlice)));
    hr1(:,iSlice) = physio.ons_secs.hr;
    hr2 = reshape(hr1', [],1); % the ' is important to create the transpose and mix row/column vectors correctly
    end
    hr3(:,prefixIDX) = hr2
end

T = table(hr3)
writetable(T, 'Heart_Rate_Data_Task_Sequence_All_Slices.csv', 'WriteVariableNames', false)

