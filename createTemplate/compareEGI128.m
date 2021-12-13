

fiffGeo = fiff_read_evoked('newGEO.fif',1);
fiffOld = fiff_read_evoked('old.fif',1);
fiffNew = fiff_read_evoked('new.fif',1);

for iElec = 1:fiffOld.info.nchan
    oldElec(iElec,1:3) = fiffOld.info.chs(iElec).eeg_loc(1:3,1)';
    newElec(iElec,1:3) = fiffNew.info.chs(iElec).eeg_loc(1:3,1)';
    geoElec(iElec,1:3) = fiffGeo.info.chs(iElec).eeg_loc(1:3,1)';
end
    
figure;hold on;
scatter3(oldElec(:,1),oldElec(:,2),oldElec(:,3))
scatter3(geoElec(:,1),geoElec(:,2),geoElec(:,3),'filled')
scatter3(newElec(:,1),newElec(:,2),newElec(:,3),'filled')