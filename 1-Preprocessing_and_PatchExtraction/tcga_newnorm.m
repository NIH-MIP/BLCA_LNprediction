 
function [] = tcga_newnorm(refimg,svs_file)
    svsname = strsplit(svs_file,'.');
    try
        if(~exist([svsname{1} '_norm.mat']))
            caseData = bfGetReader(svs_file);
            caseMeta = caseData.getMetadataStore();
            caseMag = round(caseMeta.getObjectiveNominalMagnification(0,0).doubleValue(),2);

            %we are expecting all to be scanned the same, at 40x
            % selecting other levels based on what is available from all
            if(caseMag == 40)
                refratio = 4;
            elseif(caseMag == 20)
                refratio = 1;
            end

            lvls = zeros(caseMeta.getImageCount,3);
            for k = 1:caseMeta.getImageCount
                %i=1 is the reference level
                lvls(k,1) = eval(caseMeta.getPixelsSizeX(k-1)); 
                lvls(k,2) = eval(caseMeta.getPixelsSizeY(k-1)); 
                lvls(k,3) = lvls(1,1)/lvls(k,1); %mag difference from 40x
            end

            %find mag levels within file structure
            lvls_mag = lvls(:,3)';
            lvls_mag = round(lvls_mag(mod(log(round(lvls_mag))/log(2),1)==0));
            start_lvl = find(round(lvls(:,3)) == refratio);
            caseData.setSeries(start_lvl-1);
            I1 = bfGetPlane(caseData,1);
            I2 = bfGetPlane(caseData,2);
            I3 = bfGetPlane(caseData,3);
            mimg(:,:,1) = I1;
            mimg(:,:,2) = I2;
            mimg(:,:,3) = I3;
            [ NormMM ] = Norm(mimg, refimg, 'Macenko', 255, 0.15, 1, []);


            save([svsname{1} '.mat'],'mimg','-v7.3');
            save([svsname{1} '_norm.mat'],'NormMM','-v7.3');
        else
            disp('   ...already completed')
        end
    catch ME
        disp('   ...errored')
    end
end