% load svs file and xml file, pull at low res and visualize
function [] = load_disp_tcga(xml_file, svs_file, saveDir, tcga_id, savename)
    xy = xml_parse(xml_file);
    
    caseData = bfGetReader(svs_file);
    caseMeta = caseData.getMetadataStore();
    caseMag = round(caseMeta.getObjectiveNominalMagnification(0,0).doubleValue(),2)

    %we are expecting all to be scanned the same, at 40x
    % selecting other levels based on what is available from all
    refratio = 16;
    
    if(caseMag == 40)
        refratio = 16;
    elseif(caseMag == 20)
        refratio = 8;
    end
    
    
    % PROCESS
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
    %cases  = cat(1,cases,[{patient_i} lvls_mag]);

    %once we have lowest quality --> use 16 (=5x if 40x) to
    %begin data analysis
    start_lvl = find(round(lvls(:,3)) == refratio);
    caseData.setSeries(start_lvl-1);
    I1 = bfGetPlane(caseData,1);
    I2 = bfGetPlane(caseData,2);
    I3 = bfGetPlane(caseData,3);
    test_img(:,:,1) = I1;
    test_img(:,:,2) = I2;
    test_img(:,:,3) = I3;
    imshow(test_img);

    if(~exist([saveDir filesep tcga_id filesep 'voi']))
        mkdir([saveDir filesep tcga_id filesep 'voi'])
    end

    imshow(test_img);
    for i = 1:size(xy,2)
        if(~isempty(xy{1,i}))
            jroi = xy{1,i};
            jroi = jroi./refratio;
            hold on
            plot(jroi(:,1),jroi(:,2),'r','LineWidth',3);
        end
    end

    roi_final = xy;
    
    roi_final = xy;
    export_fig C:\tmp\tmp.png -native
    movefile(['C:\tmp\tmp.png'],[saveDir filesep tcga_id filesep 'voi' filesep savename '.png']);
    close all

    save([saveDir filesep tcga_id filesep 'voi' filesep savename '_points.mat'],'roi_final')
end
