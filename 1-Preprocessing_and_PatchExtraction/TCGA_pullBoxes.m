
function [] = TCGA_pullBoxes(img_file,mainDir,patchDir,outcomeFile)
% pull in magnification
    disp(img_file)
    case_id = img_file(1:12);
    tcga_id = case_id;
    
    saveDir = [patchDir filesep 'prelim'];
    
    mkdir([saveDir filesep tcga_id]);
    mkdir([saveDir filesep tcga_id filesep 'norm']);
    mkdir([saveDir filesep tcga_id filesep 'decon']);
    mkdir([saveDir filesep tcga_id filesep 'orig']);

    save_id = img_file;

    load(outcomeFile)
    case_find = find(strcmpi(outcomedata(:,1),case_id)>0);
    outcome = outcomedata{case_find,2};
    disp([img_file ': ' outcome])

            
    svsfind = dir([mainDir filesep tcga_id filesep img_file '*.svs']);
    svs_file = svsfind(1).name;
    caseData = bfGetReader([mainDir filesep tcga_id filesep svs_file]);
    caseMeta = caseData.getMetadataStore();
    caseMag = round(caseMeta.getObjectiveNominalMagnification(0,0).doubleValue(),2);
    voifind = dir([mainDir filesep tcga_id filesep 'voi' filesep img_file '*_points.mat']);
    load([voifind(1).folder filesep voifind(1).name])
    load([mainDir filesep tcga_id filesep img_file '.mat'])

    %we are expecting all to be scanned the same, at 40x
    % selecting other levels based on what is available from all
    if(caseMag == 40)
        refratio = 4;
    elseif(caseMag == 20)
        refratio = 1;
    end
    boxSz = 1200;


    load([mainDir filesep tcga_id filesep save_id '_norm.mat'])
    stains = Deconvolve(NormMM, [], 0);
    [s1, s2, s3] = PseudoColourStains(stains, []);


    %bw_big = double(rgb2gray(NormMM))./255;
    %bw_val = prctile(bw_big,99.5,[1 2]);
    img_decon(:,:,1) = rgb2gray(NormMM);
    img_decon(:,:,2) = rgb2gray(s1);
    img_decon(:,:,3) = rgb2gray(s2);
    save([mainDir filesep tcga_id filesep save_id '_decon.mat'],'img_decon','-v7.3')
    smallsave = imresize(img_decon,0.1);
    imwrite(smallsave,[mainDir filesep tcga_id filesep save_id '_decon.jpeg']);

    % pull over each roi at a time
        for i = 1:size(roi_final,2)
            jroi = roi_final{1,i};
            jroi = jroi./refratio;
            %hold on
            %plot(jroi(:,1),jroi(:,2),'r','LineWidth',3); 

            %pull roi indices    
            y_min = round(min(jroi(:,1))); y_max = round(max(jroi(:,1)));
            x_min = round(min(jroi(:,2))); x_max = round(max(jroi(:,2)));
            if(x_max>size(NormMM,1)) x_max=size(NormMM,1); end 
            if(x_min<1) x_min=1; end
            if(y_max>size(NormMM,2)) y_max=size(NormMM,2); end 
            if(y_min<1) y_min=1; end

            [xmin, xmax, xstride] = findminmax(x_min, x_max, size(NormMM,1), boxSz);
            [ymin, ymax, ystride] = findminmax(y_min, y_max, size(NormMM,2), boxSz);

            jmask = poly2mask(jroi(:,1),jroi(:,2),size(NormMM,1),size(NormMM,2));
            x_samp = ceil((xmax-xmin)/boxSz);
            y_samp = ceil((ymax-ymin)/boxSz);


            warning('off','all')
            box_count = 1;
            for xi = 1:x_samp
                for yi = 1:y_samp
                    x_loc = xmin + (xi-1)*boxSz - (xi-1)*ceil(xstride);
                    y_loc = ymin + (yi-1)*boxSz - (yi-1)*ceil(ystride);
                    sm_mask = jmask(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1);
                    perc_in = length(find(sm_mask>0))/(boxSz*boxSz);
                    if(perc_in > 0.1)
                        norm_crop = NormMM(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1,:);
                        bw_img = double(rgb2gray(norm_crop))./(255);
                        find_white = find(bw_img > 0.8);
                        prop_white = length(find_white)/(size(bw_img,1)*size(bw_img,2)*size(bw_img,3));
                        if(prop_white < 0.9)
                            %return
                            imwrite(norm_crop,[saveDir filesep tcga_id filesep 'norm' filesep save_id '_roi' int2str(i) '_box' int2str(box_count) '-' int2str(x_loc) '-' int2str(y_loc) '_' sprintf('%.2f',perc_in) '_' sprintf('%.2f',prop_white) '_' outcome '.jpeg']);
                            decon_crop = img_decon(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1,:);
                            imwrite(decon_crop,[saveDir filesep tcga_id filesep 'decon' filesep save_id '_roi' int2str(i) '_box' int2str(box_count) '-' int2str(x_loc) '-' int2str(y_loc) '_' sprintf('%.2f',perc_in) '_' sprintf('%.2f',prop_white) '_' outcome '.jpeg']);
                            orig_crop = mimg(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1,:);
                            imwrite(orig_crop,[saveDir filesep tcga_id filesep 'orig' filesep save_id '_roi' int2str(i) '_box' int2str(box_count) '-' int2str(x_loc) '-' int2str(y_loc) '_' sprintf('%.2f',perc_in) '_' sprintf('%.2f',prop_white) '_' outcome '.jpeg']);
                            disp([tcga_id ' box ' int2str(box_count) ' ' int2str(x_loc) ' ' int2str(y_loc)])
                            box_count = box_count + 1;
                        end
                    else
                        norm_crop = NormMM(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1,:);
                        bw_img = double(rgb2gray(norm_crop))./(255);
                        find_white = find(bw_img > 0.8);
                        prop_white = length(find_white)/(size(bw_img,1)*size(bw_img,2)*size(bw_img,3));
                        if(prop_white < 0.9)
                            %return
                            imwrite(norm_crop,[saveDir filesep tcga_id filesep 'norm' filesep save_id '_roi' int2str(i) '_box' int2str(box_count) '-' int2str(x_loc) '-' int2str(y_loc) '_' sprintf('%.2f',perc_in) '_' sprintf('%.2f',prop_white) '_' outcome '.jpeg']);
                            decon_crop = img_decon(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1,:);
                            imwrite(decon_crop,[saveDir filesep tcga_id filesep 'decon' filesep save_id '_roi' int2str(i) '_box' int2str(box_count) '-' int2str(x_loc) '-' int2str(y_loc) '_' sprintf('%.2f',perc_in) '_' sprintf('%.2f',prop_white) '_' outcome '.jpeg']);
                            orig_crop = mimg(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1,:);
                            imwrite(orig_crop,[saveDir filesep tcga_id filesep 'orig' filesep save_id '_roi' int2str(i) '_box' int2str(box_count) '-' int2str(x_loc) '-' int2str(y_loc) '_' sprintf('%.2f',perc_in) '_' sprintf('%.2f',prop_white) '_' outcome '.jpeg']);
                            disp([tcga_id ' box ' int2str(box_count) ' ' int2str(x_loc) ' ' int2str(y_loc)])
                            box_count = box_count + 1;
                        end
                    end
                end              
            end
        end
      
      %make sure the above steps get saved in prelim  
      imglist = dir([saveDir filesep tcga_id filesep 'decon' filesep '*.jpeg']);
      for j = 1:size(imglist,1)
        makesmallerboxes(imglist(j).name, tcga_id,patchDir)
      end
      
      %figure out what this needs to point to
      imglist = dir([patchDir filesep '5x' filesep 'orig' filesep tcga_id filesep '*.jpeg']);
      for j = 1:size(imglist,1)
        makesmallerboxes_40x(imglist(j).name, tcga_id,patchDir,imgDir)
      end      
      
      imglist = dir([saveDir filesep tcga_id filesep 'decon' filesep '*.jpeg']);
      
end


