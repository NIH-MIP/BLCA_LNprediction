function [] = makesmallerboxes_40x(img_file, tcga_id,mainDir,svsDir)

    %img file is from 5x folder
    warning('off','all')
    svssearch = strsplit(img_file,'_');
    svssearch = svssearch{1};
    svsfind = dir([svsDir filesep tcga_id filesep svssearch '*.svs']);
    svs_file = svsfind(1).name;
    caseData = bfGetReader([svsDir filesep tcga_id filesep svs_file]);
    caseMeta = caseData.getMetadataStore();
    caseMag = round(caseMeta.getObjectiveNominalMagnification(0,0).doubleValue(),2);

    if(caseMag == 40)
        refratio = 4;
    elseif(caseMag == 20)
        refratio = 1;
    end
    
        
    %load mask, find indices and put into ref ratio space
    maskedimg = 0;
    maskfind = dir([svsDir filesep tcga_id filesep svssearch '*ink-mask*.mat']);
    if(~isempty(maskfind))
        maskedimg = 1;
    end
    
    imgname = strsplit(img_file,'_');
    warning('off','all')
    mkdir([mainDir filesep '40x/orig' filesep tcga_id])
    mkdir([mainDir filesep '40x/norm' filesep tcga_id])
    mkdir([mainDir filesep '40x/decon' filesep tcga_id])
    mkdir([mainDir filesep 'TIL/orig' filesep tcga_id])
    
    box_label = imgname{3};
    boxparts = strsplit(box_label,'-');
    box_dims = [str2num(boxparts{2}) str2num(boxparts{3})];
    
    findbig = dir([mainDir filesep 'prelim' filesep tcga_id filesep 'norm' filesep imgname{1} '_' imgname{2} '_' imgname{3} '_*']);
    refbigbox = imread([findbig.folder filesep findbig.name]);
            
   
    imgSize = 300;
    counter20 = 1;
    counter40=1;
    tilcounter = 1;
    for i = 1:4
        for j = 1:4
                box_Val = ['xxbox' int2str(counter20)];
                ref20box= refbigbox((i-1)*imgSize+1:i*imgSize,(j-1)*imgSize+1:j*imgSize,:);
                bwref = rgb2gray(ref20box);

                dims_20 = [box_dims(1)+(i-1)*imgSize box_dims(2)+(j-1)*imgSize];
                
                actual_20 = refratio.*dims_20;
                caseData.setSeries(0);
                I1 = bfGetPlane(caseData,1,actual_20(2),actual_20(1),imgSize*refratio,imgSize*refratio);
                I2 = bfGetPlane(caseData,2,actual_20(2),actual_20(1),imgSize*refratio,imgSize*refratio);
                I3 = bfGetPlane(caseData,3,actual_20(2),actual_20(1),imgSize*refratio,imgSize*refratio);
                img_20(:,:,1) = I1;
                img_20(:,:,2) = I2;
                img_20(:,:,3) = I3;
                img20 = imresize(img_20,0.5);
                clear img_20
                
                if(maskedimg==1)
                    maskimg = zeros(size(bwref));
                    maskinds = find(double(bwref)./255 > 0.85);
                    maskimg(maskinds) = 1;
                    mask_big = imresize(maskimg,2,'method','nearest'); 
                    mask_big = cat(3,mask_big,mask_big,mask_big);
                    mask_inds=find(mask_big==1);
                    if(~isempty(mask_inds))
                        img20(mask_inds) = prctile(double(rgb2gray(img20)),99,'all');
                    end
                end
                
                if(length(find(double(bwref)./255 > 0.8))/(size(bwref,1)*size(bwref,2))<0.9)
                    [ NormMM ] = Norm(img20, ref20box, 'Reinhard');
                    NormMM = uint8(NormMM.*255);
                    stains = Deconvolve(NormMM, [], 0);
                    [s1, s2, s3] = PseudoColourStains(stains, []);
                    %bw
                    newRGB(:,:,1) =  rgb2gray(NormMM);
                    newRGB(:,:,2) =  rgb2gray(s1);
                    newRGB(:,:,3) =  rgb2gray(s2);
                    rgbout = newRGB;
                    clear newRGB
                else
                    NormMM = img20;
                    NormMM = uint8(NormMM.*255);
                    stains = Deconvolve(NormMM, [], 0);
                    [s1, s2, s3] = PseudoColourStains(stains, []);
                    %bw
                    newRGB(:,:,1) =  rgb2gray(NormMM);
                    newRGB(:,:,2) =  rgb2gray(s1);
                    newRGB(:,:,3) =  rgb2gray(s2);
                    rgbout = newRGB;
                    clear newRGB
                end
                
                for x20 = 1:2
                    for y20 = 1:2
                        box40 = ['xxxbox' int2str(counter40)];
                        saveimg = img20((x20-1)*imgSize+1:x20*imgSize,(y20-1)*imgSize+1:y20*imgSize,:);
                        normsave = NormMM((x20-1)*imgSize+1:x20*imgSize,(y20-1)*imgSize+1:y20*imgSize,:);
                        rgbsave = rgbout((x20-1)*imgSize+1:x20*imgSize,(y20-1)*imgSize+1:y20*imgSize,:);
                        ws_check = double(rgb2gray(normsave))./255;
                        findwhite = find(ws_check>0.8);
                        if(length(findwhite)/(imgSize*imgSize) < 0.7)
                            ws_Val = sprintf('%.2f',length(findwhite)/(imgSize*imgSize));
                            savename40 = strjoin([imgname(1:3),box40,ws_Val,imgname(5)],'_');
                            imwrite(saveimg,[mainDir filesep '40x/orig' filesep tcga_id filesep savename40]);
                            imwrite(normsave,[mainDir filesep '40x/norm' filesep tcga_id filesep savename40]);
                            imwrite(rgbsave,[mainDir filesep '40x/decon' filesep tcga_id filesep savename40]);
                        end
                        counter40 = counter40+1;
                    end
                end

                tbxSz = 100;
                for tilli = 1:3
                    for tillj = 1:3
                            big_dims = [dims_20(1)+(tillj-1)*tbxSz dims_20(2)+(tilli-1)*tbxSz];
                            actual_dims = refratio.*big_dims;                          
                          
                            caseData.setSeries(0);
                            I1 = bfGetPlane(caseData,1,actual_dims(2),actual_dims(1),tbxSz*refratio,tbxSz*refratio);
                            I2 = bfGetPlane(caseData,2,actual_dims(2),actual_dims(1),tbxSz*refratio,tbxSz*refratio);
                            I3 = bfGetPlane(caseData,3,actual_dims(2),actual_dims(1),tbxSz*refratio,tbxSz*refratio);
                            in_img(:,:,1) = I1;
                            in_img(:,:,2) = I2;
                            in_img(:,:,3) = I3;
                            lvl_img = imresize(in_img,0.5);
                            clear in_img
                            
                            for smx = 1:2
                                for smy = 1:2
                                    saveimg = lvl_img((smx-1)*tbxSz+1:smx*tbxSz,(smy-1)*tbxSz+1:smy*tbxSz,:);
                                    ws_check = double(rgb2gray(saveimg))./255;
                                    findwhite = find(ws_check>0.8);
                                    if(length(findwhite)/(10000) < 0.9)
                                        tilname = ['til-' int2str(tilcounter)];
                                        ws_Val = sprintf('%.2f',length(findwhite)/(10000));
                                        savenameTIL = strjoin([imgname(1:3),tilname,ws_Val,imgname(5)],'_');
                                        imwrite(saveimg,[mainDir filesep 'TIL/orig' filesep tcga_id filesep savenameTIL]);
                                    end
                                tilcounter = tilcounter+1;
                                end
                            end

                    end
                end   
            counter20 = counter20 + 1;
        end
    end
end