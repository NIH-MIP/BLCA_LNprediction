ptlist = dir('/path/to/dataset/maps/dataset//matfiles/TCGA*40x_prob*');
pt_out = [];
saveDir = '/path/to/save/predictions/stats_save';
ptorig = ptlist;
ptlist = [];
for pti = 1:numel(ptorig)
    ptname = strsplit(ptorig(pti).name,'_');
    if(~exist([saveDir filesep ptname{1} '.mat']))
        ptlist = cat(1,ptlist,ptorig(pti));
    end
end

%parpool(20)

for pti = 1:numel(ptlist)
    ptname = strsplit(ptlist(pti).name,'_');
  %if(~exist([saveDir filesep ptname{1} '.mat']))
  try
    disp(['starting ' ptlist(pti).name(1:12)])
    img40 = load([ptlist(pti).folder filesep ptlist(pti).name]);
    img20 = load([ptlist(pti).folder filesep strrep(ptlist(pti).name,'40x','20x')]);
    img10 = load([ptlist(pti).folder filesep strrep(ptlist(pti).name,'40x','10x')]);
    img5  = load([ptlist(pti).folder filesep strrep(ptlist(pti).name,'40x','5x')]);
    imgTIL = load([ptlist(pti).folder filesep strrep(ptlist(pti).name,'40x','TIL')]);
    if(contains(ptlist(pti).name,'_pos'))
        outcome = 1;
    else
        outcome = 0;
    end
    
    img40 = img40.out_map_40x; img20 = img20.out_map_20x; img10 = img10.out_map_10x; img5 = img5.out_map_5x; imgTIL = imgTIL.out_map_TIL;
    avg_prob = (img40 + img20 + img10 + img5)./4;
    three_prob = (img20 + img10 + img5)./3;
    tumorsz = numel(find(avg_prob>0));
    prod_prob = img40.*img20.*img20.*img5;
    avg_weighted = (2*((img5+4*img10)./5)+((img20+img40)./2))/3;
    bin5 = img5; bin5(find(bin5<0.5))=0; bin5(find(bin5>=0.5))=1;
    bin10 = img10; bin10(find(bin10<0.5))=0; bin10(find(bin10>=0.5))=1;
    bin20 = img20; bin20(find(bin20<0.5))=0; bin20(find(bin20>=0.5))=1;
    prod_weighted = (bin5+bin10+bin20)./3;
    imgTIL(find(avg_prob==0)) = 0;
    tumorsz = numel(find(avg_prob>0));

    
    % outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   avg_probs w varying threshold: 0.4 0.5 0.6           %
    %   avg weighted probs w varying threshold: 0.4 0.5 0.6  %
    %   prod weighted probs w varying threshold: 0.2 0.3 0.4 % 
    %   product probs w varying threshold: 0.05 0.10 0.15    %
    %   5x: 0.4 0.5 0.6                                      %
    %   10x: 0.4 0.5 0.6                                     %
    %   20x: 0.4 0.5 0.6                                     %
    %   40x: 0.4 0.5 0.6                                     %
    %   TIL alone: 0.5 0.75                                  %
    %   TIL dbscan: 0.5 0.75                                 % 
    %   (avg_prob+TIL) 0.5 0.75                              %
    %   (avg_weight+TIL) 0.5 0.75                            %
    %   (prod_weight+TIL) 0.5 0.75                           %
    %   (prod_prob+TIL) 0.5 0.75                             %
    %                                                        %
    % use cluster output and calc distance metrics:          %
    %   avg_probs (0.4 0.5 0.6) + TIL 0.5)                   %
    %   hybrid (0.4 0.5 0.6) + TIL 0.5)                      %
    %   product (0.05 0.15 0.25) + TIL 0.5)                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %fill in cell array with variable name and threshold and points and idx and stats
    %aggregate all stats
    
    pt_data = {};
    counter = 1;
    
    % individual probabilities + all other probability maps
    imgsets = {'avg_prob','prod_prob','avg_weighted','prod_weighted','three_prob','img5','img10','img20','img40'};
    imgdata = [{avg_prob},{prod_prob},{avg_weighted},{prod_weighted},{three_prob},{img5},{img10},{img20},{img40}];
    tilinds = find(imgTIL>0.5);
    [til1, til2] = ind2sub(size(imgTIL),tilinds);
    tilsubs = cat(2,til1,til2);
    for imgi = 1:numel(imgsets)
        for thri = 1:3
            disp([ptlist(pti).name(1:12) '     ... starting ' imgsets{imgi} ' ' int2str(thri) ' of 3'])
            if(imgi == 2) %set different threshold for prod_weighted
                thresh = 0.05+0.05*(thri-1);
            elseif(imgi==4)
                thresh = 0+1*(thri-1)/3;
            else
                thresh = 0.4+0.1*(thri-1);
            end
            %eval(['img = ' imgsets{imgi} ';'])
            img = imgdata{imgi};
            [idx_final, coords_final] = make_patches(200,100,img,thresh,0);
            %if(thri != 1)
            [data_out] = pull_stats(img,coords_final,idx_final);
            pt_data{counter,1} = imgsets{imgi};
            pt_data{counter,2} = thresh;
            pt_data{counter,3} = idx_final;
            pt_data{counter,4} = coords_final;
            pt_data{counter,5} = data_out;
            if(imgi<6) %we dont want to compare to TILs for individual probaility maps
               disp([ptlist(pti).name(1:12)  '           + comparison to TIL'])
               [data_out] =  DL_v_TIL(img,coords_final,idx_final,tilsubs,tumorsz);
               pt_data{counter,6} = data_out;
            end
            counter=counter+1;
        end
    end

    % probabilities avergaed with TIL
    imgsets = {'avg_prob','avg_weighted','prod_weighted','prod_prob','three_prob'};
    imgdata = [{avg_prob},{avg_weighted},{prod_weighted},{prod_prob},{three_prob}];
    for imgi = 1:numel(imgsets)
        for thri = 1:2
            disp([ptlist(pti).name(1:12) '     ... starting ' imgsets{imgi} '+TIL ' int2str(thri) ' of 3'])
            if(imgi == 4) %different threshold for prod_prob
                thresh = 0.25+0.25*(thri-1);
            else
                thresh = 0.5+0.25*(thri-1);
            end
            img = imgdata{imgi};
            %eval(['img = ' imgsets{imgi} ';'])
            img_new = img + imgTIL;
            [idx_final, coords_final] = make_patches(200,100,img_new,thresh,0);
            [data_out] = pull_stats(img,coords_final,idx_final);
            pt_data{counter,1} = imgsets{imgi};
            pt_data{counter,2} = thresh;
            pt_data{counter,3} = idx_final;
            pt_data{counter,4} = coords_final;
            pt_data{counter,5} = data_out;
            counter=counter+1;
        end
    end
        
    %TIL dbscan
    for thri = 1:2
        disp([ptlist(pti).name(1:12) '     ... starting TIL ' int2str(thri) ' of 2'])
        thresh = 0.5+0.25*(thri-1);
        [idx_final, coords_final] = make_patches(200,100,imgTIL,thresh,1);
        [data_out] = til_stats(imgTIL,coords_final,idx_final,tumorsz);
        pt_data{counter,1} = 'imgTIL';
        pt_data{counter,2} = thresh;
        pt_data{counter,3} = idx_final;
        pt_data{counter,4} = coords_final;
        pt_data{counter,5} = data_out;
        counter=counter+1;
    end
    
    %TIL simple threshold
    for thri = 1:2
        disp([ptlist(pti).name(1:12) '     ... starting simple TIL ' int2str(thri) ' of 2'])
        thresh = 0.5+0.25*(thri-1);
        tilinds = find(imgTIL>thresh);
        [til1, til2] = ind2sub(size(imgTIL),tilinds);
        tilsubs = cat(2,til1,til2);
        til_idx = ones(numel(tilinds),1);
        [data_out] = til_simple_stats(imgTIL,tilsubs,til_idx,tumorsz);
        pt_data{counter,1} = 'TILsimple';
        pt_data{counter,2} = thresh;
        pt_data{counter,3} = til_idx;
        pt_data{counter,4} = tilsubs;
        pt_data{counter,5} = data_out;
        counter=counter+1;
    end
       
    ptname = strsplit(ptlist(pti).name,'_');
    callSaveFun([saveDir filesep ptname{1} '.mat'],pt_data)

  catch ME
  end
  %end
    %til_stats(img,grp_pts,idx,tumorsz)
    %[data_out] = DL_v_TIL(img,dl_pts,dl_idx,til_pts,tumorsz);
    
%     figure, imshow(avg_prob)
%     hold on
%     h = scatter(coords_final(:,2),coords_final(:,1),'r.');
%     h.MarkerEdgeAlpha = 0.25;
%     hold on
%     h = scatter(tilsubs(:,2),tilsubs(:,1),'b.');
%     h.MarkerEdgeAlpha = 0.25;

%     figure, gscatter(coords_final(:,2),coords_final(:,1),idx_final);
%     set(gca,'YDir','reverse');    
end


%add dbscan flag
function [idx_final, coords_final] = make_patches(subSize,stride,img,thresh,tilflag)
    %subSize = 200;
    imgSize = size(img);
    %stride = 100;
    
    leftoverx = round((ceil(imgSize(1)/subSize)-imgSize(1)/subSize)*subSize);
    leftovery = round((ceil(imgSize(2)/subSize)-imgSize(2)/subSize)*subSize);
    
    img = padarray(img,[leftoverx,leftovery],0,'post');
    imgSize = size(img);
    %find number we can fill space
    num_subs_updown = floor((imgSize(1)-subSize)/stride)+1;
    num_subs_leftright = floor((imgSize(2)-subSize)/stride)+1;
    x=1:num_subs_updown;
    y=1:num_subs_leftright;
    [X,Y] = meshgrid(x,y);
    
    master_idx = [];
    master_coords = [];
    for indi = 1:numel(X)
        sub_img = img(stride*(X(indi)-1)+1:stride*(X(indi)-1)+subSize,stride*(Y(indi)-1)+1:stride*(Y(indi)-1)+subSize);
        if(tilflag == 0)
            [coordsi,idxi]=run_dbscan(sub_img,thresh);
        else
            [coordsi,idxi]=til_dbscan(sub_img,thresh);
        end
        if(~isempty(coordsi))
            coords_transf = [coordsi(:,1)+stride*(X(indi)-1) coordsi(:,2)+stride*(Y(indi)-1)];
            grp_idx = find(idxi>0);
            if(isempty(master_idx))
                master_max = 0;
            else
                master_max = max(master_idx);
            end
            coords_keep = coords_transf(grp_idx,:);
            idx_keep = idxi(grp_idx,:)+master_max;
            master_coords = cat(1,master_coords,coords_keep);
            master_idx = cat(1,master_idx,idx_keep);
        end
    end
    
    %aggregate all together if regions are overlapping then join them
    if(~isempty(master_coords))
        master_ind = sub2ind(size(img),master_coords(:,1),master_coords(:,2));
        new_id = 1;
        new_coords = [];
        new_idx = [];
        already_handled = [];
        for i = 1:length(unique(master_idx))
            if(~ismember(i,already_handled))
                %find all points associated with region i
                f_i = find(master_idx == i);
                f_others = find(master_idx ~= i);
                idx_i = master_idx(f_i,:);
                idx_other = master_idx(f_others,:);
                coords_i = master_coords(f_i,:);
                coords_other = master_coords(f_others,:);

                [in_i,loci] = ismember(coords_other,coords_i,'rows');
                idx_in_other = idx_other(find(in_i));
                idx_in_other = unique(idx_in_other);
                [idx_in,loci] = ismember(idx_other,idx_in_other);

                new_matched = cat(1,coords_i,coords_other(idx_in,:));
                match_search = new_matched;
                stop_Add = 0;
                %counter = 1;
                while(stop_Add<1)
                    %disp(counter)
                    new_added = unique(match_search,'rows');
                    length_prior = size(new_added,1);
                    [search_again, search_b]= ismember(coords_other,new_added,'rows');
                    idx_again = idx_other(find(search_again));
                    idx_again = unique(idx_again);
                    [idx_in,loci] = ismember(idx_other,idx_again);
                    match_add = cat(1,match_search,coords_other(idx_in,:));
                    match_add = unique(match_add,'rows');
                    length_added = size(match_add,1);
                    if(length_added > length_prior)
                        match_search = match_add;
                    else
                        stop_Add = 1;
                    end
                end
                new_matched = match_search;
                
                new_assign = new_id*ones(size(new_matched,1),1);
                new_coords = cat(1,new_coords,new_matched);
                new_idx = cat(1,new_idx,new_assign);

                new_id = new_id + 1;
                already_handled = cat(1,already_handled,i,idx_again);
            end       
        end

        [coords_final, ai, ci] = unique(new_coords,'rows');
        idx_final = new_idx(ai);
    else
        coords_final = [];
        idx_final = [];
    end
end


function [pycoords,idx] = run_dbscan(img,thresh)
    init_inds=find(img>thresh);
    if(numel(init_inds)>0)
        [pycoords(:,1), pycoords(:,2)] = ind2sub(size(img),init_inds);
        D = pdist2(pycoords,pycoords);
        [idx, corepts] = dbscan(D,10,100,'Distance','precomputed');
    else
        idx = [];
        pycoords = [];
    end
end

function [pycoords,idx] = til_dbscan(img,thresh)
    init_inds=find(img>thresh);
    if(numel(init_inds)>0)
        [pycoords(:,1), pycoords(:,2)] = ind2sub(size(img),init_inds);
        D = pdist2(pycoords,pycoords);
        [idx, corepts] = dbscan(D,10,10,'Distance','precomputed');
    else
        idx = [];
        pycoords = [];
    end
end

%probability and distance information
function [data_out] = pull_stats(img,grp_pts,idx)
    %grp_idx = find(idx>0);
    
    if(~isempty(idx))
        tumorsz = numel(find(img>0));
        groups = unique(idx);
        %grp_pts = pycoords(grp_idx,:);
        grp_inds = sub2ind(size(img),grp_pts(:,1),grp_pts(:,2));

        group_dat = zeros(numel(groups),10);
        for i = 1:numel(groups)
            grpi = find(idx == groups(i));
            grpi_pts = grp_pts(grpi,:);
            grpi_inds = sub2ind(size(img),grpi_pts(:,1),grpi_pts(:,2));
            oi = find(idx ~= groups(i));
            oi_pts = grp_pts(oi,:);
            if(~isempty(oi))
                Dmin = pdist2(oi_pts,grpi_pts,'euclidean','Smallest',1);
            end
            Smin = pdist2(grpi_pts,grpi_pts,'euclidean','Smallest',2);
            Smin = Smin(2,:);
            Smax = pdist2(grpi_pts,grpi_pts,'euclidean','Largest',1);
            Sw = pdist2(grpi_pts,grpi_pts,'euclidean');
            group_dat(i,1) = numel(grpi);
            group_dat(i,2) = median(img(grpi_inds));
            group_dat(i,3) = std(img(grpi_inds))/mean(img(grpi_inds));
            group_dat(i,4) = std(img(grpi_inds))/sqrt(numel(grpi));
            if(~isempty(oi))
                group_dat(i,5) = min(Dmin(:));
                group_dat(i,6) = mean(Dmin(:));
            end
            group_dat(i,7) = sum(Smin,'all');
            group_dat(i,8) = sum(Smax,'all');
            group_dat(i,9) = sum(Sw,'all');
            group_dat(i,10) = max(Smax(:));
        end

        total_elements = sum(group_dat(:,1));
        burden = total_elements/tumorsz;
        max_ind = find(group_dat(:,1)==max(group_dat(:,1)));
        max_elements = max(group_dat(:,1));
        max_burd = max_elements/total_elements;
        maxprob = group_dat(max_ind,2);
        if(numel(max_ind)>1)
            maxprob = max(maxprob);
        end
        med_all = median(img(grp_inds)); %median of all cluster points
        cov_all = std(img(grp_inds))/mean(img(grp_inds));
        se_all = std(img(grp_inds))/sqrt(total_elements);
        mean_meds = mean(group_dat(:,2));
        mean_cov = mean(group_dat(:,3));
        mean_se = mean(group_dat(:,4));
        mean_minD = mean(group_dat(:,5));
        mean_meanD = mean(group_dat(:,6));
        mean_maxD = mean(group_dat(:,10));
        C_index = (sum(group_dat(:,9))-sum(group_dat(:,7)))/(sum(group_dat(:,8))-sum(group_dat(:,7)));
        
        data_out = cat(2,total_elements,burden,max_elements,max_burd,maxprob,med_all,cov_all,se_all,mean_meds,mean_cov,mean_se,mean_minD,mean_meanD,mean_maxD,C_index);
    else
        data_out = cat(2,[0],[0],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN]);
    end

end

%distance only info
function [data_out] = til_stats(img,grp_pts,idx,tumorsz)
    if(~isempty(idx))
        groups = unique(idx);
        grp_inds = sub2ind(size(img),grp_pts(:,1),grp_pts(:,2));

        group_dat = zeros(numel(groups),10);
        for i = 1:numel(groups)
            grpi = find(idx == groups(i));
            grpi_pts = grp_pts(grpi,:);
            grpi_inds = sub2ind(size(img),grpi_pts(:,1),grpi_pts(:,2));
            oi = find(idx ~= groups(i));
            oi_pts = grp_pts(oi,:);
            if(~isempty(oi))
                Dmin = pdist2(oi_pts,grpi_pts,'euclidean','Smallest',1);
            end
            Smin = pdist2(grpi_pts,grpi_pts,'euclidean','Smallest',2);
            Smin = Smin(2,:);
            Smax = pdist2(grpi_pts,grpi_pts,'euclidean','Largest',1);
            Sw = pdist2(grpi_pts,grpi_pts,'euclidean');
            group_dat(i,1) = numel(grpi);
            if(~isempty(oi))
                group_dat(i,2) = min(Dmin(:));
                group_dat(i,3) = mean(Dmin(:));
            end
            group_dat(i,4) = sum(Smin,'all');
            group_dat(i,5) = sum(Smax,'all');
            group_dat(i,6) = sum(Sw,'all');
            group_dat(i,7) = max(Smax(:));
        end

        total_elements = sum(group_dat(:,1));

        burden = total_elements/tumorsz;
        max_ind = find(group_dat(:,1)==max(group_dat(:,1)));
        max_elements = max(group_dat(:,1));
        max_burd = max_elements/total_elements;
        maxprob = group_dat(max_ind,2);
        if(numel(max_ind)>1)
            maxprob = max(maxprob);
        end
        mean_minD = mean(group_dat(:,2));
        mean_meanD = mean(group_dat(:,3));
        mean_maxD = mean(group_dat(:,7));
        C_index = (sum(group_dat(:,6))-sum(group_dat(:,4)))/(sum(group_dat(:,5))-sum(group_dat(:,4)));
        
        data_out = cat(2,total_elements,burden,max_elements,max_burd,maxprob,mean_minD,mean_meanD,mean_maxD,C_index);
    else
        data_out = cat(2,[0],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN]);
    end
end

%distance only info
function [data_out] = til_simple_stats(img,grp_pts,idx,tumorsz)
    if(numel(idx)>1)
        Smin = pdist2(grp_pts,grp_pts,'euclidean','Smallest',2);
        Smin = Smin(2,:);
        Smax = pdist2(grp_pts,grp_pts,'euclidean','Largest',1);
        Sw = pdist2(grp_pts,grp_pts,'euclidean');

        min_sum = sum(Smin,'all');
        max_sum = sum(Smax,'all');
        all_sum = sum(Sw,'all');
        
        meanD = mean(Sw(:));
        mean_minD = mean(Smin,'all');
        total_elements = numel(idx);
        burden = total_elements/tumorsz;

        C_index = (all_sum-min_sum)/(max_sum-min_sum);
        
        data_out = cat(2,total_elements,burden,meanD,mean_minD,C_index);
    else
        data_out = cat(2,[0],[NaN],[NaN],[NaN],[NaN]);
    end
end

%distance
function [data_out] = DL_v_TIL(img,dl_pts,dl_idx,til_pts,tumorsz)
    groups = unique(dl_idx);

    if(~isempty(dl_pts) && ~isempty(til_pts))
        group_dat = zeros(numel(groups),2);
        for i = 1:numel(groups)
            grpi = find(dl_idx == groups(i));
            grpi_pts = dl_pts(grpi,:);
            Dmin = pdist2(til_pts,grpi_pts,'euclidean','Smallest',1);
            Dall = pdist2(til_pts,grpi_pts,'euclidean');
            group_dat(i,1) = min(Dmin(:));
            group_dat(i,2) = mean(Dall(:));
        end

        grp_inds = sub2ind(size(img),dl_pts(:,1),dl_pts(:,2));
        til_inds = sub2ind(size(img),til_pts(:,1),til_pts(:,2));
        either_inds = cat(1,grp_inds,til_inds);
        either_inds = unique(either_inds);
        both_inds = intersect(grp_inds,til_inds);
        dl_only = setdiff(grp_inds,til_inds);
        til_only = setdiff(til_inds,grp_inds);
        burden = numel(either_inds)/tumorsz;
        burden_int = numel(both_inds)/tumorsz;
        burden_til = numel(til_only)/tumorsz;
        burden_dl = numel(dl_only)/tumorsz;
        min_minD = min(group_dat(:,1));
        mean_minD = mean(group_dat(:,1));
        min_meanD = mean(group_dat(:,2));
        mean_meanD = mean(group_dat(:,2));
        data_out = cat(2,burden,burden_int,burden_til,burden_dl,min_minD,mean_minD,min_meanD,mean_meanD);
    else
        data_out = cat(2,[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN]);
    end
end

% function out = getVarName(var)
%     out = inputname(1);
% end

function [] = callSaveFun(filename,slidedata)
    save(filename,'slidedata');
end