clear
clc
close all
warning('OFF');


% Description: calculate origin point from XYZ lead and generate image of XYZ 
% median beat based on the modified baseline.
% SYNTAX : origin_point_algorithm(file_name,source_folder,image_folder,destination_folder)
% INPUT  : file_name          - file name for the matfile containing XYZ 
%                               median beat detected R peak indexed location.
%          source_folder      - the folder location where the file is placed.
%          image_folder       - the destination folder to place generated image.
%          destination_folder - where the matfile will be placed after 
%                               calculation is done.

% Origin Point Calculation Code V.1.1
% Authors:
%       Erick Andres Perez Alday, PhD, <perezald@ohsu.edu>
%       Annabel Li-Pershing, BS, <annabel2@pdx.edu>
% Principle Investigator:
%       Larisa Tereshchenko, MD, PhD, < tereshch@ohsu.edu > 
% Last update: August 2nd, 2018



% sampling frequency
prompt='What is the sampling frequency in Hz';
dlg_title='Sample Frequency';
fs_d = inputdlg(prompt,dlg_title,1);
fs=str2num(fs_d{1});


% load mat files
[file_name,path_name] = uigetfile('*','Select file for Median beat');

%% file import
file_ID          = strsplit(file_name,'.');
file_path        = fullfile(path_name,file_name);

A= load(file_path);

%% Create Results folder
Results_folder = [path_name, '/', 'Origin_point_results', '/'];
if (exist(Results_folder,'file') ~= 7)
    mkdir (Results_folder);
end



    %% Load data from matfile
%   A=load(strcat(source_folder,'/',file_name));
    XYZ_M=A.XYZ_M_T;
    Rx=A.Rpeaks(:,1);	
    disp(file_name)
    
    %% Calculate various windows analysis
    w=floor((mean(diff(Rx))*.8)-(372.69*(fs/1000)));

    %% Is the RR interval < 600ms (0.6*fs)
    if mean(diff(Rx))<(600*(fs/1000)) && (A.Label_beat~='S')
	cent=(600-260)*(fs/1000);
    else
	cent=(600-320)*(fs/1000);
    end


    w1_start=cent-w;
    w1_end=cent+w;

    w2_start=cent-(160*(fs/1000));
    w2_end=cent+(160*(fs/1000));

    if w1_start < w2_start
       w1_start=w2_start;
    end

    if w1_end > w2_end
        w1_end=w2_end;
    end

    %% cut xyz median beat based on the window
    xyz_w1 = XYZ_M(w1_start:w1_end,:);
    xyz_w2 = XYZ_M(w2_start:w2_end,:);

    %% calculate origin point with clustering algorithm ori_clustering
    [cluster_pt,~]=ori_clustering(xyz_w1);

    if isnan(cluster_pt)
        [cluster_pt,~]=ori_clustering(xyz_w2);
        Ori_pt(1) = cluster_pt+w2_start;
    else 
        Ori_pt(1) = cluster_pt+w1_start;  % convert the calculated index to real time window
    end


    % calculate the first origin point and the modified vector magnitude
    for ii=1:length(XYZ_M(:,1))
        VecMag_M1(ii)=norm([XYZ_M(Ori_pt(1),1)-XYZ_M(ii,1),XYZ_M(Ori_pt(1),2)-XYZ_M(ii,2),XYZ_M(Ori_pt(1),3)-XYZ_M(ii,3)]);
    end


    % calculate the sum of the absolute gradient on XYZ within w1
    [abs_grad_ori,~]=abs_grad_fun(xyz_w1);

    % calcaulate the second origin point
    Ori_pt2(1)=abs_grad_ori+w1_start;

    % calculate the second origin point and the modified vector magnitude
    VecMag_M2 = nan(1,length(XYZ_M(:,1)));
    for ii=1:length(XYZ_M(:,1))
        VecMag_M2(ii)=norm([XYZ_M(Ori_pt2(1),1)-XYZ_M(ii,1),XYZ_M(Ori_pt2(1),2)-XYZ_M(ii,2),XYZ_M(Ori_pt2(1),3)-XYZ_M(ii,3)]);
    end

    % calculate the area under the curve of the given angle
    Area_m1=trapz(VecMag_M1(w1_start:w1_end));
    Area_m2=trapz(VecMag_M2(w1_start:w1_end));

    % modify the XYZ median beat with the two origon point calculated
    XYZ_M1 = repmat(XYZ_M(Ori_pt(1),:),size(XYZ_M,1),1) - XYZ_M;
    XYZ_M2 = repmat(XYZ_M(Ori_pt2(1),:),size(XYZ_M,1),1) - XYZ_M;

    % calculate the sum of the absolute gradient on modified XYZ median beats
    data_grad_M1 = data_processing(XYZ_M1,'gradsum',5);
    data_grad_M2 = data_processing(XYZ_M2,'gradsum',5);

    % calcualte the average of the sum of the absolute gradient about the
    % origin points with 3 points on either side
    Sum_M1=mean(data_grad_M1(Ori_pt-3:Ori_pt+3));
    Sum_M2=mean(data_grad_M2(Ori_pt2-3:Ori_pt2+3));

    % decide which origin point to use based on the area under the curve and
    % the average gradient about the origin points
    if Area_m1 < Area_m2 && Sum_M1 < Sum_M2 
        dat_ori=Ori_pt;
        VecMag_ML=VecMag_M1;
    else
        dat_ori=Ori_pt2;
        VecMag_ML=VecMag_M2;
    end

    %% plotting
    VCG_fig = figure('visible','off');

    subplot(211)
    hold on
    plot(VecMag_ML,'b','Linewidth',1.5);
    plot([0 600],[0 0],'k:');
    hold off



    subplot(212)
    hold on;
    plot(XYZ_M(:,1),'Linewidth',1);
    plot(XYZ_M(:,2),'Linewidth',1);
    plot(XYZ_M(:,3),'Linewidth',1);
    ylabel('XYZ');
    hold off

    %% write variable to file and move file to destination
    saveas(VCG_fig,strcat(Results_folder,file_name(1:length(file_name)-4),'_origin.jpg'));
    copyfile(file_path,Results_folder);
    save(strcat(Results_folder,file_name(1:length(file_name)-4),'.mat'),'-append','VecMag_ML','dat_ori','-mat');




%% INCLUDED FUNCTIONS



function [ori_pt_idx,ori_val]=abs_grad_fun(xyz_median)

    % calculate the sum of the absolute gradient on XYZ within w1
    diff_w1 = data_processing(xyz_median,'gradsum',2);

    iter_num=0;
    [~,min_idx]=min(diff_w1);

    % find the point at which has the least angle variation at both side of the
    % window
    for ii=1:length(diff_w1)-1
        if abs(diff_w1(length(diff_w1)-ii)-diff_w1(min_idx))<0.1
            iter_num=length(diff_w1)-ii;

            break
        end

    end

    ori_pt_idx=iter_num;
    ori_val=xyz_median(iter_num);	

end


function [ori_pt_idx,ori_pt_cluster] = ori_clustering(xyz_median)

    %% Find the clusters of possible origin point candidate 
    % SYNTAX :  [ori_pt_idx,ori_pt_cluster] = ori_clustering(xyz_median)
    % INPUT  :  xyz_median  - XYZ median beat (segmented to a window) 
    %                         size(M x 3).
    % OUTPUT :  ori_pt_idx - the index of the calculated origin point
    %           ori_pt_cluster - origin point candidate clusters index

    %% Process data
    processed_data = data_processing(xyz_median,'norm_mvvar',10);
    [N,~,bin_label] = histcounts(processed_data,'BinWidth',10);
    xyz_norm = nan(size(xyz_median,1),1);

    for ii=1:length(xyz_median(:,1))
        xyz_norm(ii)=norm([xyz_median(ii,1),xyz_median(ii,2),xyz_median(ii,3)]);
    end

    %% Find origin point candidate clusters
    idx = find(N == max(N));
    candidate_idx = find(bin_label==idx(1));
    idx_diff = diff(candidate_idx);
    diff_idx = [0; find(idx_diff>=4)];

    if length(diff_idx) == 1
        % when there's only one point found
        ori_pt_cluster = candidate_idx;
        ori_pt_idx = round(median(candidate_idx));
    else
        % look for two largest clusters
        cluster_ele_count = zeros(length(diff_idx),1);

        for ii = 1:length(diff_idx)
            if ii == length(diff_idx)
                cluster_ele_count(ii) = sum(idx_diff(diff_idx(ii)+1:...
                                            length(idx_diff)));
            else
                cluster_ele_count(ii) = sum(idx_diff(diff_idx(ii)+1:...
                                            diff_idx(ii+1)-1));
            end
        end

        [~,I]=sort(cluster_ele_count,'descend');

        if I(1) == max(I)
            cluster_group1 = candidate_idx(diff_idx(I(1))+1:end);
        else
            cluster_group1 = candidate_idx(diff_idx(I(1))+1:...
                                           diff_idx(I(1)+1)-1);
        end

        if I(2) == max(I)
                    cluster_group2 = candidate_idx(diff_idx(I(2))+1:end);
        else
            cluster_group2 = candidate_idx(diff_idx(I(2))+1:...
                                           diff_idx(I(2)+1)-1);
        end

        %% find the best cluster by total average gradient of the clusters
        data_gradient_sum = data_processing(xyz_median,'gradsum',5);
        cluster_avg_slope(1) = mean(data_gradient_sum(cluster_group1));
        cluster_avg_slope(2) = mean(data_gradient_sum(cluster_group2));

        [~,bestgroup] = min(cluster_avg_slope);

        if bestgroup == 1
            ori_pt_cluster = cluster_group1;
            ori_pt_idx = round(median(cluster_group1));
        else
            ori_pt_cluster = cluster_group2;
            ori_pt_idx = round(median(cluster_group2));
        end
    end
end

function processed_data = data_processing(data,process,process_window)
    %% find the clusters of possible origin point candidate 
    % SYNTAX : processed_data = data_processing(data,process,process_window)
    % INPUT  : data            - data matrix with size(M x 3) [X,Y,X leads].
    %          process         - the method of processing data.
    %                            'norm_var': norn of the moving variance 
    %                            with given window.
    %                            'gradsum': sum of the absolute gradient
    %                            of the smoothed XYZ signals with given window.
    %          process_window  - window for specific process
    % OUTPUT : processed_data  - data processed based on specified process
    %                             and process window
    switch process
        case 'norm_mvvar'
            data_norm = nan(size(data,1),1);
            for ii=1:length(data(:,1))
                data_norm(ii)=norm([data(ii,1),data(ii,2),data(ii,3)]);
            end
            processed_data = movvar(data_norm,process_window);
        case 'gradsum'
            processed_data = abs(gradient(smooth(data(:,1),process_window)))...
                            +abs(gradient(smooth(data(:,2),process_window)))...
                            +abs(gradient(smooth(data(:,3),process_window)));

    end
end



