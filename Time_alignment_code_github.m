%
% Median beat calculation code V.1.1
% Erick Andres Perez Alday, PhD, <perezald@ohsu.edu>
% Larisa Tereshchenko, MD, PhD, < tereshch@ohsu.edu >
% Last update:  August 1st, 2018



clear 
clc
close all
warning('OFF');

% Description: calculate the median beat from XYZ leads and generate an image of XYZ
% median and Vectormagnitude beat based on norm of the XYZ leads
% Input:
%      1. sampling frequency in HZ
%      2. matlab file containing the following:
%        - XYZ_O(:,3): 10 seconds X, Y, and Z leads
%        - Rx: index of R peak on X lead
%        - Ry: index of R peak on Y lead
%        - Rz: index of R peak on Z lead

% Output:
%      2. Mat file containing the calculated median beat
%      4. 2D jpeg of median VectorMagniute plot





%% =========================== File and user input =============================
% sampling frequency

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

matfile          = matfile(file_path)

%% Create Results folder
Results_folder = [path_name, '/', 'Median_beat_results', '/'];

if (exist(Results_folder,'file') ~= 7)
    mkdir (Results_folder);
end


%% ===================== Initialize variables ========================
XYZ_O=[];
Rpeaks=[];
Beats_x_T=[];
Beats_y_T=[];
Beats_z_T=[];
XYZ_M_T=[];
VecMag_T=[];


%% ===================== load variables from .mat file ========================
XYZ_O=matfile.XYZ_O;
  
Rpeaks(:,1) = matfile.Rx;
Rpeaks(:,2) = matfile.Ry;
Rpeaks(:,3) = matfile.Rz;

R_no=length(matfile.Rx);

%% Total number of samples
Total_samples=length(XYZ_O(:,1));
%% Length of each beat is set to 1200 ms.
Beat_length=fs*1.2;

% =========================== Variable Calculation ============-================

%% Exclude the first if we do not have a complete first beat

if (Rpeaks(1,1)<(Beat_length/2)) || (Rpeaks(1,2)<(Beat_length/2)) || (Rpeaks(1,3)<(Beat_length/2))
    c_int=2;
else
    c_int=1;
end

%% Exclude the last if we do not have a complete last beat

if abs(Rpeaks(end-1,1)-Total_samples)<(Beat_length/2) || abs(Rpeaks(end-1,2)-Total_samples)<(Beat_length/2) || abs(Rpeaks(end-1,3)-Total_samples)<(Beat_length/2)
    c_fin=2;
else
    c_fin=1;
end

%% Calculate the absolute maximum first derivative adjacent to each R peak on the X lead.

I_dv=zeros(R_no,1);

for ii=c_int:R_no-c_fin
	[~,I_dv_temp]=max(abs(diff(XYZ_O(Rx(ii)-199:Rx(ii)+200,1))));
	I_dv(ii)=Rx(ii)-(199-I_dv_temp);
end


%% Separate the each beat in the X, Y, and Z leads. Each beat is centered on Maximum dv/dt detected.

bct=1;
for ii=c_int:R_no-c_fin
                           
        Beats_x_T(:,bct)=XYZ_O(I_dv(ii)-((Beat_length/2)-1):I_dv(ii)+(Beat_length/2),1);
        Beats_y_T(:,bct)=XYZ_O(I_dv(ii)-((Beat_length/2)-1):I_dv(ii)+(Beat_length/2),2);
        Beats_z_T(:,bct)=XYZ_O(I_dv(ii)-((Beat_length/2)-1):I_dv(ii)+(Beat_length/2),3);
        
        bct=bct+1;
end

%% Create the median beat of each X, Y, and Z leads 
XYZ_M_T(:,1)=median(Beats_x_T');
XYZ_M_T(:,2)=median(Beats_y_T');
XYZ_M_T(:,3)=median(Beats_z_T');


%% Calculate the Vector Magnitude (Euclidean norm)
for ii=1:length(XYZ_M_T(:,1)) 
    VecMag_T(ii)=norm(XYZ_M_T(ii,:));
end
     
timeM=((1:length(VecMag_T))/fs)*1000;
     
%% =============================================================================
% ================================== Plotting ==================================

    
ECG3 = figure('visible','off');
set (ECG3, 'Units', 'normalized', 'outerposition', [0,0,1,1]);

subplot(211)
plot(timeM,VecMag_T)
title('Time Alignment')
subplot(212)
plot(timeM,XYZ_M_T)
title('Median X, Y, and Z leads')

saveas(ECG3,strcat(Results_folder,file_ID{1},'_Median_beat.jpg'),'jpg');

%% =============================================================================

%% ================== Write Calculated Variables to .mat file ==================
Mat_folder = [Results_folder 'Mat Files' '/'];
if (exist(Mat_folder) == 0)
	mkdir (Mat_folder);
end

% Save to a new MAT file
Label_beat='N';
save(strcat(Mat_folder,file_ID{1},'.mat'),'XYZ_O','Rpeaks','VecMag_T','XYZ_M_T','Beats_x_T','Beats_y_T','Beats_z_T','Label_beat'); 

