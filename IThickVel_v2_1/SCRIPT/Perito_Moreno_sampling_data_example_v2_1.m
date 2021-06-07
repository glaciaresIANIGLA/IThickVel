clear all
format short
clc

%% Header

% Script to use ithicknvel_v2_1.m script with Perito Moreno sampling data
% 
% Sampling data :  'Perito_Moreno_sampling_data.tif'
% Grid or band 1 of the tif file = smoothed slope of Perito Moreno glacier
% (see Zorzut and others, 2020)
% Grid or band 2 of the tif file = Mean surface velocity of Perito Moreno  glacier by Mouginot and Rignot (2015).
% Grid or band 3 of the tif file = Surface elevation of Perito Moreno
% glacier from SRTM v4.1 (Jarvis and others, 2008(

% Reference:
%Jarvis, A., H.I. Reuter, A. Nelson, E. Guevara, 2008, Hole-filled SRTM for the globe Version 4, available from the CGIAR-CSI SRTM 90m Database (http://srtm.csi.cgiar.org).
%Mouginot, J. and Rignot, E.: Ice motion of the Patagonian Icefields of South America: 1984–2014, Geophys. Res. Lett., 42, 2014GL062661, https://doi.org/10.1002/2014GL062661, 2015.
%Zorzut, V., Ruiz, L., Rivera, A., Pitte, P., Villalba, R., and Medrzycka, D.: Slope estimation influences on ice thickness inversion models: a case study for Monte Tronador glaciers, North Patagonian Andes, 66, 996–1005, https://doi.org/10.1017/jog.2020.64, 2020.


%% Define the paths to be used 

% define the path to where IThikVel_v2_1 is located
path_to_folder = 'c:\users\lcsru\OneDrive\Documents\SPI\THICK_CALC_2021\'

% Output data folder
folder_output = [path_to_folder,'IThickVel_v2_1\OUTPUT_DATA\'];

% Input data folder
folder_input = [path_to_folder,'IThickVel_v2_1\INPUT_DATA\'];


%% Read input data

% Read slope grid
[input_grid,r,bbox] = geotiffread([folder_input, 'Perito_Moreno_sampling_data.tif'],1);
% No data values as NaN
input_grid(input_grid==-99999)=NaN;

% slope of glacier suface [rad]
SLOPE_GL = double(input_grid(:,:,1)); % First grid is slope
   
% glacier suface velocity [m/yr]
VEL_GL = double(input_grid(:,:,2));  % Second grid is surface velocity
  
 % glacier surface elevation [m]
 DEM_GL = double(input_grid(:,:,3)); % Third grid is surface velocity
 clear input_grid
 
 %% Ice thickness distribution calculation
 
 % First define the param
 %param = [ro,g,n,f,A,minslope,c,a];
   ro = 900; % ro = param(1,1); ice density [kgm-3]
   g = 9.8;% g = param(1,2);  gravity acceleration [ms-2]
   n = 3; % n = param(1,3);  n = 3 flow exponential
   f = 0.8; % f = param(1,4);  shape factor =0.8
   A =3.2e-24; % A = param(1,5);  %Pa–3 s–1 creep parameter
   minslope = 3; % minslope = param(1,6); % minimum threshold for slope (controls maximum
                % thickness)
   c = 1;% c=param(1,7);  if 1 glacier is calving if 0 no calving
   a = 0.25; % a = param(1,8);  proportion of Ub/Us from 0 to 1 (used if glacier is not
%                     calving
% Create the param vector
param = [ro,g,n,f,A,minslope,c,a];
 
% Use of  ithickvel_v2_1.m (only Ice thickness grid)
 [H]=ithickvel_v2_1(VEL_GL,SLOPE_GL,DEM_GL,param);
% Use of  ithickvel_v2_1.m (also Deformation and basal velocity)
 [H,Ud,Ub]=ithickvel_v2_1(VEL_GL,SLOPE_GL,DEM_GL,param);
  
 

%% Export grid results

% Ice thickness (H)
geotiffwrite([folder_output,'Perito_Moreno_example_v2_1_ice_thickness.tif'],H,r,'CoordRefSysCode',32718);

% Deformation velocity [m/yr]
% geotiffwrite([folder_output,'Perito_Moreno_example_v2_1_deformation_vel.tif'],Ud,r,'CoordRefSysCode',32718);
 
% Basal velocity [m/yr]   
% geotiffwrite([folder_output,'Perito_Moreno_example_v2_1_basal_vel.tif'],Ub,r,'CoordRefSysCode',32718); 

    %% Create some figures 
   % first define some variables
   [m,n]=size(H);
   x = r(3,1)+r(2,1).*(1:1:(n));
   y = r(3,2)+r(1,2).*(1:1:(m));
   clear r
   [XX_GL,YY_GL] = meshgrid(x,y);
  
   scrsz = get(0,'ScreenSize');
    figure('Position',[1 1 scrsz(3)-1 scrsz(4)-1])
    subplot(2,2,1) % Surface velocity map
    if isnan(floor(nanmin((nanmin(VEL_GL)))))~=1;
     stepV = (nanmax((nanmax(VEL_GL)))-nanmin((nanmin(VEL_GL))))/10;
     contourf(XX_GL,YY_GL,VEL_GL,'LineColor',[0 0 0],'LevelStep',max(floor(stepV),1),'Fill','on');
     contourcmap([floor(nanmin((nanmin(VEL_GL)))):max(floor(stepV),1):floor(nanmax((nanmax(VEL_GL))))],...
         'jet','colorbar','on','location','vertical')
     title(['Vel m a-1'])
    end 
    subplot(2,2,2) %Surface slope map
     stepV = (nanmax((nanmax(rad2deg(SLOPE_GL))))-nanmin((nanmin(rad2deg(SLOPE_GL)))))/10;
     contourf(XX_GL,YY_GL,rad2deg(SLOPE_GL),'LineColor',[0 0 0],'LevelStep',max(floor(stepV),1),'Fill','on');
     contourcmap([floor(nanmin((nanmin(rad2deg(SLOPE_GL))))):max(floor(stepV),1):floor(nanmax((nanmax(rad2deg(SLOPE_GL)))))],...
         'jet','colorbar','on','location','vertical')
         title(['Slope ° '])
     
    subplot(2,2,3) %Surface elevation map
     stepV = (nanmax((nanmax(DEM_GL)))-nanmin((nanmin(DEM_GL))))/10;
     contourf(XX_GL,YY_GL,DEM_GL,'LineColor',[0 0 0],'LevelStep',stepV,'Fill','on');
     contourcmap([floor(nanmin((nanmin(DEM_GL)))):floor(stepV):floor(nanmax((nanmax(DEM_GL))))],...
         'jet','colorbar','on','location','vertical')
          title(['Elevation m  GL ID '])
     
    subplot(2,2,4) % Ice thickness distribution map
     stepV = (nanmax((nanmax(H)))-nanmin((nanmin(H))))/10;
     contourf(XX_GL,YY_GL,H,'LineColor',[0 0 0],'LevelStep',max(floor(stepV),1),'Fill','on');
     contourcmap([floor(nanmin((nanmin(H)))):max(floor(stepV),1):floor(nanmax((nanmax(H))))],...
         'jet','colorbar','on','location','vertical')
          title(['Ice Thickness m '])
    clear stepV
    %Export Figure
    orient portrait 
    path_file_fig = [folder_output,'Perito_Moreno_example_v_2_figure_maps'];
    print(gcf,'-r300', '-dpsc', '-append', path_file_fig)
    close

