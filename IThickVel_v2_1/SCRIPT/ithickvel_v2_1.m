function [H,Ud,Ub]=ithickvel_v2_1(vel,slope,elevation,param)
%% ------------------------------------------------------------------------
% [H,Ud,Ub]=ithickvel_v2_1(vel,slope,elevation,param) [m]
% or
% [H]=ithickvel_v2_1(vel,slope,elevation,param) [m]
%
%% OUTPUT DATA
% H =ice thickness [m]
% Ub = basal velocity [m/a]
% Ud = deformation velocity [m/a]
%
%% INPUT DATA 
% vel = surface velocity [m/a]
% slope =  smooth slope [rad]
% surface = elevation of glacier surface [m]
%
%param = [ro,g,n,f,A,minslope,c,a];
   % ro = param(1,1); ice density [kgm-3]
   % g = param(1,2);  gravity acceleration [ms-2]
   % n = param(1,3);  n = 3 flow exponential
   % f = param(1,4);  shape factor =0.8
   % A = param(1,5);  %Pa–3 s–1 creep parameter
   % minslope = param(1,6); % minimum threshold for slope (controls maximum
                % thickness)
   % calving=param(1,7);  if 1 glacier is calving if 0 no calving
   % a = param(1,8);  proportion of Ub/Us from 0 to 1 (used if glacier is not
%                     calving
%
%% ICE THICKNESS CALC METHOD 
% It is based on the method discribed by:
% Gantayat P, Kulkarni AV and Srinivasan J (2014)
% Estimation of ice thickness using surface velocities and slope: case study
% at Gangotri Glacier, India. Journal of Glaciology 60(220), 277–282
%(doi:10.3189/2014JoG13J078)
% Calculate ice thickness from the equation of laminar flow. 
%
%Us=Ud+Ub
%Ud=((2*A)/n+1).*taob.^n.*H;
%taob=f*ro*g.*H.*sin(slope);
%Ub=0.25.*Us (they use this simplification I could calculated).
%H=(1.5.*Us./(A.*(f*ro*g.*sin(slope)).^3))^(1/4)
% In order to use a better solution of Ub I propose to use amore general
% form of your equation:
%
%H=((Us-Ub)*(n+1)/(2*A)*1./(f*ro*g.*sin(slope))^n)^(1/(n+1))
% 
% Lucas Ruiz (lruiz@mendoza-conicet.gob.ar)
%


%% Definition of variables from param

 ro = param(1,1);
 g = param(1,2);
 n = param(1,3);
 f = param(1,4);
 A = param(1,5);
 minslope = param(1,6);
 calving=param(1,7);
 a = param(1,8);

%% grid size
 [xsize,ysize] = size(elevation);

%% Slope input
alfa = slope; %slope
% correction of minimum slope
%in case minimum slope is lower than threshold apply correction
if min(min(slope))<deg2rad(minslope);
    display('Minimum slope threshold applied')
    display(['Min slope =',num2str(rad2deg(min(min(slope)))),'°',' Threshold = ',num2str((minslope)),'°'])
    alfa(alfa<=deg2rad(minslope))=deg2rad(minslope); %apply minimum treshold.
end

%% surface ice velocity [m a-1]
 
Us = vel./31536000; 

%% subglacial velocity [m a-1]

if calving==0;
   Ub=a.*Us; % Ub if there is no calving
   
elseif calving==1; % if there is calving
    aa =zeros(size(vel)); %prelocation of ub/us param
    znorm =(nanmax(nanmax(elevation))-elevation)./(nanmax(nanmax(elevation))-nanmin(nanmin(elevation))); %normalization of Elevation data
    
    % Relation between ub/us and znorm from Stuefer  (1999) estimations for Perito Moreno
    %Stuefer, M. 1999. Investigations on mass balance and dynamics of
    %Moreno Glacier based on field measurements and satellite
    %imagery. (PhD thesis, University of Innsbruck.)
     
    % a depends on how close to the front or to the minum elevation it is
    table_norm = [0,0.25;0.4,0.25;0.6,0.26;0.7,0.28;0.8,0.32;0.9,0.4;0.91,0.5;0.92,0.64;0.937,0.83;0.97,0.97;0.98,0.996;1,0.997];
    
    % creat anorm and aa
    x = 0:0.01:1;
    y = interp1(table_norm(:,1),table_norm(:,2),x,'cubic','extrap');
    anorm = smooth(y,5);
    for j = 1:1:xsize;
        for k = 1:1:ysize;
            for l = 1:1:length(x)-1;
                if znorm(j,k)>=x(l) && znorm(j,k)<x(l+1);
                    aa(j,k) =y(l); %Ub/us for the elevation
                else if znorm(j,k)>=x(l+1)
                    aa(j,k) = y(end);
                end
            end
        end
    end

  Ub = aa.*Us;  % Ub if there is calving
    end
end

%% Ice thickness calc

H = real(((((Us-Ub).*(n+1))/(2*A))*1./(f*ro*g.*sin(alfa)).^n).^(1/(n+1))); %ice thickness [m]

%% Zero ice thickness at the contour of the glacier (un-comment if you want to use it)
%Glacier contours have H = zero)
%H(isnan(gl))=0;

%% Simple filter to smooth the final grid (un-comment if you want to use it)
%K = 0.125*ones(3); %define 3x3 kernel
%Hsmooth = conv2(H,K,'same'); %smoothed thickness [m]

%H=Hsmooth;

%% Deformation velocity and Subglacial velocity calculation
if nargout >  1
    % Deformation velocity
     taob=f*ro*g.*H.*sin(alfa);
    Ud=((2*A)/(n+1)).*H.*taob.^3;
    Ud=Ud.*31536000;
    % Subglacial velocity
     Ub=Ub.*31536000;
end

end

