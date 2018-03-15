%%%%Global Gravity change due to earthquake impact using Okubo 1992%%%%%%%%
%This code allow to convert mass distribution, deformation and change in elevation triggered by earthquake to a change of gravity.
%Tape "globalcmt_create" to download the catalog
%Tape "globalcmt_update" in command windows to update last earthquake from CMT
%% Clean up
clc; close all; clear all;
%% Earthquake catalogue
cd('C:\Users\mikpo\My Cloud\Stage_2017_Michael_Pons\seizmo\event')
load('globalcmt_full.mat');
%% Rlistfunction
% Input model parameters
OPEN=0;     % OPEN   : dislocation in tensile component [m]
RHO=2670;   % RHO    : density of the medium [kg/m^3]
RHOP=2670;  % RHOP   : optional density of the cavity-filling matter [kg/m^3]; takes RHO value if not specified.
mu=30e9;    % Shear modulus (Pa) for Leonard2010 model
BETA=0; %No gravity  air correction
NU=0.25;%Poisson coef
coef_elag=5;
coef_tile=90;
%% Seismic moment convertion to N.m
Moexpo=scalarmoment.*power(10,exponent)*1e-7;
Mw_catalog = (2/3)*log10(Moexpo*1e7)-10.7;
%% Define your area
prompt = {'longlimwest','lonlimest','latlimdown','latlimup'};
dlg_title = 'Define Area';
num_lines = 1;
defaultans ={'-180','180','-90','90'};%{'60','120','-25','25'}%{'120','122','21.5','26'};%{'-180','180','-90','90'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
lonlimwest =str2num(answer{1});
lonlimest  =str2num(answer{2});
latlimdown =str2num(answer{3});
latlimup   =str2num(answer{4});
%% Define date
dateRquake=datetime(year,month,day);
%datelist=datetime(year,month,day);
prompt4={'Starting Year','Starting Month','Starting Day','Ending Year','Ending Month','Ending Day'};
dlg_title= 'Define date';
num_lines=1;
defaultans={'2002','04','17','2014','10','11'};%{'2004','12','24','2004','12','26'};%{'2002','1','1','2012','1','1'};
answer4=inputdlg(prompt4,dlg_title,num_lines,defaultans);
year1=str2num(answer4{1});month1=str2num(answer4{2});day1=str2num(answer4{3});
year2=str2num(answer4{4});month2=str2num(answer4{5});day2=str2num(answer4{6});
datein=datetime(year1,month1,day1);dateout=datetime(year2,month2,day2);
if datein<=min(dateRquake)
    datein=min(dateRquake);
elseif datein>max(dateRquake)
    disp('No earthquake found');
end
if dateout<min(dateRquake)
    disp('No earthquake found');
elseif dateout>max(dateRquake)
    dateout=max(dateRquake);
end
%%Look for earthquakes defined in Rlistfunction and find smaller Magnitude
k=find(longitude>=lonlimwest & longitude<=lonlimest & dateRquake>=datein & dateRquake<=dateout);
%find all the earthquake in longitude limit, latitude limit and Mo limit and Give the index to keep for the concerned zone
for T=k
    ind1=find(latitude(T)>=latlimdown & latitude(T)<=latlimup & dateRquake(T)>=datein & dateRquake(T)<=dateout);
end
if isempty(ind1)==1
    disp('No earthquake found in this zone at this time');
elseif isempty(ind1)==0
    disp([num2str(numel(ind1)),' earthquakes found in this zone at this time']);
end
Rt_all=k(ind1)';
Mo_all=Moexpo(Rt_all);
Mw_all = (2/3)*log10(Mo_all*1e7)-10.7; % do the conversion of each Mo of each earthquake in Mw.
minMw=min(Mw_all);
maxMw=max(Mw_all);
%% Magnitude range
    M(1,:)=[minMw maxMw];
%% Found all earthquakes for defined ranges
for Mnum=1:size(M,1)
    disp([num2str(Mnum),'/',num2str(size(M,1)),' Magnitude range defined'])
    Mwspace=M(Mnum,:);
    Mospacemin=10^(1.5*Mwspace(1)+9.105);
    Mospacemax=10^(1.5*Mwspace(end)+9.105);
    %% Find all earthquakes
    k=find(longitude>=lonlimwest & longitude<=lonlimest & dateRquake>=datein & dateRquake<=dateout);
    %find all the earthquake in longitude limit, latitude limit and Mo limit and Give the index to keep for the concerned zone
    %date filter
    for T=k
        ind1=find(latitude(T)>=latlimdown & latitude(T)<=latlimup & Moexpo(T)>=Mospacemin & Moexpo(T)<=Mospacemax & Moexpo(T)~=0 & dateRquake(T)>=datein & dateRquake(T)<=dateout);
    end
    if isempty(ind1)==1
        disp('No earthquake found in this zone at this time for this magnitude range');
    elseif isempty(ind1)==0
        disp([num2str(numel(ind1)),' earthquakes found in this zone at this time for this magnitude range ']);
    end
    %% Extracted parameters for earthquakes define
    Rt=k(ind1)';
    Earthquake_number(Mnum)=numel(ind1)%Earthquake number found for the range of magnitude
    Yo=centroidlat(Rt);Xo=centroidlon(Rt);dip1a=dip1(Rt);dip2a=dip2(Rt);strike1a=strike1(Rt);strike2a=strike2(Rt);rake1a=rake1(Rt);rake2a=rake2(Rt);depth1a=centroiddep(Rt).*1e3;mb1a=mb(Rt);Mo=Moexpo(Rt);Rname=name(Rt);%depth(Rt)
    for i=1:numel(ind1)
        [X(i),Y(i),zone(i)]=ll2utm(Yo(i),Xo(i)); % Convert from decimal degree to meters
    end
    Mw = (2/3)*log10(Mo*1e7)-10.7;

 
end 

%% Plot and projection
utmlonlim=[lonlimwest lonlimest];
utmlatlim=[latlimdown latlimup];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currDir=cd;
cd('C:\Users\mikpo\Desktop\StageM2Rennes\Programmes\Okubo\Okubo2\Results\Tendancecum12ansL60ddk3')
slope=importdata('slopeb.mat');
lat1=importdata('latb.mat');
lon1=importdata('lonb.mat');
cd(currDir)

markerWidth = (Mw-3).^3;
 m_proj('Miller Cylindrical','long',[-179.94 179.94],'lat',[-70 80]); 
%m_proj('Miller Cylindrical','long',utmlonlim,'lat',utmlatlim);
set(figure,'Position',[1 1 2000 1500]) % [coin_x coin_y hauteur(px) largeur(px)]
set(gcf,'PaperPositionMode','auto');hold on
m_pcolor(lon1,lat1,slope),shading interp;set(gca, 'CLim',[-30, 30]);
%colorbar;h = colorbar;ylabel(h, '\muGal','FontSize',12);colormap(jet);hold off;

hold on;
Sc = m_scatter(Xo,Yo,markerWidth,'k','filled');
currentunits = get(gca,'Units');
set(gca, 'Units', 'Points');
axpos = get(gca,'Position');
set(gca, 'Units', currentunits);
set(Sc, 'SizeData', markerWidth);
m_coast('color','k','LineWidth',1);
m_grid('box','fancy','tickdir','out');
set(gca,'FontSize',20);
colorbar;h = colorbar;ylabel(h, '\muGal','FontSize',20);colormap(jet)




