%%%%Global Gravity change due to earthquake impact using Okubo 1992%%%%%%%%

%This code can calculate change of gravity or elevation using CMT catalog and analytic solution of Okubo 1992
%Download Seizmo then
%Tape "globalcmt_create" to download the catalog
%Tape "globalcmt_update" in command windows to update last earthquake from CMT

%% Clean up
clc; close all; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parametrisation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Earthquake catalogue
load('globalcmt_full.mat');

%% list for rupture model known
Sumatraname='M122604A';Tohokuname='M201103110546A';Chile10name='C201002270634A';Iquique14name='C201404012346A';
Balochistan13name='C201309241129A';Carlsberg03name='C071503F';Bhuj01name='M012601A';Jaya04name='C020704A';Ca04name='C092804H';
Haiti10name='C201001122153A';EMC10name='C201004042240A';S_Sumatra07name='C200709121110A';Kash13name='C201304161044A';
C_Sumatra07name='C200709122348A';Kuril07name='C200701130423A';Kuril06name='C200612151659A';Mexico12name='C201203201802A';
Padang09name='C200909301016A';Peru07name='C200708152340A';Sichuan08name='C200805120628A';Simeulue08name='C200802200808A';
Sulawesi08name='C200811161702A';Sumatra05name='C200503281609A';Tocopilla07name='C200711141540A';Vanuatu09name='C200910072203A';

%% Input model parameters
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

%% Choice for calcultating dG and/or dH
options.Interpreter = 'tex';
% Include the desired Default answer
options.Default = 'dG';
% Use the TeX interpreter in the question
qstring = 'Calculate dG or dH ?';
buttonvar = questdlg(qstring,'Boundary Condition',...
    'dG','dH','Both',options);
if strcmp(buttonvar,'dG')
    outdata=0;
elseif strcmp(buttonvar,'dH')
     outdata=1;
     elseif strcmp(buttonvar,'Both')
     outdata=2;
end

%% Free Air correction or not
options.Interpreter = 'tex';
% Include the desired Default answer
options.Default = 'Yes';
% Use the TeX interpreter in the question
qstring = 'Apply free air correction ?';
button = questdlg(qstring,'Boundary Condition',...
    'Yes','No',options);

%% Define your area
prompt = {'longlimwest','lonlimest','latlimdown','latlimup'};
dlg_title = 'Define Area';
num_lines = 1;
defaultans ={'130','150','25','50'};%{'60','120','-25','25'}%{'120','122','21.5','26'};%{'-180','180','-90','90'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
lonlimwest =str2num(answer{1});
lonlimest  =str2num(answer{2});
latlimdown =str2num(answer{3});
latlimup   =str2num(answer{4});

%% Define model
prompt2 = {'Model : 1-Leo10a, 2-Leo10b, 3-Leo10c, 4-YenMa11a, 5-YenMa11b, 6-YenMa11c, 7-WC94a, 8-WC94b '};
dlg_title = 'Define Model';
num_lines = 1;
defaultans2 = {'5'};
answer2 = inputdlg(prompt2,dlg_title,num_lines,defaultans2);
Model=str2num(answer2{1});

%% Nodal plan
prompt4 = inputdlg('Define nodal plan 1 or 2 or 3 if mixt','Nodal Plan choice:',[1 40]);
nodal = str2num(prompt4{:});
if any(nodal==3)
    options.Interpreter = 'tex';
    % Include the desired Default answer
    options.Default = 'Max';
    % Use the TeX interpreter in the question
    qstring1 = 'Calcul max(abs(dG)) for each earthquake or Min(abs(dG)) ?';
    button1 = questdlg(qstring1,'Boundary Condition',...
        'Max','Min',options);
end

%% Define date
dateRquake=datetime(year,month,day);
%datelist=datetime(year,month,day);
prompt4={'Starting Year','Starting Month','Starting Day','Ending Year','Ending Month','Ending Day'};
dlg_title= 'Define date';
num_lines=1;
defaultans={'2011','3','10','2011','3','12'}%{'2004','12','24','2004','12','26'};%{'2002','1','1','2012','1','1'};
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

%% Look for earthquakes defined in Rlistfunction and find smaller Magnitude
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
prompt3a={'Preconfig range of magnitude ','or specify space-separated magnitude if not leave empty'};
dlg_title='Define Range magnitude';
num_lines=1;
defaultans3a={'Preconfig',num2str(minMw)};
answer3a=inputdlg(prompt3a,dlg_title,num_lines,defaultans3a);
choice2=str2num(answer3a{2});
if isempty(choice2==1)
    M(1,:)=[minMw 4.999];
    M(2,:)=[5 5.999];
    M(3,:)=[6 6.999];
    M(5,:)=[7 7.999];
    M(6,:)=[8 8.999];
    M(7,:)=[9 10];
else
    M(1,:)=str2num(answer3a{2});
end

%% Check for earthquake in Rlistfunction
name2=name(Rt_all);
for Rlist=1:numel(name2)
    if strcmp(name2(Rlist),Sumatraname)==1||strcmp(name2(Rlist),Tohokuname)==1||strcmp(name2(Rlist),Chile10name)==1 ||strcmp(name2(Rlist),Iquique14name)==1 ...
            ||strcmp(name2(Rlist),Balochistan13name)==1||strcmp(name2(Rlist),Carlsberg03name)==1||strcmp(name2(Rlist),Bhuj01name)==1||strcmp(name2(Rlist),Jaya04name)==1 ...
            ||strcmp(name2(Rlist),Ca04name)==1||strcmp(name2(Rlist),Haiti10name)==1||strcmp(name2(Rlist),EMC10name)==1||strcmp(name2(Rlist),S_Sumatra07name)==1 ...
            ||strcmp(name2(Rlist),Kash13name)==1||strcmp(name2(Rlist),C_Sumatra07name)==1||strcmp(name2(Rlist),Kuril07name)==1||strcmp(name2(Rlist),Kuril06name)==1 ...
             ||strcmp(name2(Rlist),Mexico12name)==1||strcmp(name2(Rlist),Padang09name)==1||strcmp(name2(Rlist),Peru07name)==1||strcmp(name2(Rlist),Sichuan08name)==1 ...
             ||strcmp(name2(Rlist),Simeulue08name)==1||strcmp(name2(Rlist),Sulawesi08name)==1||strcmp(name2(Rlist),Sumatra05name)==1||strcmp(name2(Rlist),Tocopilla07name)==1||strcmp(name2(Rlist),Vanuatu09name)==1;
        disp('fault geometry proposed'); RR(Rlist)=1;
    else
        RR(Rlist)=0;
    end
end
if sum(RR)~=0
        str={'Select earthquake(s)parameters or cancel'};
        S={'Sumatra2004';'Tohoku2011';'Chile2010';'Iquique2014';'Balochistan2013';...
            'Carlsberg2003';'Bhuj2001';'Jaya2004';'Ca2004';'Haiti2010';'EMC2010';'S_Sumatra2007';'Kash2013';...
            'C_Sumatra2007';'Kuril2007';'Kuril2006';'Mexico2012';'Padang2009';'Peru2007';'Sichuan2008';...
            'Simeulue2008';'Sulawesi2008';'Sumatra2005';'Tocopilla2007';'Vanuatu2009'};
        result=listdlg('PromptString',str,'ListSize',[200 100],'ListString',S,'SelectionMode','multiple');
elseif sum(RR)==0
    result=[];
end

%% Coefficient to lower the resolution
prompt1a={'Lower the Optimal resolution','or choose a dec deg resolution. If not leave empty'};
dlg_title='Resolution';
num_lines=1;
defaultans1a={'10',''};
answer1a=inputdlg(prompt1a,dlg_title,num_lines,defaultans1a);
coef_resol=str2num(answer1a{1});
defdx=str2num(answer1a{2});

%% Define path for results folder
dname = uigetdir('C:\Users\mikpo\Desktop\StageM2Rennes\Programmes\Okubo\Okubo2\Results');

%% Parameters to define dx, x
dx_bestDSREV=importdata('dx_bestDSREV.mat');EDGDSREV=importdata('EDGmatDSREV.mat');
dx_bestDSNORM=importdata('dx_bestDSNORM.mat');EDGDSNORM=importdata('EDGmatDSNORM.mat');
dx_bestSS=importdata('dx_bestSS.mat');EDGSS=importdata('EDGmatSS.mat');
Mwint=0:0.25:9 ;%used later for the interpolation of dx depending on the Magnitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% defining the resolution for each earthquake calculation%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    %% Calcul of rupture geometry
    for R=1:numel(ind1)
        disp([num2str(R),'/',num2str(numel(ind1)),' Geometry'])
        Mwvec(R)=Mw(R) ;

        %% Scaling Laws
        [D,Z,L,W]=loi(Mwvec(R),rake1a(R),mu,dip1a(R),Mo(R),Model);
        %% Define the grid size and resolution for each earthquake
        % Change dataset depending on the Rake
        if Mw(R)<=9
            if (rake1a(R) <=135 && rake1a(R) >=45)% Dip-slip |reverse fault
                EDGmat=EDGDSREV;disp ('la1');dx_best=dx_bestDSREV;
            elseif (rake1a(R) <=-45 && rake1a(R) >=-135)% Dip-slip |normal fault
                EDGmat=EDGDSNORM;disp ('la2');dx_best=dx_bestDSNORM;
            elseif (rake1a(R) >-45 && rake1a(R)<=45)%Strike slip
                EDGmat=EDGSS;disp ('la3');dx_best=dx_bestSS;
            elseif (rake1a(R)>135 && rake1a(R)<=180)%Strike slip
                EDGmat=EDGSS;disp ('la4');dx_best=dx_bestSS;
            elseif (rake1a(R) <-135 && rake1a(R)>=-180)%Strike slip
                EDGmat=EDGSS;disp ('la5');dx_best=dx_bestSS;
            end
            %% interpolation of x and dx
            for j=Model
                xint1(R)=interp1(Mwint,EDGmat(j,:),Mwvec(R),'linear');
                dxint(R)=interp1(Mwint,dx_best(j,:),Mwvec(R),'linear');
                xx(R)=xint1(R)<L(j);
                if isempty (xx(R))==1
                    xint(R)=xint1(R)+L(j)*coef_elag;
                    xmin(R)=-round(xint(R));
                    xmax(R)=round(xint(R));
                elseif isempty (xx(R))==0
                    xint(R)=xint1(R)+L(j)*coef_elag;
                    xmin(R)=-round(xint(R));
                    xmax(R)=round(xint(R));
                end
            end
        elseif Mw(R)>9
            if (rake1a(R) <=135 && rake1a(R) >=45)% Dip-slip |reverse fault
                EDGmat=EDGDSREV(:,end);disp ('la1');dx_best=dx_bestDSREV(:,end);
            elseif (rake1a(R) <=-45 && rake1a(R) >=-135)% Dip-slip |normal fault
                EDGmat=EDGDSNORM(:,end);disp ('la2');dx_best=dx_bestDSNORM(:,end);
            elseif (rake1a(R) >-45 && rake1a(R)<=45)%Strike slip
                EDGmat=EDGSS(:,end);disp ('la3');dx_best=dx_bestSS(:,end);
            elseif (rake1a(R)>135 && rake1a(R)<=180)%Strike slip
                EDGmat=EDGSS(:,end);disp ('la4');dx_best=dx_bestSS(:,end);
            elseif (rake1a(R) <-135 && rake1a(R)>=-180)%Strike slip
                EDGmat=EDGSS(:,end);disp ('la5');dx_best=dx_bestSS(:,end);
            end
            for j=Model
                dxint(R)=dx_best(j,end);
                xmin(R)=-EDGmat(j,end)-L(j)*coef_elag;
                xmax(R)=EDGmat(j,end)+L(j)*coef_elag;
            end
        end
        if isempty(defdx)==1
            %% Define Est North-vector for Okubo
            E=xmin(R):dxint(R):xmax(R);
            N=xmin(R):dxint(R):xmax(R);
            %Georeferencement
            lon1=E+X(R);
            lat1=N+Y(R);
            %Re-convert from meters to decimal degree using utm2ll fonction
            [lon2,lat2]=meshgrid(lon1,lat1);
            %[latg,long]=utm2ll(lon1,lat1,zone(R));
            [latg,long] = utm2ll(lon2(1:numel(lon2)),lat2(1:numel(lat2)),zone(R));
            long = wrapTo180(long);
            %Calcul of lat-long limits for each earthquakes and taking account of all earthquakes
            %for each Earthquakes
            minlat(R)=min(latg);
            minlong(R)=min(long);
            maxlat(R)=max(latg);
            maxlong(R)=max(long);
        end
    end

    %% Define the smallest the step of calculation dx (in decimal degree)
    if isempty(defdx)==1
        %For all Earthquakes
        minlat2=min(minlat);
        minlong2=min(minlong);
        maxlat2=max(maxlat);
        maxlong2=max(maxlong);
        mindx=find(dxint==min(dxint),1);
        for T = mindx
            Emindx=xmin(T)+dxint(T):dxint(T):xmax(T)+dxint(T);
            Nmindx=xmin(T)+dxint(T):dxint(T):xmax(T)+dxint(T);
            lon1mindx=Emindx+X(T);
            lat1mindx=Nmindx+Y(T);
            %Re-convert from meters to decimal degree using utm2ll fonction
            [latg2,long2]=utm2ll(lon1mindx,lat1mindx,zone(T));
        end
        dxgrid=min(diff(latg2))*coef_resol; %decimal degree
        % Building the grid for all earthquakes
        longT{Mnum}=minlong2:dxgrid:maxlong2;
        latT{Mnum}=minlat2:dxgrid:maxlat2;
    else
        dxgrid1=defdx;% in decimal degree

        %%  Building the grid for all earthquakes
        longT{Mnum}=-180:dxgrid1:180;
        latT{Mnum}=-90:dxgrid1:90;
    end
    [longT2,latT2]=meshgrid(longT{Mnum},latT{Mnum});
    % construction of the total grid where all the local grid of each earthuake will be interpolated
    matot=zeros(numel(latT{Mnum}),numel(longT{Mnum}));
    if outdata==2;
        matot2=zeros(numel(latT{Mnum}),numel(longT{Mnum}));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Calculation of gravity change %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Okubo calculation
    for R=1:numel(Rt)
        Mwvec(R)=Mw(R) ;
        %% Rupture model known
        if isempty(defdx)==0
            dxgrid=dxgrid1;
        elseif isempty(defdx)==1
            dxgrid1=dxgrid;
        end
        if strcmp(Rname(R),Sumatraname)==1 && any(result==1)%Sumatra 2004 Banerjee
            disp('Sumatra function');[Gq,Gq2]=Sumatra(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Chile10name)==1 && any(result==3)
            disp('Chile function');[Gq,Gq2]=Chile2010(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Iquique14name)==1 && any(result==4)
            disp('Iquique function');[Gq,Gq2]=Iquique2014(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Balochistan13name)==1 && any(result==5)
            disp('Balochistan function');[Gq,Gq2]=Balochistan2013(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Carlsberg03name)==1 && any(result==5)
            disp('Carlsberg2004 function');[Gq,Gq2]=Carlsberg2004(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Bhuj01name)==1 && any(result==6)
            disp('Bhuj2001 function');[Gq,Gq2]=Bhuj2001(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Jaya04name)==1 && any(result==7)
            disp('Jaya2004 function');[Gq,Gq2]=Jaya2004(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Ca04name)==1 && any(result==8)
            disp('Ca2004 function');[Gq,Gq2]=Ca2004(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Haiti10name)==1 && any(result==9)
            disp('Haiti2010 function');[Gq,Gq2]=Haiti2010(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),EMC10name)==1 && any(result==10)
            disp('EMC2010 function');[Gq,Gq2]=EMC2010(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),S_Sumatra07name)==1 && any(result==11)
            disp('S_Sumatra2007 function');[Gq,Gq2]=S_Sumatra2007(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Kash13name)==1 && any(result==12)
            disp('Kash2013 function');[Gq,Gq2]=Kash2013(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),C_Sumatra07name)==1 && any(result==13)
            disp('CSumatra2007 function');[Gq,Gq2]=CSumatra2007(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Kuril07name)==1 && any(result==14)
            disp('Kuril2007 function');[Gq,Gq2]=Kuril2007(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Kuril06name)==1 && any(result==15)
            disp('Kuril2006 function');[Gq,Gq2]=Kuril2006(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Mexico12name)==1 && any(result==16)
            disp('Mexico2012 function');[Gq,Gq2]=Mexico2012(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Padang09name)==1 && any(result==17)
            disp('Padang2009 function');[Gq,Gq2]=Padang2009(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Peru07name)==1 && any(result==18)
            disp('Peru2007 function');[Gq,Gq2]=Peru2007(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Sichuan08name)==1 && any(result==19)
            disp('Sichuan2008 function');[Gq,Gq2]=Sichuan2008(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Simeulue08name)==1 && any(result==20)
            disp('Simeulue2008 function');[Gq,Gq2]=Simeulue2008(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Sulawesi08name)==1 && any(result==21)
            disp('Sulawesi2008 function');[Gq,Gq2]=Sulawesi2008(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Sumatra05name)==1 && any(result==22)
            disp('Sumatra2005 function');[Gq,Gq2]=Sumatra2005(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Tocopilla07name)==1 && any(result==23)
            disp('Tocopilla2007 function');[Gq,Gq2]=Tocopilla2007(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Vanuatu09name)==1 && any(result==24)
            disp('Vanuatu2009 function');[Gq,Gq2]=Vanuatu2009(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        elseif  strcmp(Rname(R),Tohokuname)==1 && any(result==2)
            disp('Tohoku function');[Gq,Gq2]=Tohoku(Model,dxgrid,longT2,latT2,coef_tile,button,outdata);
        else

            %% NormalCase_Scaling Laws
            [D,Z,L,W]=loi(Mwvec(R),rake1a(R),mu,dip1a(R),Mo(R),Model);
            disp([num2str(R),'/',num2str(numel(ind1)),' Okubo calculation'])
            %% Define Est North-vector for Okubo
            E=xmin(R):dxint(R):xmax(R);
            N=xmin(R):dxint(R):xmax(R);
            [E1,N1]=meshgrid(E,N);
            % Georeferencement of the earthquake
            lon1=E+X(R);
            lat1=N+Y(R);
            % Okubo function
            tic

            %% Calculation depending on nodal 1 2 or both maximum or minim signal dG
            % For Nodal 1
            if any(nodal==1)||isempty(nodal==0)
                disp([num2str(R),'/',num2str(numel(ind1)),' Okubo calculation for Nodal plan 1'])
                if strcmp(button,'No')
                    [dG,dH]=okubo92(E1,N1,depth1a(R),strike1a(R),dip1a(R),L(j),W(j),rake1a(R),D(j),OPEN,RHO,RHOP,BETA,NU);
                    S(R)=1e8.*max(max(abs(dG)));
                else
                    [dG,dH]=okubo92(E1,N1,depth1a(R),strike1a(R),dip1a(R),L(j),W(j),rake1a(R),D(j),OPEN,RHO,RHOP);
                    S(R)=1e8.*max(max(abs(dG)));
                end
                 if outdata==1
                    dvar1=dH;
                elseif outdata==0
                    dvar1=dG;
                elseif outdata==2
                    dvar1=dG;
                    dvar2=dH;
                end
                if Mw(R)<=6.5 %interpolation option 1 is faster when distorsion is minimized
                    [latg,long]=utm2ll(lon1,lat1,zone(R));
                    long = wrapTo180(long);
                    [longR1,latgR1]=meshgrid(long,latg);
                    Gq=interp2(longR1,latgR1,dvar1,longT2,latT2,'linear');
                    Gq(isnan(Gq)==1)=0;
                    if outdata==2
                        Gq2=interp2(longR1,latgR1,dvar2,longT2,latT2,'linear');
                        Gq2(isnan(Gq2)==1)=0;
                    end
                else
                    %Re-convert from meters to decimal degree using utm2ll fonction
                    [lon2,lat2]=meshgrid(lon1,lat1);
                    [LA2,LO2] = utm2ll(lon2(1:numel(lon2)),lat2(1:numel(lat2)),zone(R));
                    indlocallon=find(wrapTo360(longT2(1,:))>=wrapTo360(min(LO2)) & wrapTo360(longT2(1,:))<=wrapTo360(max(LO2)));
                    indlocallat=find(latT2(:,1)>=min(LA2) & latT2(:,1)<=max(LA2));
                    LO2 = wrapTo180(LO2);
                    Gqlocal = gridfit(LO2,LA2,dvar1(1:numel(dvar1)),longT2(1,indlocallon),latT2(indlocallat,1),'tilesize',round(numel(longT{Mnum})/coef_tile),'overlap',0.25);
                    Gqlocal(isnan(Gqlocal)==1)=0;
                    Gq=zeros(size(longT2));
                    Gq(indlocallat,indlocallon)=Gq(indlocallat,indlocallon)+Gqlocal;
                    if outdata==2
                        Gqlocal2 = gridfit(LO2,LA2,dvar2(1:numel(dvar2)),longT2(1,indlocallon),latT2(indlocallat,1),'tilesize',round(numel(longT{Mnum})/coef_tile),'overlap',0.25);
                        Gqlocal2(isnan(Gqlocal)==1)=0;
                        Gq2=zeros(size(longT2));
                        Gq2(indlocallat,indlocallon)=Gq2(indlocallat,indlocallon)+Gqlocal2;
                    end
                end
                % For Nodal 2
            elseif any(nodal==2)%nodal 2
                disp([num2str(R),'/',num2str(numel(ind1)),' Okubo calculation for Nodal plan 2'])
                if strcmp(button,'No')
                    [dG,dH]=okubo92(E1,N1,depth1a(R),strike2a(R),dip2a(R),L(j),W(j),rake2a(R),D(j),OPEN,RHO,RHOP,BETA,NU);
                    S(R)=1e8.*max(max(abs(dG)));
                else
                    [dG,dH]=okubo92(E1,N1,depth1a(R),strike2a(R),dip2a(R),L(j),W(j),rake2a(R),D(j),OPEN,RHO,RHOP);
                    S(R)=1e8.*max(max(abs(dG)));
                end
                if outdata==1
                    dvar1=dH;
                elseif outdata==0
                    dvar1=dG;
                elseif outdata==2
                    dvar1=dG;
                    dvar2=dH;
                end
                if Mw(R)<=6.5 %interpolation option 1 is faster when distorsion is minimized
                    [latg,long]=utm2ll(lon1,lat1,zone(R));
                    long = wrapTo180(long);
                    [longR1,latgR1]=meshgrid(long,latg);
                    Gq=interp2(longR1,latgR1,dvar1,longT2,latT2,'linear');
                    Gq(isnan(Gq)==1)=0;
                    if outdata==2
                    Gq2=interp2(longR1,latgR1,dvar2,longT2,latT2,'linear');
                    Gq2(isnan(Gq2)==1)=0;
                    end
                else
                    %Re-convert from meters to decimal degree using utm2ll fonction
                    [lon2,lat2]=meshgrid(lon1,lat1);
                    [LA2,LO2] = utm2ll(lon2(1:numel(lon2)),lat2(1:numel(lat2)),zone(R));
                    indlocallon=find(wrapTo360(longT2(1,:))>=wrapTo360(min(LO2)) & wrapTo360(longT2(1,:))<=wrapTo360(max(LO2)));
                    indlocallat=find(latT2(:,1)>=min(LA2) & latT2(:,1)<=max(LA2));
                    LO2 = wrapTo180(LO2);
                    Gqlocal = gridfit(LO2,LA2,dvar1(1:numel(dvar1)),longT2(1,indlocallon),latT2(indlocallat,1),'tilesize',round(numel(longT{Mnum})/coef_tile),'overlap',0.25);
                    Gqlocal(isnan(Gqlocal)==1)=0;
                    Gq=zeros(size(longT2));
                    Gq(indlocallat,indlocallon)=Gq(indlocallat,indlocallon)+Gqlocal;
                    if outdata==2
                    Gqlocal2 = gridfit(LO2,LA2,dvar2(1:numel(dvar2)),longT2(1,indlocallon),latT2(indlocallat,1),'tilesize',round(numel(longT{Mnum})/coef_tile),'overlap',0.25);
                    Gqlocal2(isnan(Gqlocal2)==1)=0;
                    Gq2=zeros(size(longT2));
                    Gq2(indlocallat,indlocallon)=Gq2(indlocallat,indlocallon)+Gqlocal2;
                    end
                end
                % For Nodal 3
            elseif any(nodal==3)%nodal mixte
                disp([num2str(R),'/',num2str(numel(ind1)),' Okubo calculation for mixte Nodal plan 1 and 2'])
                if strcmp(button,'No')
                    [dG1,dH1]=okubo92(E1,N1,depth1a(R),strike1a(R),dip1a(R),L(j),W(j),rake1a(R),D(j),OPEN,RHO,RHOP,BETA,NU);
                    S2(R)=1e8.*max(max(abs(dG1)));%look at the max
                    [dG2,dH2]=okubo92(E1,N1,depth1a(R),strike2a(R),dip2a(R),L(j),W(j),rake2a(R),D(j),OPEN,RHO,RHOP,BETA,NU);
                    S3(R)=1e8.*max(max(abs(dG2)));%look at the max
                else
                    [dG1,dH1]=okubo92(E1,N1,depth1a(R),strike1a(R),dip1a(R),L(j),W(j),rake1a(R),D(j),OPEN,RHO,RHOP);
                    S2(R)=1e8.*max(max(abs(dG1)));%look at the max
                    [dG2,dH2]=okubo92(E1,N1,depth1a(R),strike2a(R),dip2a(R),L(j),W(j),rake2a(R),D(j),OPEN,RHO,RHOP);
                    S3(R)=1e8.*max(max(abs(dG2)));%look at the max
                end
                if strcmp(button1,'Max')
                    Sa(R)=S2(R);
                    Sb(R)=S3(R);
                else
                    Sa(R)=S3(R);
                    Sb(R)=S2(R);
                end
                if Sa(R)>=Sb(R) % Choose the max
                    if outdata==1
                        dvar1=dH1;
                    elseif outdata==0
                        dvar1=dG1;
                    elseif outdata==2
                        dvar1=dG1;
                        dvar2=dH1;
                    end
                elseif Sa(R)<Sb(R)
                    if outdata==1
                        dvar1=dH2;
                    elseif outdata==0
                        dvar1=dG2;
                    elseif outdata==2
                        dvar1=dG2;
                        dvar2=dH2;
                    end
                end
                %Here we imposed interp2 for smaller earthquake and grid fit for biggest because of the distortion undergone
                %interp2 is faster
                if Mw(R)<=6.5
                    disp('Mw(R)<=6.5')
                    if Mw(R)<=6.5 && Mw(R)>=6
                    end
                    [latg,long]=utm2ll(lon1,lat1,zone(R));
                    [longR1,latgR1]=meshgrid(long,latg);
                    long = wrapTo180(long);
                        Gq=interp2(longR1,latgR1,dvar1,longT2,latT2,'linear');
                        Gq(isnan(Gq)==1)=0;
                        if outdata==2
                        Gq2=interp2(longR1,latgR1,dvar2,longT2,latT2,'linear');
                        Gq2(isnan(Gq2)==1)=0;
                        end
                else
                        [lon2,lat2]=meshgrid(lon1,lat1);
                        %Re-convert from meters to decimal degree using utm2ll fonction
                        [LA2,LO2] = utm2ll(lon2(1:numel(lon2)),lat2(1:numel(lat2)),zone(R));
                        indlocallon=find(wrapTo360(longT2(1,:))>=wrapTo360(min(LO2)) & wrapTo360(longT2(1,:))<=wrapTo360(max(LO2)));
                        %indlocallon=find(wrapTo360(longT2(1,:))>=min(LO2) & wrapTo360(longT2(1,:))<=max(LO2));
                        indlocallat=find(latT2(:,1)>=min(LA2) & latT2(:,1)<=max(LA2));
                        LO2 = wrapTo180(LO2);
                        Gqlocal = gridfit(LO2,LA2,dvar1(1:numel(dvar1)),longT2(1,indlocallon),latT2(indlocallat,1),'tilesize',round(numel(longT{Mnum})/coef_tile),'overlap',0.25);
                        Gqlocal(isnan(Gqlocal)==1)=0;
                        Gq=zeros(size(longT2));
                        Gq(indlocallat,indlocallon)=Gq(indlocallat,indlocallon)+Gqlocal;
                    if outdata==2
                    Gqlocal2 = gridfit(LO2,LA2,dvar2(1:numel(dvar2)),longT2(1,indlocallon),latT2(indlocallat,1),'tilesize',round(numel(longT{Mnum})/coef_tile),'overlap',0.25);
                    Gqlocal2(isnan(Gqlocal2)==1)=0;
                    Gq2=zeros(size(longT2));
                    Gq2(indlocallat,indlocallon)=Gq2(indlocallat,indlocallon)+Gqlocal2;
                    end
                end
            end
            toc
        end
        %% Sum everything in one matrix
        matot=matot+Gq;
        if outdata==2
         matot2=matot2+Gq2;
        end
    end
    [Xq,Yq]=meshgrid(longT{1},latT{1});
    matnew{Mnum}=interp2(longT{Mnum},latT{Mnum},matot,Xq,Yq,'linear');
     if outdata==2
    matnew2{Mnum}=interp2(longT{Mnum},latT{Mnum},matot2,Xq,Yq,'linear');
     end
end
catmatnew=cat(3,matnew{:});
matotal=sum(catmatnew,3);
 if outdata==2
catmatnew2=cat(3,matnew2{:});
matotaldH=sum(catmatnew2,3);
 end

%% Plot and projection
utmlonlim=[lonlimwest lonlimest];
utmlatlim=[latlimdown latlimup];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fig points
if strcmp(buttonvar,'dG')
 m_proj('Miller Cylindrical','long',utmlonlim,'lat',utmlatlim);
%m_proj('Miller Cylindrical','long',[119 123],'lat',[21.5 25.4]);%Taiwan
%m_proj('Miller Cylindrical','long',utmlonlim,'lat',[latlimdown latlimup-5]);
% fig points
set(figure,'Position',[1 1 2000 1500]) % [coin_x coin_y hauteur(px) largeur(px)]
set(gcf,'PaperPositionMode','auto')
m_pcolor(longT{1},latT{1},1e8.*matotal),shading interp;%flat;
% [cs,h]=m_contour(longT{1},latT{1},1e8.*matotal)
% clabel(cs,h,'fontsize',6);
% set(gca, 'CLim', [-10, 10]);
colormap(jet);
%colormap(b2r(min(min(1e8.*matotal))),max(max(1e8.*matotal)));
m_coast('color','k','LineWidth',1);
m_grid('box','fancy','tickdir','out');
colorbar;h = colorbar;ylabel(h, '\delta\muGal','FontSize',12)
dim = [.2 .5 .3 .3];
Models = {'Leo10a','Leo10b','Leo10c','YenMa11a','YenMa11b','YenMa11c','WC94a','WC94b'};
Nodals = {'Nodale 1','Nodale 2','Mixte'};
str = {['Mod�le :' ' ' cell2mat(Models(Model)) ',' ' ' 'Nodale :' ' ' cell2mat(Nodals(nodal)) '-' button1]};
ah=annotation('textbox',dim,'units','normalized','position',[.32 .867 .1 .05],'String',...
    str,'FitBoxToText','on','backgroundcolor',[0.9 0.9 0.9]);
title('','FontSize',13);
end


if strcmp(buttonvar,'dH')
m_proj('Miller Cylindrical','long',utmlonlim,'lat',[latlimdown latlimup-5]);
% fig points
set(figure(1),'Position',[1 1 2000 1500]) % [coin_x coin_y hauteur(px) largeur(px)]
set(gcf,'PaperPositionMode','auto')
m_pcolor(longT{1},latT{1},matotal),shading interp;%flat;
% [cs,h]=m_contour(longT{1},latT{1},matotal)
% clabel(cs,h,'fontsize',6);
%set(gca, 'CLim', [-10, 10]);
colormap(jet);
%colormap(b2r(min(min(matotal))),max(max(matotal)));
m_coast('color','k','LineWidth',1);
m_grid('box','fancy','tickdir','out');
colorbar;h = colorbar;ylabel(h, '\deltam','FontSize',12)
dim = [.2 .5 .3 .3];
Models = {'Leo10a','Leo10b','Leo10c','YenMa11a','YenMa11b','YenMa11c','WC94a','WC94b'};
Nodals = {'Nodale 1','Nodale 2','Mixte'};
str = {['Mod�le :' ' ' cell2mat(Models(Model)) ',' ' ' 'Nodale :' ' ' cell2mat(Nodals(nodal)) '-' button1]};
ah=annotation('textbox',dim,'units','normalized','position',[.32 .867 .1 .05],'String',...
    str,'FitBoxToText','on','backgroundcolor',[0.9 0.9 0.9]);
title('','FontSize',13);
end

if strcmp(buttonvar,'Both')
m_proj('Miller Cylindrical','long',[-180 180],'lat',[-70 75]);
%m_proj('Miller Cylindrical','long',utmlonlim,'lat',[latlimdown latlimup-5]);
% fig points
set(figure,'Position',[1 1 2000 1500]) % [coin_x coin_y hauteur(px) largeur(px)]
set(gcf,'PaperPositionMode','auto')
m_pcolor(longT{1},latT{1},1e8.*matotal),shading interp;%flat;
% [cs,h]=m_contour(longT{1},latT{1},1e8.*matotal)
% clabel(cs,h,'fontsize',6);
%set(gca, 'CLim', [-10, 10]);
colormap(jet);
%colormap(b2r(min(min(1e8.*matotal))),max(max(1e8.*matotal)));
m_coast('color','k','LineWidth',1);
m_grid('box','fancy','tickdir','out');
colorbar;h = colorbar;ylabel(h, '\delta\muGal','FontSize',12)
dim = [.2 .5 .3 .3];
Models = {'Leo10a','Leo10b','Leo10c','YenMa11a','YenMa11b','YenMa11c','WC94a','WC94b'};
Nodals = {'Nodale 1','Nodale 2','Mixte'};
str = {['Mod�le :' ' ' cell2mat(Models(Model)) ',' ' ' 'Nodale :' ' ' cell2mat(Nodals(nodal)) '-' button1]};
ah=annotation('textbox',dim,'units','normalized','position',[.32 .867 .1 .05],'String',...
    str,'FitBoxToText','on','backgroundcolor',[0.9 0.9 0.9]);
title('','FontSize',13);

%m_proj('Miller Cylindrical','long',utmlonlim,'lat',[latlimdown latlimup-5]);
m_proj('Miller Cylindrical','long',[-180 180],'lat',[-70 80]);
% fig points
set(figure,'Position',[1 1 2000 1500]) % [coin_x coin_y hauteur(px) largeur(px)]
set(gcf,'PaperPositionMode','auto')
m_pcolor(longT{1},latT{1},matotaldH),shading interp;%flat;
% [cs,h]=m_contour(longT{1},latT{1},matotaldH)
% clabel(cs,h,'fontsize',6);
%set(gca, 'CLim', [-10, 10]);
colormap(jet);
%colormap(b2r(min(min(matotaldH))),max(max(matotaldH)));
m_coast('color','k','LineWidth',1);
m_grid('box','fancy','tickdir','out');
colorbar;h = colorbar;ylabel(h, '\deltam','FontSize',12)
dim = [.2 .5 .3 .3];
Models = {'Leo10a','Leo10b','Leo10c','YenMa11a','YenMa11b','YenMa11c','WC94a','WC94b'};
Nodals = {'Nodale 1','Nodale 2','Mixte'};
str = {['Mod�le :' ' ' cell2mat(Models(Model)) ',' ' ' 'Nodale :' ' ' cell2mat(Nodals(nodal)) '-' button1]};
ah=annotation('textbox',dim,'units','normalized','position',[.32 .867 .1 .05],'String',...
    str,'FitBoxToText','on','backgroundcolor',[0.9 0.9 0.9]);
title('','FontSize',13);

end

load handel %gong%chirp %handel
sound(y,Fs)

%% Save
if outdata==0
    matnewdG=matnew;
    matotaldG=matotal;
    save([dname,'\matnewdG.mat'],'matnewdG')
    save([dname,'\longT.mat'],'longT')
    save([dname,'\latT.mat'],'latT')
    save([dname,'\matotaldG.mat'],'matotaldG')
end

if outdata==1
    matnewdH=matnew;
    matotaldH=matotal;
    save([dname,'\matnewdH.mat'],'matnewdH')
    save([dname,'\longT.mat'],'longT')
    save([dname,'\latT.mat'],'latT')
    save([dname,'\matotaldH.mat'],'matotaldH')
end

if outdata==2
    matnewdG=matnew;
    matotaldG=matotal;
    matnewdH=matnew2;
    save([dname,'\matnewdH.mat'],'matnewdH')
    save([dname,'\matotaldH.mat'],'matotaldH')
    save([dname,'\matnewdG.mat'],'matnewdG')
    save([dname,'\matotaldG.mat'],'matotaldG')
    save([dname,'\longT.mat'],'longT')
    save([dname,'\latT.mat'],'latT')
end

%% Plot Area Extension
% m_proj('Miller Cylindrical','long',utmlonlim,'lat',utmlatlim);
%     set(figure,'Position',[1 1 2000 1500]) % [coin_x coin_y hauteur(px) largeur(px)]
%     set(gcf,'PaperPositionMode','auto')
%     m_pcolor(longT{1},latT{1},1e8.*matnew{i}),shading flat;
%     colormap(b2r(min(min(1e8.*matnew{i})),max(max(1e8.*matnew{i}))));
%     m_coast('color','k','LineWidth',1);
%     h=colorbar('h');
% %m_grid('box','fancy','tickdir','in');
% hold on;
% for R=1:numel(ind1)
%     disp([num2str(R),'/',num2str(numel(ind1)),'Earthquake area'])
%     m_plot([minlong(R),minlong(R),maxlong(R),maxlong(R),minlong(R)],[minlat(R),maxlat(R),maxlat(R),minlat(R),minlat(R)]);
%     hold on
%     m_grid
%     %plot(minlong(2),minlong(2),maxlong(2),maxlong(2),minlong(2),minlat(2),maxlat(2),maxlat(2),minlat(2),minlat(2));
% end
%% References
% Leonard, M. (2010). Earthquake fault scaling: Self-consistent relating of rupture length, width, average displacement, and moment release. Bulletin of the Seismological Society of America, 100(5A), 1971-1988.
% Wells, D. L., & Coppersmith, K. J. (1994). New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement. Bulletin of the seismological Society of America, 84(4), 974-1002.
% Yen, Y. T., & Ma, K. F. (2011). Source-scaling relationship for M 4.6�8.9 earthquakes, specifically for earthquakes in the collision zone of Taiwan. Bulletin of the Seismological Society of America, 101(2), 464-481.
% Okubo, S. (1992). Gravity and potential changes due to shear and tensile faults in a half?space. Journal of Geophysical Research: Solid Earth, 97(B5), 7137-7144.
