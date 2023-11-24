%% initialize
clc; clearvars; close all;

set(0,'defaultfigureposition',[500 500 800 600]') % x, y , w , h
set(0,'defaultAxesFontSize',20)
set(0,'defaultAxesFontName','Arial')
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultFigureRenderer', 'painters');

newsavefigs = 1;

% folder name should end with '/'
% name of subfolder in output folder that contains the data of interest
%datatype = 'norm_membratio_BdytoLumen/';

% only boundary signal
% datatype = 'NOlumen_60plperp_SAVEALL_v2022-11-22/';

% boundary:Lumen ratio, but not membrane marker normalized
% datatype = 'NOmnorm_60plperp_SAVEALL_v2022-11-07/';

% elife paper
% datatype = 'elife_NOmnorm_60plperp_v2022-11-07/';
% datatype = 'elifeg2_NOlumen_60plperp_SAVEALL_v2022-12-27/';

% use for local guvs bc boundary error prone, so 100plperp used
% for PCA analysis use boundary:lumen ratio but not membran marker normalized "NOmnorm"
% for other boundary analysis use only boundary signal "NOLumen
%datatype = 'NOmnorm_100plperp_SAVEALL_v2022-11-27/';
datatype = 'NOlumen_100plperp_SAVEALL_v2023-01-07/';

% list of kymographs data folder for each GUV
%S = dir('elife*testout'); % list of Global output folders
%S = dir('elifeg2*testout'); % list of Global output folders
%S = dir('G*output'); % list of Global output folders
S = dir('L*output'); % list of local output folders

%S = dir(strcat('G*output/',datatype,'*data.mat')); 

% folder name should end with '/'
newoutputfolder = strcat('analysis19_LOC_ActinRate_',datatype(1:end-1),'_75pthr_v1/');
%newoutputfolder = 'analysis16_elife_NOmnorm_60plperp_v2022-11-07_vary_75pthr_v1/';

%%%%%% CHECK IF NAME MATCHES INPUTS !!!!!!!!!!!!!!
%-========================
%Gidx_analysis = [2:4,6:8,10,11,14,18]';% NON-DEFORMED GUVs
%Gidx_analysis = [1,5,9,12,13,15:17]'; % DEFORMED GUVS

%Gidx_analysis = [1,3,4,5,6,8,15:17]';
%Gidx_analysis = [1,5,15:17]'; % DEFORMED GUVS

%Gidx_analysis = [1:18]'; % for correlations
%Gidx_analysis = [2:18]'; % FOR PCA and signal corr
%Gidx_analysis = [1:3]; % FOR local and eLife analysis


%Gidx_analysis = [1:18]';
%Gidx_analysis = [2,4,11:18]';
Gidx_analysis = [3]';
%Gidx_analysis = [3,5,6]';

outputname = newoutputfolder(9:end-1);
tstamp_analysis = datestr(now,'dd_mmm_yyyy_HH_MM_SS');

%% %%% SETTTINGS
% default is 1, to use fixed kymos (360 points) so larger guvs do not get weighted more
do_fixedspcorr = 1; % fixedkymos use 360 spatial bins
do_fixedsigcorr = 1; % fixedkymos use 360 spatial bins
do_norm_corr = 0; % 0: no normalization-
                % 1: mean norm; 
                % 2: unit norm. % requires using fixed kymo for L2 norm of columns
                % 3: max norm.
                % 4: max norm every frame
do_partial_corr = 0;
nthresh = 180; % usually 100, just make sure it's at least half # of spatial bins
pthresh = 0.75; % usu 0.75; set to 1 if doing partial_spcorr!
% default 1: use top 100% of membrane marker points, (no exclusion of high membrane marker)
remove_membnorm = 0; % default is 0: no post-hoc multiplication by membrane marker kymograph to cancel out memb.norm.
add_membnorm = 0; % default is 0: no post-hoc division by membrane marker kymo
do_absvals = 0;
all_smooth = 0; % smooth all data initially before analysis
sfilt_all = 5; % half of window size for mov avg smooth all data initially before analysis

show_signalcorr = 0;
show_spatialcorr = 0;
show_alignment = 0; % alignment of mean-norm ActA and actin at a single frame pre/post rapa
show_entropy = 0;
show_variance = 0;
show_pca = 0;
show_eccent_actin_corr = 0;

show_actinrate = 1; % rate of actin growth radially in Local rapamycin expts



%% ==========================

n_analysis = length(Gidx_analysis); % # of GUVs to be analyzed


fileparams = v2struct;
if exist(newoutputfolder,'dir') == 0
    mkdir(newoutputfolder);
end

%% extract data from multiple GUVs




df_list = cell(n_analysis,1); % list of .mat files for SELECT GUVs
df = cell(n_analysis,1); % variables in each .mat file 

% compile quantities from multiple GUVs
dParams = cell(n_analysis,1);
dTrapa = zeros(n_analysis,1); % time after which Rapamycin is present (for alignment of plots)
dNframes = zeros(n_analysis,1); % total frames for each guv
dguvnames = strings(n_analysis,1);

dKBmemb = cell(n_analysis,1);
dKBActA = cell(n_analysis,1);
dKBactin = cell(n_analysis,1);
dKBvel = cell(n_analysis,1);
dKBdc = cell(n_analysis,1);
dKBcurv = cell(n_analysis,1);

dKVmemb = cell(n_analysis,1);
dKVActA = cell(n_analysis,1);
dKVactin = cell(n_analysis,1);
dKVvel = cell(n_analysis,1);
dKVdc = cell(n_analysis,1);
dKVcurv = cell(n_analysis,1);

% global shape parameters
dEccent = cell(n_analysis,1);
dAspratio_ml = cell(n_analysis,1);
dAspratio_ind = cell(n_analysis,1);
dAreas = cell(n_analysis,1);
dMAlengths = cell(n_analysis,1); % major axis lengths

for i = 1:n_analysis
    tmpfolder = fullfile(S(Gidx_analysis(i)).name,datatype);
    tmpfile = dir([tmpfolder,'*.mat']);
    df_list{i} = fullfile(tmpfolder,tmpfile.name);
    % df{i} = load(df_list{i},'dataKymoBin','dataKymoVar','dataParams');
    df{i} = load(df_list{i},'dataKymoBin','dataKymoVar','dataParams','dataRaw');

%     tmpfolder = S(Gidx_analysis(i)).folder;
%     tmpfile = S(Gidx_analysis(i)).name;
%     df_list{i} = fullfile(tmpfolder,tmpfile);
%     df{i} = load(df_list{i});

    dKBmemb{i} = df{i}.dataKymoBin.kymofluor{1};
    dKBActA{i} = df{i}.dataKymoBin.kymofluor{3}; %ActA is 3rd channel!
    dKBactin{i} = df{i}.dataKymoBin.kymofluor{2};
    dKBvel{i} = df{i}.dataKymoBin.kymovel;
    dKBdc{i} = df{i}.dataKymoBin.kymodcumul;
    dKBcurv{i} = df{i}.dataKymoBin.kymocurv;

    dKVmemb{i} = df{i}.dataKymoVar.kymofluor{1};
    dKVActA{i} = df{i}.dataKymoVar.kymofluor{3};%ActA is 3rd channel!
    dKVactin{i} = df{i}.dataKymoVar.kymofluor{2};
    dKVvel{i} = df{i}.dataKymoVar.kymovel;
    dKVdc{i} = df{i}.dataKymoVar.kymodcumul;
    dKVcurv{i} = df{i}.dataKymoVar.kymocurv;

    if remove_membnorm
        dKBActA{i} = dKBActA{i}.*dKBmemb{i};
        dKBactin{i} = dKBactin{i}.*dKBmemb{i};
        dKVActA{i} = dKVActA{i}.*dKVmemb{i};
        dKVactin{i} = dKVactin{i}.*dKVmemb{i};
    elseif add_membnorm
        dKBActA{i} = dKBActA{i}./dKBmemb{i};
        dKBactin{i} = dKBactin{i}./dKBmemb{i};
        dKVActA{i} = dKVActA{i}./dKVmemb{i};
        dKVactin{i} = dKVactin{i}./dKVmemb{i};
        
    end

    if all_smooth
        dKBActA{i} = mavgsmoothboundary(sfilt_all,dKBActA{i});
        dKBactin{i} = mavgsmoothboundary(sfilt_all,dKBactin{i});
        dKVActA{i} = mavgsmoothboundary(sfilt_all,dKVActA{i});
        dKVactin{i} = mavgsmoothboundary(sfilt_all,dKVactin{i});
        dKBcurv{i} = mavgsmoothboundary(sfilt_all,dKBcurv{i});
        dKVcurv{i} = mavgsmoothboundary(sfilt_all,dKVcurv{i});
    end



    dParams{i} = df{i}.dataParams;
    dTrapa(i) = df{i}.dataParams.t_rapa;
    dNframes(i) = df{i}.dataParams.nframes;
    dguvnames(i) = df{i}.dataParams.guvname;

    dEccent{i} = df{i}.dataRaw.allEccents;
    dAspratio_ml{i} = df{i}.dataRaw.allAspRatio_ml;
    dAspratio_ind{i} = df{i}.dataRaw.allAspRatio_ind;
    dAreas{i} = df{i}.dataRaw.allAreas;
    dMAlengths{i} = df{i}.dataRaw.allMajAxL;
end

% clear the last variables collected
%clear dataKymoBin dataKymoVar dataParams dataRaw

dKymoBin = struct('memb',dKBmemb,'ActA',dKBActA,'actin',dKBactin,...
    'vel',dKBvel,'dc',dKBdc,'curv',dKBcurv);
dKymoVar = struct('memb',dKVmemb,'ActA',dKVActA,'actin',dKVactin,...
    'vel',dKVvel,'dc',dKVdc,'curv',dKVcurv);

% % rows are different GUVs, columns are the 6 kymos
% allKymoBin = struct2cell(dKymoBin)'; 
% allKymoVar = struct2cell(dKymoVar)';

% % these ranges needed for signal entropy which only is based on KymoBin
% psort = 0.005;
% 
% mb0 = cat(2,dKBmemb{:});
% mb0_size = length(mb0(:));
% mb0_sort = sort(mb0(:),'descend','MissingPlacement','last');
% mb0_max = mb0_sort(round(psort*mb0_size)); 
% 
% ac0 = cat(2,dKBActA{:});
% ac0_size = length(ac0(:));
% ac0_sort = sort(ac0(:),'descend','MissingPlacement','last');
% ac0_max = ac0_sort(round(psort*ac0_size)); 
% 
% an0 = cat(2,dKBactin{:});
% an0_size = length(an0(:));
% an0_sort = sort(an0(:),'descend','MissingPlacement','last');
% an0_max = an0_sort(round(psort*an0_size)); 
% 
% vl0 = abs(cat(2,dKBvel{:}));
% vl0_size = length(vl0(:));
% vl0_sort = sort(vl0(:),'descend','MissingPlacement','last');
% vl0_max = vl0_sort(round(psort*vl0_size)); 
% 
% dc0 = abs(cat(2,dKBdc{:}));
% dc0_size = length(dc0(:));
% dc0_sort = sort(dc0(:),'descend','MissingPlacement','last');
% dc0_max = dc0_sort(round(psort*dc0_size)); 
% 
% cr0 = abs(cat(2,dKBcurv{:}));
% cr0_size = length(cr0(:));
% cr0_sort = sort(cr0(:),'descend','MissingPlacement','last');
% cr0_max = cr0_sort(round(psort*cr0_size)); 
% 
% 
% % these ranges needed for signal entropy which only is based on KymoBin
% mb1 = mb0_max;
% ac1 = ac0_max;
% an1 = an0_max;
% vl1 = vl0_max;
% dc1 = dc0_max;
% cr1 = cr0_max;


% these ranges needed for signal entropy which only is based on KymoBin
mb1 = max(max(cat(2,dKBmemb{:})));
ac1 = max(max(cat(2,dKBActA{:})));
an1 = max(max(cat(2,dKBactin{:})));
vl1 = max(max(abs(cat(2,dKBvel{:}))));
dc1 = max(max(abs(cat(2,dKBdc{:}))));
cr1 = max(max(abs(cat(2,dKBcurv{:}))));

% this creates fixed edges for all guv entropy analysis
% histograms and entropies now comparable

rmemb = [0,mb1];
rActA  = [0,ac1];
ractin  = [0,an1];
rvel = [-vl1,vl1];
rdc = [-dc1,dc1];
rcurv = [-cr1,cr1];
pranges = v2struct(rmemb,rActA,ractin,rvel,rdc,rcurv);

[maxf,idmaxframes] = max(dNframes); %GUV with most frames
tickstep = df{idmaxframes}.dataParams.tickstep;
[maxt,~] = max(dTrapa); %max time after rapa is on
maxd = max(maxt-dTrapa); % max time shift needed for before rap
maxf1 = maxf + maxd;

pstatalign = v2struct(dNframes,dTrapa,dguvnames,pranges,maxf,maxf1,maxt,maxd,tickstep);

%lcmframes = lcms(dNframes)
%---------------

%% ACTIN RATE
if show_actinrate

do_smoothAR = 1;
mfilt_AR = 10;

for i = 1:n_analysis
    
    Tmpactin00 = dKBactin{i};
    
    
    I = mavgsmoothboundary(mfilt_AR,Tmpactin00);
    I = imgaussfilt(I,1);
    Tmpactin0 = medfilt2(I,[1 1]);

    figAR(1) = figure(1); imagesc(Tmpactin0); 
    axis fill; 
   title(strcat('L0',num2str(Gidx_analysis)));
   xlabel('Frames (15 sec)');
   ylabel('Angular bins (deg)')    
    drawnow;

    BWactin0 = imbinarize(rescale(Tmpactin0),0.1);
    figure(2); imagesc(BWactin0); colormap gray; axis fill; drawnow;

    % 10 iterations
    BWactin1 = activecontour(Tmpactin0,BWactin0,10,'Chan-Vese','SmoothFactor',1,'ContractionBias',0);

    % Get rid of small holes (noticeable at initial actin polymerization)
    BWactin2 = imfill(BWactin1, 'holes');
    figure(3); imagesc(BWactin2); colormap gray; axis fill; drawnow;

    [bwb_actin0,lab_actin0,n_actin0,~] = bwboundaries(BWactin2,'noholes');
    actinpatches = regionprops(lab_actin0, 'area');
	% Get all the areas
	allActinAreas = [actinpatches.Area];
    [sortedAreas, sortIndexes] = sort(allActinAreas, 'descend');
    % find largest actin patch
    bwb_big = ismember(lab_actin0, sortIndexes(1)); % integer-labeled image

    bwb_tmp = bwb_big > 0; % convert to binary from integer-label of largest actin patch
    bwb_bdy = bwb_actin0{sortIndexes(1)}; % boundary points of largest actin patch

    % savitzky-golay smoothing of boundary
    windowWidth = 2*(20)+1; % must be odd; 
    polynomialOrder = 2; % keep this as 2 for quadratic smoothing of boundary
    bdy_sg = sgsmoothboundary(bwb_tmp,bwb_bdy,windowWidth,polynomialOrder);


    figAR(2) = figure(4); imagesc(bwb_tmp); hold on;
    %plot(bwb_bdy(:,2),bwb_bdy(:,1),'g*','LineWidth',2)
    plot(bdy_sg(:,2),bdy_sg(:,1),'b*','LineWidth',1)
    % title('Original (green), Savitzky-Golay smoothed (blue)',...
    %     'FontSize',10);
    % f50.Position(3:4) = [500 500];
   colormap gray; 
   title(strcat('L0',num2str(Gidx_analysis)));
   xlabel('Frames (15 sec)');
   ylabel('Angular bins (deg)')
    
    % find series of points to fit to based on rows(y values) and columns (x values)
    xrA = [39,74]; % column range
    yrA = [39,201]; % row range
    p0 = bdy_sg(:,2)>xrA(1) & bdy_sg(:,2)<xrA(2) ...
        & bdy_sg(:,1)>yrA(1) & bdy_sg(:,1)<yrA(2);

    x = bdy_sg(p0,2);
    y = bdy_sg(p0,1);
    % Do the fit to a line:
    coefficients = polyfit(x,-y,1);
    slope = coefficients(1);
    intercept = coefficients(2);
    message = sprintf('|Slope| = %.3f deg/frame = %.3f deg/min.',abs(slope), abs(slope*4));

    % Make a fitted line
    x1 = min(x);
    x2 = max(x);
    y1 = min(y);
    y2 = max(y);
    text(x1,y1,message,'Color','r','FontSize',16)
    xFit = x1:x2;
    yFit = polyval(coefficients, xFit);
    plot(xFit, -yFit, 'r-', 'LineWidth', 2);


    % find series of points to fit to based on rows(y values) and columns (x values)
    xrB = [39,62]; %xrA; % column range
    yrB = [250,349]; % row range
    p0 = bdy_sg(:,2)>xrB(1) & bdy_sg(:,2)<xrB(2) ...
        & bdy_sg(:,1)>yrB(1) & bdy_sg(:,1)<yrB(2);

    x = bdy_sg(p0,2);
    y = bdy_sg(p0,1);
    % Do the fit to a line:
    coefficients = polyfit(x,-y,1);
    slope = coefficients(1);
    intercept = coefficients(2);
    message = sprintf('|Slope| = %.3f deg/frame = %.3f deg/min.',abs(slope), abs(slope*4));

    % Make a fitted line
    x1 = min(x);
    x2 = max(x);
    y1 = min(y);
    y2 = max(y);
    text(x1,y2,message,'Color','r','FontSize',16)
    xFit = x1:x2;
    yFit = polyval(coefficients, xFit);

    plot(xFit, -yFit, 'r-', 'LineWidth', 2);
    axis fill; hold off; drawnow;



    % BW_actinpatch = poly2mask(bdy_sg(:,2),bdy_sg(:,1),size(bwb_tmp,1),size(bwb_tmp,2));
    % figure(5); imshow(BW_actinpatch); axis fill; drawnow;
    % Actinedge = edge(BW_actinpatch);    
    % figure(6); imshow(Actinedge); axis fill; drawnow;

    if newsavefigs
        saveas(figAR(1),strcat(newoutputfolder,'actin_rate_',...
            tstamp_analysis,'initkymo_L0',num2str(Gidx_analysis),'.svg'));
        saveas(figAR(2),strcat(newoutputfolder,'actin_rate_',...
            tstamp_analysis,'actinrate_L0',num2str(Gidx_analysis),'.svg'));

        savefig(figAR,strcat(newoutputfolder,'actin_rate_',...
            tstamp_analysis,'_L0',num2str(Gidx_analysis),'.fig'));
    end

end

end

%% ECCENTRICITY and Actin
afs = 14; 
if show_eccent_actin_corr
    
stime = cell(n_analysis,1);
rawsig = cell(n_analysis,1);
dMSignal = cell(n_analysis,1);
dShapeE = cell(n_analysis,1);
lag = 0; % 0 works best; useful to look for delayed correlations


flag_actin = 0; % 0 means use ActA

for i = 1:n_analysis
    % pre or post rapa time frames
    %stime{i} = dTrapa(i):dNframes(i);
    %stime{i} = dTrapa(i)+20:dTrapa(i)+54;
    stime{i} = 1:dTrapa(i);
    if flag_actin
        rawsig{i} = dKVactin{i};
    else
        rawsig{i} = dKVActA{i};
    end
end




for i = 1:n_analysis
    rawky = rawsig{i};


    if do_norm_corr == 1
        rawky = rawky./mean(rawky,'omitnan');
    elseif do_norm_corr == 2
        rawky = rawky./vecnorm(rawky,'omitnan');
    elseif do_norm_corr == 3
        rawky = rawky./max(rawky(:));
    elseif do_norm_corr == 4
        rawky = rawky./max(rawky,[],1);
    end

    tmpky = rawky(:,stime{i});

    tmpmeanS = mean(tmpky,'omitnan');
    dMSignal{i} = tmpmeanS';

    tmpE = dEccent{i};
    dShapeE{i} = tmpE(stime{i}+lag);
end


corrEAtype = 'Pearson';
%corrEAtype = 'Spearman';

dMActAList = cat(1,dMSignal{:});
dSEList = cat(1,dShapeE{:});

f1EActA = figure; 
binscatter(dMActAList,dSEList,60); % 60 bins
%f1EA = scatter(dSignal,dShape,'filled');
[corrEActA,pvalEActA] = corr(dMActAList,dSEList,'type',corrEAtype);
% xL=xlim;
% yL=ylim;
% xT=xL(1)-(xL(2)-xL(1))/6.5;
% yT=yL(2)+(yL(2)-yL(1))/17;
% text(xT,yT,num2str(corrEA));

title(["Mean boundary ActA v Eccentricity",strcat("(Pearson corr: ",num2str(corrEActA),")")])


ax = gca; ax.FontSize = afs;


%================================



flag_actin = 1; % 0 means use ActA

for i = 1:n_analysis

    if flag_actin
        rawsig{i} = dKVactin{i};
    else
        rawsig{i} = dKVActA{i};
    end
end

for i = 1:n_analysis
    rawky = rawsig{i};

    if do_norm_corr == 1
        rawky = rawky./mean(rawky,'omitnan');
    elseif do_norm_corr == 2
        rawky = rawky./vecnorm(rawky,'omitnan');
    elseif do_norm_corr == 3
        rawky = rawky./max(rawky(:));
    elseif do_norm_corr == 4
        rawky = rawky./max(rawky,[],1);
    end

    tmpky = rawky(:,stime{i});
    tmpmeanS = mean(tmpky,'omitnan');
    dMSignal{i} = tmpmeanS';

    tmpE = dEccent{i};
    dShapeE{i} = tmpE(stime{i}+lag);
end

lag = 0; % 0 works best; useful to look for delayed correlations
corrEAtype = 'Pearson';
%corrEAtype = 'Spearman';


dMActinList = cat(1,dMSignal{:});


f1EActin = figure; 
binscatter(dMActinList,dSEList,60);
%f1EA = scatter(dSignal,dShape,'filled');
[corrEActin,pvalEActin] = corr(dMActinList,dSEList,'type',corrEAtype);
% xL=xlim;
% yL=ylim;
% xT=xL(1)-(xL(2)-xL(1))/6.5;
% yT=yL(2)+(yL(2)-yL(1))/17;
% text(xT,yT,num2str(corrEA));

title(["Mean boundary actin v Eccentricity",strcat("(Pearson corr: ",num2str(corrEActin),")")])

ax = gca; ax.FontSize = afs;

%f1EA = ploteccentactin(figEAidx,dSignal,dShapes,pstatalign,fileparams);

if newsavefigs
    print(f1EActA,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_eccentActA'),'-dpng')
    savefig(f1EActA,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_eccentActA.fig'),'compact');

    print(f1EActin,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_eccentActin'),'-dpng')
    savefig(f1EActin,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_eccentActin.fig'),'compact');
end

end


%% SETTINGS for pca 
if show_pca == 1

do_PCAonfullkymo = 1; %
do_norm_pca = 3; % 0: no normalization- USE THIS IF BINARIZED data especially
                % 1: mean norm; 
                % 2: unit norm.
                % 3: max norm.
do_binarized = 0; % binarize the kymograph
do_spectimes = 1; % 0: use frames 1 : end-3
                % 1: use select frames
do_align_style = 2; % 0: no second round of alignment
                    % 1: center the max at each frame
                    % 2: center the max at frame trapa + halfway to the end
use_actin_ref = 1; % use actin as a reference kymo for alignment
do_std_clims = 1; % use standardized coloraxis limits 
do_smooth = 1; % moving average smoothing
do_meanbin = 1;
sfilt = 10; 

which_qty = "actin"; 
%which_qty = "curvature";
%which_qty = "ActA";

if strcmp(which_qty,"actin")
    use_actin_ref = 0; 
end

pcadim = 6; % 6 for global

%% alignment for pca

times = cell(n_analysis,1);
idmax = cell(n_analysis,1);



Xdata_final=cell(n_analysis,1);
Xdata_scaled=cell(n_analysis,1);
Xdata_orig = cell(n_analysis,1);

Xdata_final2=cell(n_analysis,1);
Xdata_scaled2=cell(n_analysis,1);

for i = 1:n_analysis

    if strcmp(which_qty,"actin")
        Xtmpinit = dKBactin{i};
    elseif strcmp(which_qty,"curvature")
        Xtmpref = dKBactin{i}; % for alignment
        Xtmpinit = dKBcurv{i};
    elseif strcmp(which_qty,"ActA")
        Xtmpref = dKBactin{i}; % for alignment
        Xtmpinit = dKBActA{i};
    end



    if do_meanbin
        Xtmp0 = nan(180,size(Xtmpinit,2));
        for j = 1:180
            Xtmp0(j,:) = mean(Xtmpinit(2*j-1:2*j,:));
        end
        % halfway length
        cnt1 = 90;
    else
        Xtmp0 = Xtmpinit;
        % halfway length
        cnt1 = 180;
    end

    if do_smooth
        Xtmp00 = mavgsmoothboundary(sfilt,Xtmp0);
    else
        Xtmp00 = Xtmp0;
    end


    if strcmp(which_qty,"ActA") || strcmp(which_qty,"curvature")
        % do binning and smoothing on the reference data from Actin 
        if do_meanbin
            Xtmp0ref = nan(180,size(Xtmpref,2));
            for j = 1:180
                Xtmp0ref(j,:) = mean(Xtmpref(2*j-1:2*j,:));
            end
            % halfway length
            cnt1 = 90;
        else
            Xtmp0ref = Xtmpref;
            % halfway length
            cnt1 = 180;
        end
    
        if do_smooth
            Xtmp00ref = mavgsmoothboundary(sfilt,Xtmp0ref);
        else
            Xtmp00ref = Xtmp0ref;
        end
        

    end

%     figure; imagesc(Xtmpinit); colorbar;
%     figure; imagesc(Xtmp0); colorbar


    %t = dt:dt:maxtime-dt;
%     t = 1:dNframes(i)-1;
%     n = length(t); 
%     fhat = fft(Xtmp0(:,1),n); % Compute the fast Fourier transform
%     PSD = fhat.* conj(fhat)/n; % Power spectrum (power per freq)
%     freq = 1/(1*n)*(0:n); % Create x-axis of frequencies in Hz
%     L = 1:floor(n/2); % Only plot the first half of freqs
% 
%     figure(10); semilogy(freq(L),PSD(L)); 



    % select time points for PCA

    if do_spectimes % use select frames for pca
        %times{i} = [1:dTrapa(i)]';
        %times{i} = [(dTrapa(i)+5):(dNframes(i)-3)]';     
        %times{i} = [(dTrapa(i)+10):(dTrapa(i)+19)]'; 
        %times{i} = [(dTrapa(i)-5):(dTrapa(i)+50)]'; 

        % for global
        times{i} = [(dTrapa(i)+1):(dTrapa(i)+50)]'; 

        % for local
        %times{i} = [(dTrapa(i)+1):(dTrapa(i)+250)]'; 

    else % use all frames for pca
        times{i} = [1:(dNframes(i)-3)]';
    end



    % normalization of target kymo 
    % to ensure features are comparable across guvs
    if do_norm_pca == 1
        Xtmp = Xtmp00./mean(Xtmp00);
    elseif do_norm_pca == 2
        Xtmp = Xtmp00./vecnorm(Xtmp00);
    elseif do_norm_pca == 3
        % here normalize over the specific time window selected
        % WORKS BETTER
        xxtmp = Xtmp00(:,times{i});        
        Xtmp = Xtmp00./max(xxtmp(:));
    else
        Xtmp = Xtmp00;
    end
    % for normalization of curvature, use 1/(0.5*major axis length of first frame)
    if strcmp(which_qty,"curvature")
        tmpMAlen = dMAlengths{i};
        tt = times{i}(1);
        tmpcurvinit = 1/(0.5*tmpMAlen(tt));
        Xtmp = Xtmp00./tmpcurvinit;
%        Xtmp = Xtmp00;
    end


    % normalization over full time window - NOT AS GOOD
%     if do_norm_pca == 1
%         Xtmp = Xtmp00./mean(Xtmp00);
%     elseif do_norm_pca == 2
%         Xtmp = Xtmp00./vecnorm(Xtmp00);
%     elseif do_norm_pca == 3
%         %xxtmp = Xtmp00;        
%         Xtmp = Xtmp00./max(Xtmp00(:));
%     else
%         Xtmp = Xtmp00;
%     end


    % normalize the ref actin kymograph if doing pca on other kymos
    if strcmp(which_qty,"ActA") || strcmp(which_qty,"curvature")
        % normalization of reference kymo (actin) to ensure features are comparable across guvs
        if do_norm_pca == 1
            XtmpAref = Xtmp00ref./mean(Xtmp00ref);
        elseif do_norm_pca == 2
            XtmpAref = Xtmp00ref./vecnorm(Xtmp00ref);
        elseif do_norm_pca == 3
            % here normalize over the specific time window selected
            % WORKS BETTER
            xxtmpAref = Xtmp00ref(:,times{i}); 
            XtmpAref = Xtmp00ref./max(xxtmpAref(:));
        else
            XtmpAref = Xtmp00ref;
        end
        % turn off if doing binarization
        if do_binarized
            XtmpAref = Xtmp00ref;
        end

    end


    % PCA works better with alignment
    % align boundary quantity at each frame so that max is in middle
    tmpncol = size(Xtmp,2);   
    Xtmp_align = nan(size(Xtmp));

    if do_align_style  == 1
        tmpid = nan(tmpncol,1); % store this for rotating the boundary signal

        for icol = 1:tmpncol
            if use_actin_ref
                [~,tmpid(icol)] = max(XtmpAref(:,icol));
            else
                [~,tmpid(icol)] = max(Xtmp(:,icol));
            end
        
            Xtmp_align(:,icol) = circshift(Xtmp(:,icol),cnt1-tmpid(icol));
    
        end
    elseif do_align_style == 2
        % here align based on actin kymo max at a specific frame

        %tmptime = dTrapa(i)+round((dNframes(i)-dTrapa(i))/2);

        tmptime = dTrapa(i)+35; % for global guvs
        %tmptime = dTrapa(i)+200; % for local
        
        if use_actin_ref
            [~,tmpid] = max(XtmpAref(:,tmptime));
        else
            [~,tmpid] = max(Xtmp(:,tmptime));

        end



        Xtmp_align = circshift(Xtmp,cnt1-tmpid);
    else
        tmpid = nan;
        Xtmp_align = Xtmp;
    end




    XtmpB_norm = Xtmp; % unaligned target kymo
    XtmpB_align = Xtmp_align; % aligned target kymo




    Xdata_binsm = Xtmp00(:,times{i}); % binned, smoothed
    Xdata_norm = XtmpB_norm(:,times{i}); % above + normalized
    Xdata_align = XtmpB_align(:,times{i}); % above + aligned

    



    idmax{i} = tmpid; % store this for rotating any other signals

    Xdata_final{i} = Xdata_align;
    Xdata_scaled{i} = Xdata_norm;
    Xdata_orig{i} = Xdata_binsm;

    Xdata_final2{i} = Xdata_align(:);
    Xdata_scaled2{i} = Xdata_norm(:);
end


%% pca analysis  
afs = 14;

if do_PCAonfullkymo

Zscaled = cat(2,Xdata_scaled2{:});
Zfinal = cat(2,Xdata_final2{:});

Zscaled1 = cat(2,Xdata_scaled{:});
Zfinal1 = cat(2,Xdata_final{:});

else

Zscaled = cat(2,Xdata_scaled{:});
Zfinal = cat(2,Xdata_final{:});

end

Zorig = cat(2,Xdata_orig{:});


f1pca(1) = figure;
imagesc(Zorig)
    if do_std_clims
        tmpmaxval = max(Zorig(:));
        tmpminval = min(Zorig(:));
        clim([tmpminval, tmpmaxval])
    end
ax = gca; ax.FontSize = afs;
% if do_binarized 
%     clrs = [0.2422 0.1504 0.6603; 0.9290 0.6940 0.1250];
%     ax.CLim = [0,1];
%     colormap(ax,clrs);
%     colorbar;
% else
%     colorbar;
% end
colorbar
title('Original kymographs')


if do_PCAonfullkymo

 f1pca(2) = figure;
    imagesc(Zscaled1)
    if do_std_clims
        tmpmaxval = max(Zscaled1(:));
        tmpminval = min(Zscaled1(:));
        clim([tmpminval, tmpmaxval])
    end
    ax = gca; ax.FontSize = afs;
    colorbar
    title('Scaled (NOT vectorized) kymographs')


f1pca(3) = figure;
    imagesc(Zfinal1)
    if do_std_clims
        tmpmaxval = max(Zfinal1(:));
        tmpminval = min(Zfinal1(:));
        clim([tmpminval, tmpmaxval])
    end
    ax = gca; ax.FontSize = afs;
    colorbar
    title('Scaled, aligned (NOT vectorized) kymographs')


end


f1pca(4) = figure;
imagesc(Zscaled)
    if do_std_clims
        tmpmaxval = max(Zscaled(:));
        tmpminval = min(Zscaled(:));
        clim([tmpminval, tmpmaxval])
    end
ax = gca; ax.FontSize = afs;
colorbar;
if do_PCAonfullkymo
title('Scaled (vectorized) kymographs')
else
title('Scaled kymographs')
end

f1pca(5) = figure;
imagesc(Zfinal)
    if do_std_clims
        tmpmaxval = max(Zfinal(:));
        tmpminval = min(Zfinal(:));
        clim([tmpminval, tmpmaxval])
    end
ax = gca; ax.FontSize = afs;
colorbar
if do_PCAonfullkymo
title('Scaled, aligned (vectorized) kymographs')
else
title('Scaled, aligned kymographs')
end


 % PCA computed and plotted
Zm = mean(Zfinal,2);
[Zproj,Zsvectors,percentvar,Zsvalues] = pcabases(Zfinal,pcadim);


f1pca(6) = figure;
Zfinalrel = Zfinal-Zm;
imagesc(Zfinalrel)
    if do_std_clims
        tmpmaxval = max(Zfinalrel(:));
        tmpminval = min(Zfinalrel(:));
        clim([tmpminval, tmpmaxval])
    end
ax = gca; ax.FontSize = afs;
colorbar
if do_PCAonfullkymo
title('Scaled, aligned, mean-subtracted (vectorized) kymographs')
else
title('Scaled, aligned,mean-subtracted kymographs')
end

f1pca(7) = figure;
imagesc(Zsvectors(:,1:2))
colorbar
ax = gca; ax.FontSize = afs;
title('Principal components')

f1pca(8) = figure;
plot(Zsvectors(:,1:2),'LineWidth',3)
ax = gca; ax.FontSize = afs;
hold on
plot(Zm);
hold off
legend(["PC1";"PC2"])

if do_PCAonfullkymo
f1pca(9) = figure;
Zmreshape = reshape(Zm,size(Xdata_align));
imagesc(Zmreshape)
    if do_std_clims 
        tmpmaxval = max(Zfinal(:));
        tmpminval = min(Zfinal(:));
        clim([tmpminval, tmpmaxval])
    end
ax = gca; ax.FontSize = afs;
colorbar
title('Mean kymograph')

f1pca(10) = figure;
Z1reshape = reshape(Zsvectors(:,1),size(Xdata_align));
imagesc(Z1reshape)
    if do_std_clims
        tmpmaxval = max(Zsvectors(:,1:3),[],'all');
        tmpminval = min(Zsvectors(:,1:3),[],'all');
        clim([tmpminval, tmpmaxval])
    end
ax = gca; ax.FontSize = afs;
colorbar
title('PC1 kymograph')

f1pca(11) = figure;
Z2reshape = reshape(Zsvectors(:,2),size(Xdata_align));
imagesc(Z2reshape)
    if do_std_clims
        tmpmaxval = max(Zsvectors(:,1:3),[],'all');
        tmpminval = min(Zsvectors(:,1:3),[],'all');
        clim([tmpminval, tmpmaxval])
    end
ax = gca; ax.FontSize = afs;
colorbar
title('PC2 kymograph')

f1pca(12) = figure;
Z3reshape = reshape(Zsvectors(:,3),size(Xdata_align));
imagesc(Z3reshape)
    if do_std_clims
        tmpmaxval = max(Zsvectors(:,1:3),[],'all');
        tmpminval = min(Zsvectors(:,1:3),[],'all');
        clim([tmpminval, tmpmaxval])
    end
ax = gca; ax.FontSize = afs;
colorbar
title('PC3 kymograph')

f1pca(13) = figure;
Z4reshape = reshape(Zsvectors(:,4),size(Xdata_align));
imagesc(Z4reshape)
    if do_std_clims
        tmpmaxval = max(Zsvectors(:,1:3),[],'all');
        tmpminval = min(Zsvectors(:,1:3),[],'all');
        clim([tmpminval, tmpmaxval])
    end
ax = gca; ax.FontSize = afs;
colorbar
title('PC4 kymograph')

f1pca(14) = figure;
Z5reshape = reshape(Zsvectors(:,5),size(Xdata_align));
imagesc(Z5reshape)
    if do_std_clims
        tmpmaxval = max(Zsvectors(:,1:3),[],'all');
        tmpminval = min(Zsvectors(:,1:3),[],'all');
        clim([tmpminval, tmpmaxval])
    end
ax = gca; ax.FontSize = afs;
colorbar
title('PC5 kymograph')

f1pca(15) = figure;
Z6reshape = reshape(Zsvectors(:,6),size(Xdata_align));
imagesc(Z6reshape)
    if do_std_clims
        tmpmaxval = max(Zsvectors(:,1:3),[],'all');
        tmpminval = min(Zsvectors(:,1:3),[],'all');
        clim([tmpminval, tmpmaxval])
    end
ax = gca; ax.FontSize = afs;
colorbar
title('PC6 kymograph')
end

% figure;
% plot(Zsvectors(:,4:6),'LineWidth',3)
% ax = gca; ax.FontSize = afs;
% hold on
% plot(Zm);
% hold off
% legend(["PC4";"PC5";"PC6"])

% figure;
% imagesc(5*Zsvectors(:,1:pcadim)+Zm)
% ax = gca; ax.FontSize = afs;


if newsavefigs
    print(f1pca(9),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_mean'),'-dpng')
    print(f1pca(10),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_PC1'),'-dpng')
    print(f1pca(11),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_PC2'),'-dpng')
    print(f1pca(12),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_PC3'),'-dpng')
    savefig(f1pca(9:12),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_PCtop.fig'),'compact');
end


end

%% VARIANCE and shape

if show_variance
% varshapes = v2struct(dEccent,dAspratio_ml,dAspratio_ind,dAreas);
dVariance = variancecomp(dParams,dKymoBin,pranges,fileparams);

%dEntrp: 18 guvs (rows) x 4 biochem entropy (ActA, actin, dcumul, curv)

figVidxA = 3;
f1aV = plotvariancelocalshape(figVidxA,dVariance,pstatalign,fileparams);

figVidxB = 4;
%dEntrpChem: 18 guvs (rows) x 2 biochem entropy (ActA, actin)
dVarChem = dVariance(:,1:2);

dShapes = cell(n_analysis,1); % 4 global shape features
dShapes(:,1) = dEccent;
% dShapes(:,2) = dAspratio_ml;
% dShapes(:,3) = dAspratio_ind;
% 
% normdAreas = cell(size(dAreas));
% % normalize areas to area during first frame
% for i = 1:n_analysis
%     tmparea = dAreas{i};
%     tmparea = tmparea/tmparea(1); 
%     normdAreas{i} = tmparea;
% end
% dShapes(:,4) = normdAreas;



f1bV = plotvarianceglobalshape(figVidxB,dVarChem,dShapes,pstatalign,fileparams);

if newsavefigs
    print(f1aV,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_variancelocalshape'),'-dpng')
    savefig(f1aV,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_variancelocalshape.fig'),'compact');
    print(f1bV,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_varianceglobalshape'),'-dpng')
    savefig(f1bV,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_varianceglobalshape.fig'),'compact');
end

end

%% ENTROPY and shape

if show_entropy == 1
% varshapes = v2struct(dEccent,dAspratio_ml,dAspratio_ind,dAreas);
dEntrp = entropycomp(dParams,dKymoBin,pranges,fileparams);

%dEntrp: 18 guvs (rows) x 4 biochem entropy (ActA, actin, dcumul, curv)

figEidxA = 8;
f1a = plotentropylocalshape(figEidxA,dEntrp,pstatalign,fileparams);

figEidxB = 9;
%dEntrpChem: 18 guvs (rows) x 2 biochem entropy (ActA, actin)
dEntrpChem = dEntrp(:,1:2);

dShapes = cell(n_analysis,1); % 4 global shape features
dShapes(:,1) = dEccent;
% dShapes(:,2) = dAspratio_ml;
% dShapes(:,3) = dAspratio_ind;
% 
% normdAreas = cell(size(dAreas));
% % normalize areas to area during first frame
% for i = 1:n_analysis
%     tmparea = dAreas{i};
%     tmparea = tmparea/tmparea(1); 
%     normdAreas{i} = tmparea;
% end
% dShapes(:,4) = normdAreas;



[f1b,Eall,Emed,Sall,Smed,framenums] = plotentropyglobalshape(figEidxB,dEntrpChem,dShapes,pstatalign,fileparams);


if newsavefigs
    print(f1a,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_entropylocalshape'),'-dpng')
    savefig(f1a,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_entropylocalshape.fig'),'compact');
    print(f1b,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_entropyglobalshape'),'-dpng')
    savefig(f1b,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_entropyglobalshape.fig'),'compact');
end

end

%% SETTINGS alignment setup
if show_alignment == 1

do_meannorm = 1;
tbefore = cell(n_analysis,1);
tafter = cell(n_analysis,1);
idmax0 = cell(n_analysis,1);
idmax1 = cell(n_analysis,1);
lentbefore = zeros(n_analysis,1);
lentafter = zeros(n_analysis,1);
for i = 1:n_analysis
    %tbefore{i} = [1:(dTrapa(i)-1)]';
    %tafter{i} = [(dTrapa(i) + round((1/8)*dNframes(i))):5:(dNframes(i)-3)]';
    tbefore{i} = [dTrapa(i)-2];
    %tafter{i} = [(dTrapa(i) + 3):10:(dNframes(i)-3)]';
    tafter{i} = [dTrapa(i)+10]';
    lentbefore(i) = length(tbefore{i});
    lentafter(i) = length(tafter{i});
end
cnt1 = 180;

maxpre = max(lentbefore);
maxpost = max(lentafter);


%% alignment of pre-rapamycin signal

faidx = 1;
f0a = figure(faidx); 

tiledlayout(2,4,'TileSpacing','compact','Padding','compact');

nexttile(3,[2 2]); hold on;

kactin0_all=cell(n_analysis,1);
for i = 1:n_analysis
    kfactin = dKymoBin(i).actin;
    kactin_pre = kfactin(:,tbefore{i});
    if do_meannorm
        kactin_pre = kactin_pre./mean(kactin_pre);
    end

    %kactin_pre_align = nan(size(kactin_pre));
    kactin_pre_align = nan(360,maxpre);
    % align actin at each frame
    ncol = size(kactin_pre,2);   

    tmpid = nan(ncol,1); % store this for rotating the ActA signal

    for icol = 1:ncol
        
        [max1,tmpid(icol)] = max(kactin_pre(:,icol));
        kactin_pre_align(:,icol) = circshift(kactin_pre(:,icol),cnt1-tmpid(icol));

    end

    idmax0{i} = tmpid; % store this for rotating the ActA signal

    kactin0_all{i} = kactin_pre_align;

    p1 = plot(kactin_pre_align,'Color',[0.8500 0.3250 0.0980 0.5]);
%     p1.Color(1:3) = [0.8500 0.3250 0.0980]; 
%     p1.Color(4) = 0.5;
end
mxkactin0_all = cat(2,kactin0_all{:});
med_actin0 = median(mxkactin0_all,2);
%med_actin0 = mean(mxkactin0_all,2);
plot(med_actin0,'k','LineWidth',2)
%legend(dguvnames,'Location','bestoutside')
ylim([0 4]);
ax = gca; ax.FontSize = 14;
if do_meannorm
title('Actin intensity (mean-normalized) pre-rapamycin','FontSize',17)
else
title('Actin intensity pre-rapamycin','FontSize',17)
end

hold off



nexttile(1,[2 2]); hold on;


kActA0_all=cell(n_analysis,1);
for i = 1:n_analysis
    kfActA = dKymoBin(i).ActA;
    kActA_pre = kfActA(:,tbefore{i});
    if do_meannorm
        kActA_pre = kActA_pre./mean(kActA_pre);
    end

    %kActA_pre_align = nan(size(kActA_pre));
    kActA_pre_align = nan(360,maxpre);
    % align ActA at each frame following the actin alignment 
    ncol = size(kActA_pre,2);  
    
    tmpid2 = idmax0{i};
    for icol = 1:ncol
        kActA_pre_align(:,icol) = circshift(kActA_pre(:,icol),cnt1-tmpid2(icol));
    end


    kActA0_all{i} = kActA_pre_align;

    p1 = plot(kActA_pre_align,'Color',[0 0.4470 0.7410]);
%     p1.Color(1:3) = [0 0.4470 0.7410]; 
%     p1.Color(4) = 0.5;
end
mxkActA0_all = cat(2,kActA0_all{:});
med_ActA0 = median(mxkActA0_all,2);
%med_ActA0 = mean(mxkActA0_all,2);
plot(med_ActA0,'k','LineWidth',2)
%legend(dguvnames,'Location','bestoutside')
ylim([0 4]);
ax = gca; ax.FontSize = 14;
if do_meannorm
title('ActA intensity (mean-normalized) pre-rapamycin','FontSize',17)
else
title('ActA intensity pre-rapamycin','FontSize',17)
end


hold off



f0a.Position(3:4) = [1200,600];
drawnow;

if newsavefigs
    print(f0a,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_actin_align'),'-dpng')
    savefig(f0a,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_actin_align.fig'),'compact');
end





%% alignment of post rapamycin signals

fbidx = 2;
f0b = figure(fbidx); 

nbins = [1:360]';
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');



nexttile(3,[2 2]); hold on;


kactin1_all=cell(n_analysis,1);
for i = 1:n_analysis
    kfactin = dKymoBin(i).actin;
    kactin_post = kfactin(:,tafter{i});
    if do_meannorm
        kactin_post = kactin_post./mean(kactin_post);
    end
    %kactin_post_align = nan(size(kactin_post));
    kactin_post_align = nan(360,maxpost);
    % align actin at each frame
    ncol = size(kactin_post,2);   

    tmpid = nan(ncol,1); % store this for rotating the ActA signal

    for icol = 1:ncol
        
        [max1,tmpid(icol)] = max(kactin_post(:,icol));
        kactin_post_align(:,icol) = circshift(kactin_post(:,icol),cnt1-tmpid(icol));

    end

    idmax1{i} = tmpid; % store this for rotating the ActA signal

    kactin1_all{i} = kactin_post_align;

    p1 = plot(kactin_post_align,'Color',[0.8500 0.3250 0.0980 0.5]);
%     p1.Color(1:3) = [0.8500 0.3250 0.0980]; 
%     p1.Color(4) = 0.5;
end

% postt_kactin1_all = cell(maxpost,1);
% med_postt_actin1 = cell(maxpost,1);
% for j = 1:maxpost
%     mx_tmptime = nan(360,18);
%     for i = 1:n_analysis
%     guvkactin = kactin1_all{i};
%     guvkactin_time = guvkactin(:,j);
%     mx_tmptime(:,i) = guvkactin_time;
% 
%     end
%     postt_kactin1_all{j} = mx_tmptime;
%     med_postt_actin1{j} = median(mx_tmptime,2,'omitnan');
%     plot(med_postt_actin1{j},'k','LineWidth',2)
% end
mxkactin1_all = cat(2,kactin1_all{:});
med_actin1 = median(mxkactin1_all,2);
%med_actin1 = mean(mxkactin1_all,2);
plot(med_actin1,'k','LineWidth',2)

%legend(dguvnames,'Location','bestoutside')



ylim([0 4]);
ax = gca; ax.FontSize = 14;
if do_meannorm
title('Actin intensity (mean-normalized) post-rapamycin','FontSize',17)
else
    title('Actin intensity post-rapamycin','FontSize',17)

end

hold off






nexttile(1,[2 2]); hold on;

kActA1_all=cell(n_analysis,1);
for i = 1:n_analysis
    kfActA = dKymoBin(i).ActA;
    kActA_post = kfActA(:,tafter{i});    
    if do_meannorm
        kActA_post = kActA_post./mean(kActA_post);
    end
    %kActA_post_align = nan(size(kActA_post));
    kActA_post_align = nan(360,maxpost);
    % align ActA at each frame following the actin alignment 
    ncol = size(kActA_post,2);  
    
    tmpid2 = idmax1{i};
    for icol = 1:ncol
        kActA_post_align(:,icol) = circshift(kActA_post(:,icol),cnt1-tmpid2(icol));
    end


    kActA1_all{i} = kActA_post_align;

    p1 = plot(kActA_post_align,'Color',[0 0.4470 0.7410 0.5]);
%     p1.Color(1:3) = [0 0.4470 0.7410]; 
%     p1.Color(4) = 0.5;
end
mxkActA1_all = cat(2,kActA1_all{:});
med_ActA1 = median(mxkActA1_all,2);
%med_ActA1 = mean(mxkActA1_all,2);
plot(med_ActA1,'k','LineWidth',2)
%legend(dguvnames,'Location','bestoutside')
ylim([0 4]);
ax = gca; ax.FontSize = 14;
if do_meannorm
title('ActA intensity (mean-normalized) post-rapamycin','FontSize',17)
else
title('ActA intensity post-rapamycin','FontSize',17)

end
hold off



f0b.Position(3:4) = [1200,600];
drawnow;

if newsavefigs
    print(f0b,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_fig_ActA_align'),'-dpng')
    savefig(f0b,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_ActA_align.fig'),'compact');
end




end


%% SETTINGS spatial correlation
if show_spatialcorr

corrtypeSp = 'Pearson';
%corrtypeSp = 'Spearman';
if do_fixedspcorr
    dKtmp = dKymoBin;
    kt = 'fix';
else
    dKtmp = dKymoVar;
    kt = 'var';
end

%% spatial ActA correlation

%%%%%% CHOOSE PROPER CHANNEL HERE
refidx = 2; % memb = 1, ActA = 2, actin = 3
figidx = 10;

[alignR_ActA, alignC_ActA,alignRR_ActA,alignCC_ActA,spcorrparams_ActA] =...
    spatialcorr(refidx,dParams,dKtmp,corrtypeSp,pstatalign,fileparams);
[spc_ActA,Rmed_ActA]= plotspatialcorr(figidx,alignR_ActA,alignC_ActA,alignRR_ActA,alignCC_ActA,...
    spcorrparams_ActA,pstatalign,fileparams);

%actamemb = cat(2,alignR_ActA{:,1});

if newsavefigs
print(spc_ActA(1),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_memb_SPcorr_',corrtypeSp(1:4)),'-dpng')
print(spc_ActA(2),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_actin_SPcorr_',corrtypeSp(1:4)),'-dpng')
print(spc_ActA(3),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_vel_SPcorr_',corrtypeSp(1:4)),'-dpng')

if do_absvals
print(spc_ActA(4),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_ABSdcumul_SPcorr_',corrtypeSp(1:4)),'-dpng')
print(spc_ActA(5),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_ABScurv_SPcorr_',corrtypeSp(1:4)),'-dpng')
else
print(spc_ActA(4),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_dcumul_SPcorr_',corrtypeSp(1:4)),'-dpng')
print(spc_ActA(5),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_curv_SPcorr_',corrtypeSp(1:4)),'-dpng')
end

print(spc_ActA(6),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_npointsused_SPcorr_',corrtypeSp(1:4)),'-dpng')
savefig(spc_ActA,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_ActA_spatialcorrs_',corrtypeSp(1:4),'.fig'),'compact');
end

% 
%% spatial actin correlation

%%%%%% CHOOSE PROPER CHANNEL HERE
refidx = 3; % memb = 1, ActA = 2, actin = 3
figidx = 20; % number to start new set of figs

[alignR_actin, alignC_actin,alignRR_actin,alignCC_actin,spcorrparams_actin] = ...
    spatialcorr(refidx,dParams,dKtmp,corrtypeSp,pstatalign,fileparams);
[spc_actin,Rmed_actin] = plotspatialcorr(figidx,alignR_actin,alignC_actin,alignRR_actin,alignCC_actin,...
    spcorrparams_actin,pstatalign,fileparams);

if newsavefigs
print(spc_actin(1),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_memb_SPcorr_',corrtypeSp(1:4)),'-dpng')
print(spc_actin(2),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_ActA_SPcorr_',corrtypeSp(1:4)),'-dpng')
print(spc_actin(3),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_vel_SPcorr_',corrtypeSp(1:4)),'-dpng')

if do_absvals
print(spc_actin(4),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_ABSdcumul_SPcorr_',corrtypeSp(1:4)),'-dpng')
print(spc_actin(5),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_ABScurv_SPcorr_',corrtypeSp(1:4)),'-dpng')
else
print(spc_actin(4),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_dcumul_SPcorr_',corrtypeSp(1:4)),'-dpng')
print(spc_actin(5),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_curv_SPcorr_',corrtypeSp(1:4)),'-dpng')
end

print(spc_actin(6),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_npointsused_SPcorr_',corrtypeSp(1:4)),'-dpng')
savefig(spc_actin,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_actin_spatialcorrs_',corrtypeSp(1:4),'.fig'),'compact');
end

end

if show_signalcorr
%% SETTINGS signal correlation

corrtypeSig = 'Pearson';
%corrtypeSig = 'Spearman';
if do_fixedsigcorr
    dKtmp = dKymoBin;
    kt = 'fix';
else
    dKtmp = dKymoVar;
    kt = 'var';
end

%% all signal ActA correlation

%%%%%% CHOOSE PROPER CHANNEL HERE
refidx = 2; % memb = 1, ActA = 2, actin = 3
figidx = 30; % number to start new set of figs

[R1_ActA,C1_ActA,pval_ActA,refAgg_ActA,compAgg_ActA,scp_ActA] = ...
    signalcorr(refidx,dParams,dKtmp,corrtypeSig,pstatalign,fileparams);

allc_ActA = ...
    plotsignalcorr(figidx,R1_ActA,C1_ActA,refAgg_ActA,compAgg_ActA,scp_ActA,pstatalign,fileparams);

if newsavefigs
print(allc_ActA(1),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_memb_allcorr',corrtypeSig(1:4)),'-dpng')
print(allc_ActA(2),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_actin_allcorr',corrtypeSig(1:4)),'-dpng')
print(allc_ActA(3),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_vel_allcorr',corrtypeSig(1:4)),'-dpng')

if do_absvals
print(allc_ActA(4),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_ABSdcumul_allcorr',corrtypeSig(1:4)),'-dpng')
print(allc_ActA(5),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_ABScurv_allcorr',corrtypeSig(1:4)),'-dpng')
else
print(allc_ActA(4),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_dcumul_allcorr',corrtypeSig(1:4)),'-dpng')
print(allc_ActA(5),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_curv_allcorr',corrtypeSig(1:4)),'-dpng')
end

print(allc_ActA(6),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_ActA_npointsused_allcorr_',corrtypeSig(1:4)),'-dpng')
savefig(allc_ActA,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_ActA_allcorr_',corrtypeSig(1:4),'.fig'),'compact');

end

%% all signal actin correlation
% 
% %%%%%% CHOOSE PROPER CHANNEL HERE
% refidx = 3; % memb = 1, ActA = 2, actin = 3
% figidx = 40; % number to start new set of figs
% 
% [R1_actin,C1_actin,pval_actin,refAgg_actin,compAgg_actin,scp_actin] = ...
%     signalcorr(refidx,dParams,dKtmp,corrtypeSig,pstatalign,fileparams);
% 
% allc_actin = ...
%     plotsignalcorr(figidx,R1_actin,C1_actin,refAgg_actin,compAgg_actin,scp_actin,pstatalign,fileparams);
% 
% if newsavefigs
% print(allc_actin(1),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_memb_allcorr',corrtypeSig(1:4)),'-dpng')
% print(allc_actin(2),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_ActA_allcorr',corrtypeSig(1:4)),'-dpng')
% print(allc_actin(3),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_vel_allcorr',corrtypeSig(1:4)),'-dpng')
% 
% if do_absvals
% print(allc_actin(4),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_ABSdcumul_allcorr',corrtypeSig(1:4)),'-dpng')
% print(allc_actin(5),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_ABScurv_allcorr',corrtypeSig(1:4)),'-dpng')
% else
% print(allc_actin(4),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_dcumul_allcorr',corrtypeSig(1:4)),'-dpng')
% print(allc_actin(5),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_curv_allcorr',corrtypeSig(1:4)),'-dpng')
% end
% 
% print(allc_actin(6),strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_fig_actin_npointsused_allcorr_',corrtypeSig(1:4)),'-dpng')
% savefig(allc_actin,strcat(newoutputfolder,outputname,'_',tstamp_analysis,'_',kt,'_actin_allcorr_',corrtypeSig(1:4),'.fig'),'compact');
% 
% end

end
%% spatiotemporal patches
% corrtype = 'Spearman';
% corrtype_sptmp = 'Pearson';
% nbins = 360;
% patches = genpatch(nbins,0.5);
% plotspatiotempcorr(dParams,dKymoBin,patches,corrtype_sptmp,fileparams);


%% save m file
if newsavefigs

    % save mfile if saving figs 
    currentfile = strcat(mfilename,'.m');
    destinationfile = strcat(newoutputfolder,mfilename,'_',tstamp_analysis,'.m');
    copyfile(currentfile,destinationfile);   
end


%---------------------------------------------------------------

%---------------------------------------------------------------
function bdy_sg = sgsmoothboundary(I,bdy,windowWidth,polynomialOrder)
%% smooth boundary using savitsky golay filter

% make equal spacing before sgsmooth to have each 
n = size(bdy,1); % circular
Borig1 = interpboundary(bdy,n); % circular


% get unique points, non-circular
x = Borig1(1:end-1,2); y = Borig1(1:end-1,1);

% overlap curve of unique boundary points
% workaround edge effects problems with sgolayfilt
overlap = 2*windowWidth;
xplus = [x; x(1:overlap)];
yplus = [y; y(1:overlap)];

% Now smooth with a Savitzky-Golay sliding polynomial filter
xs = sgolayfilt(xplus, polynomialOrder, windowWidth);
ys = sgolayfilt(yplus, polynomialOrder, windowWidth);
%xs = [xs; xs(1)]; ys = [ys; ys(1)];

Borigs1 = [ys,xs];

% keep region with good smoothing, and keep length same as unique boundary points
Borigs2 = Borigs1(windowWidth+1:end-windowWidth,:);

% circular shift to match indices to original list of boundary points
Borigs3 = circshift(Borigs2,windowWidth);

% make circular bc currently NONE of the lists are circular
bdy_sg = [Borigs3;Borigs3(1,:)];

%Borigs4 = interpboundary(Borigs3c,length(Borig));
Borigs4 = bdy_sg;

% f50 = figure(50); %clf(f50);
% imshow(I);
% axis equal, hold on;
% plot(bdy(1,2),bdy(1,1),'r*','LineWidth',10)
% plot(bdy(end,2),bdy(end,1),'wx','LineWidth',10)
% plot(bdy(:,2),bdy(:,1),'g*','LineWidth',2)
% plot(Borigs4(1,2),Borigs4(1,1),'r*','LineWidth',10)
% plot(Borigs4(end,2),Borigs4(end,1),'wx','LineWidth',10)
% plot(Borigs4(:,2),Borigs4(:,1),'b*','LineWidth',2)
% title('Original (green), Savitzky-Golay smoothed (blue)',...
%     'FontSize',10);
% f50.Position(3:4) = [500 500];
% drawnow

% % original boundary
% plot(bdy(1,2),bdy(1,1),'r*','LineWidth',10)
% plot(bdy(end,2),bdy(end,1),'wx','LineWidth',10)
% plot(bdy(:,2),bdy(:,1),'g*','LineWidth',2)
% % interpolated boundary for equal spacing
% plot(Borig1(1,2),Borig1(1,1),'r*','LineWidth',10)
% plot(Borig1(end,2),Borig1(end,1),'wx','LineWidth',10)
% plot(Borig1(:,2),Borig1(:,1),'g*','LineWidth',2)
% title('Original boundary (green), equal spaced boundary (blue)')
% hold off
% drawnow


end

function varargout = interpboundary(bdy,n,varargin)
    % This function takes a boundary and, possibly, data and a number n and 
    % interpolates the boundary and any data so that it is now of length n
    % for circular data: n = (# of desired unique points) + 1
    % output is circular if input is circular
    itype = 'makima'; % avoids excessive undulations
    %itype = 'pchip';
    if nargin~=nargout+1
        disp('Error: Number of input and output data sets must equal')
    else 
        d       = diff(bdy); % difference between rows of bdy
        dist    = sqrt(d(:,1).^2+d(:,2).^2);
        cumdist = [0;cumsum(dist)]; % (x) cumulative distances around boundary
        perim   = max(cumdist);
        s       = linspace(0,1,n)*perim; % (xq) n query points from 0 to perim
        % vq = interp1(x,v,xq,method)
        % interpolate (row,col) ie (y,x) coordinates across cumulative sum of boundary arc length
        row       = interp1(cumdist,bdy(:,1),s,itype)';
        col       = interp1(cumdist,bdy(:,2),s,itype)';
        B = [row col];

        varargout = cell(length(varargin)+1,1);
        varargout{1}=B;
        for i=1:length(varargin)
            % if nans just return nans
            if all(isnan(varargin{i}))
                varargout{i+1} = nan(length(row),1);
%                 disp('Error: Found NaNs in list')
            else
%                 interpolate fluorescent signal,velocities,other quantities
%                 across cumulative sum of boundary arc length
                varargout{i+1}  = interp1(cumdist,varargin{i},s,itype)';
            end

        end
     end
end

function varargout = mavgsmoothboundary(mfilter,varargin)
    %%
    % Moving average; Should be an odd number
    s = mfilter;


        n = length(varargin); % # boundary quantities
        varargout = cell(1,n);
    if s>0 % then do movavg
        for k = 1:n % each boundary quantity
            bdy = varargin{k};
            sbdy = bdy;
            for i=1:s % 
                sbdy = sbdy+circshift(bdy,i)+circshift(bdy,-i);
            end
            sbdy = sbdy./(2*s+1);
            varargout{k} = sbdy;
        end    
        % figure()
        %plot(bdy(1:200,1),bdy(1:200,2),'x',sbdy(1:200,1),sbdy(1:200,2),'.r')
    else
        varargout = varargin;
    end
end

function [patches] = genpatch(nbins,prop)
   
if prop > 0.5
    
    patches = 1:nbins;
   
elseif prop == 0.5
    
    binbounds = floor(0.5*prop*nbins);
    nearpoint = ceil(nbins/2);
    farpoint = nbins;
   
    
    l0 = 1:nbins;
    p0 = [nearpoint,farpoint];
    tmp = cell(1,2); % bin indices of each patch
    npts = cell(1,2); % number of points in each patch
    
    for g = 1:2
        % bring p0(g) index to front of list of bins
        l1 = circshift( l0, -(p0(g)-1) );
        % bring binbounds to the left of it
        l2 = circshift( l1, binbounds);
        % select patch with 2*binbounds+1 points
        l3 = l2(1:2*binbounds+1); % must be +1
    
        % regions to consider
        tmp{g} = l3';
    
        % number of perimeter points in region
        npts{g} = length(l3');
    end
    
    c = mod(nbins,4);


    switch c
        case 0
            patches{1} = tmp{1}(1:end-1);
            patches{2} = tmp{2}(1:end-1);
        case 1
            patches{1} = tmp{1}(1:end-1);
            patches{2} = tmp{2}(1:end);
        case 2
            patches{1} = tmp{1}(1:end);
            patches{2} = tmp{2}(1:end);
        case 3
            patches{1} = tmp{1};
            patches{2} = [tmp{2};tmp{2}(end)+1];    
    end

else
    binbounds = floor(0.5*prop*nbins);
    nearpoint = ceil(nbins/2);
    farpoint = nbins;
   
    
    l0 = 1:nbins;
    p0 = [nearpoint,farpoint];
    patches = cell(1,2); % bin indices of each patch
    npts = cell(1,2); % number of points in each patch
    
    for g = 1:2
        % bring p0(g) index to front of list of bins
        l1 = circshift( l0, -(p0(g)-1) );
        % bring binbounds to the left of it
        l2 = circshift( l1, binbounds);
        % select patch with 2*binbounds+1 points
        l3 = l2(1:2*binbounds+1); % must be +1
    
        % regions to consider
        patches{g} = l3';
    
        % number of perimeter points in region
        npts{g} = length(l3');
    end

end

end

function dVariance = variancecomp(dP,dK,pranges,fileparams)
% (dParams,dKymoBin,varshapes,pranges,fileparams);

v2struct(fileparams);
v2struct(pranges);


% recall: pranges = v2struct(rmemb,rActA,ractin,rvel,rdc,rcurv);

% need fixed bins to compare entropy so only do when kymofixed


% kymos or bdyquantities to consider
nky = 5;

% singke 
% compile entropy (guv #, bdyquant #)
dVariance = cell(n_analysis,nky); 
% each element is a column vector of entropy time series (1 x nframes)

for k = 1:n_analysis  % for each GUV

    nbins = dP{k}.nbins;
    nframes = dP{k}.nframes;
    
    tmpKymo = dK(k);
    
    % 4 kymos to consider
    tmpkActA  = tmpKymo.ActA; % kymos for fluorescent channels
    tmpkactin = tmpKymo.actin;
    tmpkdc = tmpKymo.dc; % kymo cumulative disp
    tmpkcurv  = tmpKymo.curv; % kymo curvatures

    tmpkmemb = tmpKymo.memb;
    
    % make cumulative displacement first frame all nans
    % tmpkdc(:,1) = nan;
    
    % for variance - always use absolute values of dc and curv
    % compute corr to absolute values - NOT MUCH DIFFERENT
    %if do_absvals
        tmpkdc = abs(tmpkdc);
        tmpkcurv = abs(tmpkcurv);
    %end
    
    %% temporal region to consider
    tmpframes = 1:nframes;

%% signal variance


     % spatial regions to consider
    tmppatch = 1:nbins;


    % number of perimeter points in region
    npts = length(tmppatch);


    ky = cell(1,nky);
    varky = cell(1,nky);

    ky{1} = tmpkActA(tmppatch,tmpframes); % ActA
    ky{2}  = tmpkactin(tmppatch,tmpframes); % actin
    ky{3} = tmpkdc(tmppatch,tmpframes); % cumulative displacemennt    
    ky{4} = tmpkcurv(tmppatch,tmpframes); % curvature
    ky{5} = tmpkmemb(tmppatch,tmpframes);

    tmpmemb = tmpkmemb(tmppatch,tmpframes); % membrane marker

    % compute edge of bins for signal discretizations for each kymograph
    uedges = [rActA(2),ractin(2),rdc(2),rcurv(2),rmemb(2)];
    ledges = [rActA(1),ractin(1),rdc(1),rcurv(1),rmemb(1)];



    umemb = rmemb(2); % upper range of memb marker
    lmemb = rmemb(1); % lower range of memb marker

    nbins = 256; % # of bins
    %edges = -u: 2*u/nbins: u;
    edges = cell(1,nky);
    for id = 1:nky
        edges{id} = ledges(id): (uedges(id) - ledges(id))/nbins: uedges(id);
    end

   % edges for membrane marker
    medge = lmemb : (umemb - lmemb)/nbins : umemb;





    % compute signal entropy 

    for i = 1:nky %ith boundary quantity
        % initialize vector of entropy (scalar for each frame)
        varky{i} = nan(length(tmpframes),1); 


        
        % edges for ith kymograph
        tmpedge = edges{i};
        % to allow for bins to include all values below 0 potentially in biochem kymos
        tmpedge(1) = -100; 
        tmpedge(end) = 1000;

        medge(1) = -100; 
        medge(end) = 1000;

        % all values in ith kymograph
        tmpky = ky{i};
        jj = 0;
        for j = tmpframes  % this is a vector from 1:nframes usually

            % if tmpframes has skipped frames, keep jj
            jj = jj + 1;
            % extract column vector of frame j from kymograph i
            tmpdata = tmpky(:,j);
            % extract column vector of frame j from memb kymograph
             tmpmb = tmpmemb(:,j);
            
%             % some spatial bins in a given frame of a fixed kymo may be nans
%             % these will have less than 301 points so ignore to ensure proper comparison
%             % want fixed signal and spatial discretization
%             % e.g. guv so small that # points < 301 spatial bins 
             if ~any(isnan(tmpdata)) 
%                 % compute entropy of frame j in kymograph i
% 
%                 a1 = tmpdata;
%                 na1 = histcounts(a1,tmpedge)'; 
%                 Y = discretize(a1,tmpedge);    
%                 pixelCount = na1;
%                 grayLevels = [0:255]';
% 
%                 a2 = tmpmb;
%                 na2 = histcounts(a2,medge)'; 
%                 Y2 = discretize(a2,medge);    
%                 pixelCount2 = na2;
%                 grayLevels2 = [0:255]';
% 
%                 if any(isnan(Y)) == 1 % some data points are outside the min/max edges
%                       disp('Some points are outside min/max edges of histogram')
% 
%                 else
%                     a1i = uint8(grayLevels(Y));
%                     %a1i = im2uint8(grayLevels(Y));
%                     a1h = na1;
%                     a1h(a1h==0) = [];
%                     a1p = a1h./numel(a1);
%                     a1e = -sum(a1p.*log2(a1p));
%                     a1e_mat = entropy(a1i);
% 
% 
%                     a2i = uint8(grayLevels(Y2));
%                     a2e_mat = entropy(a2i);


                    %tmpent = a1e/a2e_mat;
                    %tmpent = a1e;

                    tmpvar = var(tmpdata,'omitnan');
                    tmpmean = mean(tmpdata,'omitnan');


                    % store entropy
                    varky{i}(jj) = sqrt(tmpvar)/tmpmean;
    
    %                 f1 = figure(); 
    %                 ax1 = imagesc(a1i); yline(0.5:length(a1)-0.5);
    %                 colormap(parula(256))
    %                 ax1.Parent.CLim = [0 255];
    %                 cc = colorbar('FontSize',fs);
    %                 cc.Label.String = ('8-bit pixel intensity');
    %                 ylabel('Spatial Points','FontSize',fs)
    %                 yticks(1:length(a1i))
    %                 yticklabels(1:length(a1i))
    %                 f1.Position(3) = 200;
    %                 
    %                 figure;
    %                 imhist(a1i,parula);grid on
    %                 ylim([0,5]);
    %                 
    %                 figure;
    %                 [pixelCount, grayLevels] = imhist(a1i,parula);
    %                 bar(grayLevels, pixelCount); % Plot it as a bar chart.
    %                 grid on;
    %                 ylim([0,5])
    %                 title('Histogram of pixel values', 'FontSize', fs);
    %                 xlabel('Re-scaled pixel values','FontSize',fs);
    %                 ylabel('Counts','FontSize',fs);
    %                 
    %                 figure;
    %                 bar(grayLevels, pixelCount/sum(pixelCount),BaseValue=3e-3); % Plot it as a bar chart.
    %                 set (gca,'yscale','log');
    %                 ylim([3e-3,1])
    %                 grid on;
    %                 title('Prob of pixel values', 'FontSize', fs);
    %                 xlabel('Re-scaled pixel values','FontSize',fs);
    %                 ylabel('Probability','FontSize',fs);
                    

%                end
               
            else

                disp('There are some nans in the kymo')
                disp('guv #:'); k
                disp('kymo #:'); i
                disp('frame #:'); j
   
            end

        end % jth frame

    end % ith kymo

    dVariance(k,:) = varky;

end % end of kth guv

end


function f1 = plotvariancelocalshape(fignum,dVariance,pstatalign,fileparams)
v2struct(fileparams);
v2struct(pstatalign);


tlsize = 14;
dcolors = [ 0      0.4470 0.7410; 
            0.8500 0.3250 0.0980;
            0.4660 0.6740 0.1880;
            0.4940 0.1840 0.5560;
            0.3010 0.7450 0.9330];

% [maxf,idmaxf] = max(dNframes); %GUV with most frames
% tickstep = df{idmaxf}.dataParams.tickstep;
% [maxt,~] = max(dTrapa); %GUV with latest time after rapa is on
% maxd = max(maxt-dTrapa);
% maxf1 = maxf + maxd;

tmpframes = [1:maxf1]';
tmpframes = tmpframes - maxt;

transp = 0.5; %transparency
nky = size(dVariance,2); % # of kymographs for which entropy is computed

%    f1 = figure(fignum); clf(f1);
f1 = figure(fignum); 
hold on

tmpVallguv = nan(maxf1,n_analysis);
Kmed = cell(nky,1);


for i=1:nky % entropy for each kymo
    
    varky = dVariance(:,i);%entropy of ith kymo,all GUVs

    for k = 1:n_analysis % each guv


        % align time series of entropy based on largest t_rapa
    
        tmpVky = nan(maxf1,1);
        ell  = maxt-dTrapa(k);


        tmpVky(ell+1:ell+dNframes(k)) = varky{k}; %kth GUV 
        tmpVallguv(:,k) = tmpVky;

        ax7 = plot(tmpframes,tmpVky,'-','LineWidth',0.5,'Color',[dcolors(i,:),transp]);
   

    end
    Kmed{i} = median(tmpVallguv,2,'omitnan');
    %ax8 = plot(Kmed,'LineWidth',2,'Color',[dcolors(i,:),1]);

end



ax8 = gobjects(nky,1);
for i = 1:nky
    ax8(i) = plot(tmpframes,Kmed{i},'LineWidth',2,'Color',[dcolors(i,:),1]);
end
xline(0,'LineWidth',2);

if do_absvals
legend(ax8(1:nky),["ActA";"Actin";"Abs.Val. Cumulative disp.";"Abs.Val. Curvature"],'Location','southeast',...
    'FontSize',tlsize,'Location','eastoutside')
else
legend(ax8(1:nky),["ActA";"Actin";"Cumulative disp.";"Curvature"],'Location','southeast',...
    'FontSize',tlsize,'Location','eastoutside')
end
%     legend('Velocity','Curvature','Cumul. displacement','Membrane marker','Actin','ActA',...
%         'Location','southeast')

% legend('ActA','Actin','Cumulative disp.','Curvature','Location','southeast')
%tmpframes = tmpframes - maxt;
xlim([tmpframes(1)-1,tmpframes(end)+1]); 
xticks([tmpframes(1):tickstep:tmpframes(end)])

%ylim([0,8])
ax7.Parent.YLabel.String = 'Variance';
ax7.Parent.YLabel.FontName = 'Arial';
%ax7.Parent.YLabel.FontSize = tlsize;
ax7.Parent.XLabel.String = 'Frames';
ax7.Parent.XLabel.FontName = 'Arial';
%ax7.Parent.XLabel.FontSize = tlsize;
% title(strcat(outputname,{' '},': entropy in boundary quantities'),'Interpreter','none',...
%     'FontSize',tlsize)

ax7.Parent.FontSize = tlsize;
title({outputname;'Variance in boundary quantities'},'Interpreter','none',...
    'FontSize',tlsize+3)
 
hold off
f1.Position(3:4) = [800,500];

drawnow;
%     if newsavefigs
%         print(f1,strcat(newoutputfolder,guvname,'_',tstamp1,'_fig_globalentropy'),'-dpng')
%     end



end % end plotvariancelocalshape


function f1 = plotvarianceglobalshape(fignum,dVChem,dShapes,pstatalign,fileparams)
v2struct(fileparams);
v2struct(pstatalign);
% v2struct(varshapes); % global shape features
% % recall: varshapes = v2struct(dEccent,dAspratio_ml,dAspratio_ind,dAreas);


nshapes = size(dShapes,2);


tlsize = 14;
dcolors = [0 0.4470 0.7410;  % ActA entropy
        0.8500 0.3250 0.0980; % actin entropy

        0.4660 0.6740 0.1880; % Eccent
        0.4940 0.1840 0.5560; % Aspratio_ml
        0.3010 0.7450 0.9330; % Aspration_ind
        0.6350 0.0780 0.1840]; % Areas


% [maxf,idmaxf] = max(dNframes); %GUV with most frames
% tickstep = df{idmaxf}.dataParams.tickstep;
% [maxt,~] = max(dTrapa); %GUV with latest time after rapa is on
% maxd = max(maxt-dTrapa);
% maxf1 = maxf + maxd;

tmpframes = [1:maxf1]';
tmpframes = tmpframes - maxt;

transp = 0.5; %transparency
nky = size(dVChem,2); % # of kymographs for which entropy is computed

%    f1 = figure(fignum); clf(f1);
f1 = figure(fignum); 
hold on


Vallguv = cell(nky,1);
Kmed = cell(nky,1);
ax8 = gobjects(nky+nshapes,1);

for i=1:nky % entropy for each kymo
    
    varky = dVChem(:,i);%entropy of ith kymo,all GUVs
    tmpVallguv = nan(maxf1,n_analysis);

    for k = 1:n_analysis % each guv


        % align time series of entropy based on largest t_rapa
    
        tmpEky = nan(maxf1,1);
        ell  = maxt-dTrapa(k);


        tmpEky(ell+1:ell+dNframes(k)) = varky{k}; %kth GUV 
        tmpVallguv(:,k) = tmpEky;

        ax7 = plot(tmpframes,tmpEky,'-','LineWidth',0.5,'Color',[dcolors(i,:),transp]);
   

    end
    Vallguv{i} = tmpVallguv;
    Kmed{i} = median(tmpVallguv,2,'omitnan');
    %ax8 = plot(Kmed,'LineWidth',2,'Color',[dcolors(i,:),1]);

end


for i = 1:nky
    ax8(i) = plot(tmpframes,Kmed{i},'LineWidth',3,'Color',[dcolors(i,:),1]);
end
xline(0,'LineWidth',2);

% legend(ax8(1:nky),["ActA";"Actin"],'FontSize',tlsize,'Location','eastoutside')

%ylim([0,8])
ax7.Parent.YLabel.String = 'Variance';
ax7.Parent.YLabel.FontName = 'Arial';

%% now plot shape parameters on right y-axis

yyaxis right

Sallguv = cell(nshapes,1);
KmedS = cell(nshapes,1);

for i = 1:nshapes

    shapeID= dShapes(:,i); % ith shape parameter of all GUVs
    tmpSallguv = nan(maxf1,n_analysis);


    for k = 1:n_analysis
        % align time series of entropy based on largest t_rapa
    
        tmpshape = nan(maxf1,1);
        ell  = maxt-dTrapa(k);


        tmpshape(ell+1:ell+dNframes(k)) = shapeID{k}; %kth GUV 
        tmpSallguv(:,k) = tmpshape;

        ax7 = plot(tmpframes,tmpshape,'-','LineWidth',0.5,'Color',[dcolors(i+nky,:),transp]);
   

    end

    Sallguv{i} = tmpSallguv;
    KmedS{i} = median(tmpSallguv,2,'omitnan');
    %ax8 = plot(Kmed,'LineWidth',2,'Color',[dcolors(i,:),1]);


end

% ax9 = gobjects(nshapes,1);
for i = 1:nshapes
    ax8(i+nky) = plot(tmpframes,KmedS{i},'-','LineWidth',3,'Color',[dcolors(i+nky,:),1]);
end


legend(ax8(1:end),["ActA";"Actin";"Eccent.";"AspRatio-Mat";"AspectRatio-Ind";"Norm. Area"],...
    'FontSize',tlsize,'Location','eastoutside')

% legend(ax8(nky+1:end),["Eccent.";"AspRatio-Mat";"AspectRatio-Ind";"Norm. Area"],...
%     'FontSize',tlsize,'Location','eastoutside')

% if do_absvals
% legend(ax8(1:nky),["ActA";"Actin";"Abs.Val. Cumulative disp.";"Abs.Val. Curvature"],'Location','southeast',...
%     'FontSize',tlsize,'Location','eastoutside')
% else
% legend(ax8(1:nky),["ActA";"Actin";"Cumulative disp.";"Curvature"],'Location','southeast',...
%     'FontSize',tlsize,'Location','eastoutside')
% end
%     legend('Velocity','Curvature','Cumul. displacement','Membrane marker','Actin','ActA',...
%         'Location','southeast')

% legend('ActA','Actin','Cumulative disp.','Curvature','Location','southeast')
%tmpframes = tmpframes - maxt;
xlim([tmpframes(1)-1,tmpframes(end)+1]); 
xticks([tmpframes(1):tickstep:tmpframes(end)])

ylim([0,0.85])
ax7.Parent.YLabel.String = 'Shape quantities';
ax7.Parent.YLabel.FontName = 'Arial';
%ax7.Parent.YLabel.FontSize = tlsize;
ax7.Parent.XLabel.String = 'Frames';
ax7.Parent.XLabel.FontName = 'Arial';
%ax7.Parent.XLabel.FontSize = tlsize;
% title(strcat(outputname,{' '},': entropy in boundary quantities'),'Interpreter','none',...
%     'FontSize',tlsize)

ax7.Parent.FontSize = tlsize;
title({outputname;'Variance in boundary quantities'},'Interpreter','none',...
    'FontSize',tlsize+3)
 
hold off
f1.Position(3:4) = [800,500];

drawnow;
%     if newsavefigs
%         print(f1,strcat(newoutputfolder,guvname,'_',tstamp1,'_fig_globalentropy'),'-dpng')
%     end



end % end plotvarianceglobalshape


function dEntrp = entropycomp(dP,dK,pranges,fileparams)
% signalentropy(dParams,dKymoBin,varshapes,pranges,fileparams);

v2struct(fileparams);
v2struct(pranges);
mthresh = nan(n_analysis,1);

% recall: pranges = v2struct(rmemb,rActA,ractin,rvel,rdc,rcurv);

% need fixed bins to compare entropy so only do when kymofixed


% kymos or bdyquantities to consider
nky = 5;

% singke 
% compile entropy (guv #, bdyquant #)
dEntrp = cell(n_analysis,nky); 
% each element is a column vector of entropy time series (1 x nframes)

for k = 1:n_analysis  % for each GUV

    nbins = dP{k}.nbins;
    nframes = dP{k}.nframes;
    
    tmpKymo = dK(k);
    
    % 4 kymos to consider
    tmpkActA  = tmpKymo.ActA; % kymos for fluorescent channels
    tmpkactin = tmpKymo.actin;
    tmpkdc = tmpKymo.dc; % kymo cumulative disp
    tmpkcurv  = tmpKymo.curv; % kymo curvatures

    tmpkmemb = tmpKymo.memb;
    
    % make cumulative displacement first frame all nans
    % tmpkdc(:,1) = nan;
    
    
    % compute corr to absolute values - NOT MUCH DIFFERENT
    if do_absvals
        tmpkdc = abs(tmpkdc);
        tmpkcurv = abs(tmpkcurv);
    end
    
    %% temporal region to consider
    tmpframes = 1:nframes;

%% signal entropy


     % spatial regions to consider
    tmppatch = 1:nbins;


    % number of perimeter points in region
    npts = length(tmppatch);


    ky = cell(1,nky);
    entky = cell(1,nky);

    ky{1} = tmpkActA(tmppatch,tmpframes); % ActA
    ky{2}  = tmpkactin(tmppatch,tmpframes); % actin
    ky{3} = tmpkdc(tmppatch,tmpframes); % cumulative displacemennt    
    ky{4} = tmpkcurv(tmppatch,tmpframes); % curvature
    ky{5} = tmpkmemb(tmppatch,tmpframes);

    tmpmemb = tmpkmemb(tmppatch,tmpframes); % membrane marker

    % compute edge of bins for signal discretizations for each kymograph
    uedges = [rActA(2),ractin(2),rdc(2),rcurv(2),rmemb(2)];
    ledges = [rActA(1),ractin(1),rdc(1),rcurv(1),rmemb(1)];



    umemb = rmemb(2); % upper range of memb marker
    lmemb = rmemb(1); % lower range of memb marker

    nbins = 256; % # of bins
    %edges = -u: 2*u/nbins: u;
    edges = cell(1,nky);
    for id = 1:nky
        edges{id} = ledges(id): (uedges(id) - ledges(id))/nbins: uedges(id);
    end

   % edges for membrane marker
    medge = lmemb : (umemb - lmemb)/nbins : umemb;





    % compute signal entropy 

    for i = 1:nky %ith boundary quantity
        % initialize vector of entropy (scalar for each frame)
        entky{i} = nan(length(tmpframes),1); 


        
        % edges for ith kymograph
        tmpedge = edges{i};
        % to allow for bins to include all values below 0 potentially in biochem kymos
        tmpedge(1) = -100; 
        tmpedge(end) = 70000;

        medge(1) = -100; 
        medge(end) = 70000;

        % all values in ith kymograph
        tmpky = ky{i};
        jj = 0;
        for j = tmpframes  % this is a vector from 1:nframes usually

            % if tmpframes has skipped frames, keep jj
            jj = jj + 1;
            % extract column vector of frame j from kymograph i
            tmpdata = tmpky(:,j);
            % extract column vector of frame j from memb kymograph
            tmpmb = tmpmemb(:,j);
            
            % some spatial bins in a given frame of a fixed kymo may be nans
            % these will have less than 301 points so ignore to ensure proper comparison
            % want fixed signal and spatial discretization
            % e.g. guv so small that # points < 301 spatial bins 
            if ~any(isnan(tmpdata)) 
                % compute entropy of frame j in kymograph i

                a1 = tmpdata;
                na1 = histcounts(a1,tmpedge)'; 
                Y = discretize(a1,tmpedge);    
                pixelCount = na1;
                grayLevels = [0:255]';

                a2 = tmpmb;
                na2 = histcounts(a2,medge)'; 
                Y2 = discretize(a2,medge);    
                pixelCount2 = na2;
                grayLevels2 = [0:255]';

                if any(isnan(Y)) == 1 % some data points are outside the min/max edges
                      disp('Some points are outside min/max edges of histogram')

                else
                    a1i = uint8(grayLevels(Y));
                    %a1i = im2uint8(grayLevels(Y));
                    a1h = na1;
                    a1h(a1h==0) = [];
                    a1p = a1h./numel(a1);
                    a1e = -sum(a1p.*log2(a1p));
                    a1e_mat = entropy(a1i);


                    a2i = uint8(grayLevels(Y2));
                    a2e_mat = entropy(a2i);


                    %tmpent = a1e/a2e_mat;
                    tmpent = a1e;
                    % store entropy
                    entky{i}(jj) = tmpent';
    
    %                 f1 = figure(); 
    %                 ax1 = imagesc(a1i); yline(0.5:length(a1)-0.5);
    %                 colormap(parula(256))
    %                 ax1.Parent.CLim = [0 255];
    %                 cc = colorbar('FontSize',fs);
    %                 cc.Label.String = ('8-bit pixel intensity');
    %                 ylabel('Spatial Points','FontSize',fs)
    %                 yticks(1:length(a1i))
    %                 yticklabels(1:length(a1i))
    %                 f1.Position(3) = 200;
    %                 
    %                 figure;
    %                 imhist(a1i,parula);grid on
    %                 ylim([0,5]);
    %                 
    %                 figure;
    %                 [pixelCount, grayLevels] = imhist(a1i,parula);
    %                 bar(grayLevels, pixelCount); % Plot it as a bar chart.
    %                 grid on;
    %                 ylim([0,5])
    %                 title('Histogram of pixel values', 'FontSize', fs);
    %                 xlabel('Re-scaled pixel values','FontSize',fs);
    %                 ylabel('Counts','FontSize',fs);
    %                 
    %                 figure;
    %                 bar(grayLevels, pixelCount/sum(pixelCount),BaseValue=3e-3); % Plot it as a bar chart.
    %                 set (gca,'yscale','log');
    %                 ylim([3e-3,1])
    %                 grid on;
    %                 title('Prob of pixel values', 'FontSize', fs);
    %                 xlabel('Re-scaled pixel values','FontSize',fs);
    %                 ylabel('Probability','FontSize',fs);
                    

                end
               
            else

                disp('There are some nans in the kymo')
                disp('guv #:'); k
                disp('kymo #:'); i
                disp('frame #:'); j
   
            end

        end % jth frame

    end % ith kymo

    dEntrp(k,:) = entky;

end % end of kth guv

end

function f1 = plotentropylocalshape(fignum,dEntrp,pstatalign,fileparams)
v2struct(fileparams);
v2struct(pstatalign);


tlsize = 14;
dcolors = [ 0      0.4470 0.7410; 
            0.8500 0.3250 0.0980;
            0.4660 0.6740 0.1880;
            0.4940 0.1840 0.5560;
            0.3010 0.7450 0.9330];

% [maxf,idmaxf] = max(dNframes); %GUV with most frames
% tickstep = df{idmaxf}.dataParams.tickstep;
% [maxt,~] = max(dTrapa); %GUV with latest time after rapa is on
% maxd = max(maxt-dTrapa);
% maxf1 = maxf + maxd;

tmpframes = [1:maxf1]';
tmpframes = tmpframes - maxt;

transp = 0.5; %transparency
nky = size(dEntrp,2); % # of kymographs for which entropy is computed

%    f1 = figure(fignum); clf(f1);
f1 = figure(fignum); 
hold on

tmpEallguv = nan(maxf1,n_analysis);
Kmed = cell(nky,1);


for i=1:nky % entropy for each kymo
    
    entky = dEntrp(:,i);%entropy of ith kymo,all GUVs

    for k = 1:n_analysis % each guv


        % align time series of entropy based on largest t_rapa
    
        tmpEky = nan(maxf1,1);
        ell  = maxt-dTrapa(k);


        tmpEky(ell+1:ell+dNframes(k)) = entky{k}; %kth GUV 
        tmpEallguv(:,k) = tmpEky;

        ax7 = plot(tmpframes,tmpEky,'-','LineWidth',0.5,'Color',[dcolors(i,:),transp]);
   

    end
    Kmed{i} = median(tmpEallguv,2,'omitnan');
    %ax8 = plot(Kmed,'LineWidth',2,'Color',[dcolors(i,:),1]);

end



ax8 = gobjects(nky,1);
for i = 1:nky
    ax8(i) = plot(tmpframes,Kmed{i},'LineWidth',2,'Color',[dcolors(i,:),1]);
end
xline(0,'LineWidth',2);

if do_absvals
legend(ax8(1:nky),["ActA";"Actin";"Abs.Val. Cumulative disp.";"Abs.Val. Curvature"],'Location','southeast',...
    'FontSize',tlsize,'Location','eastoutside')
else
legend(ax8(1:nky),["ActA";"Actin";"Cumulative disp.";"Curvature";"Membrane marker"],'Location','southeast',...
    'FontSize',tlsize,'Location','eastoutside')
end
%     legend('Velocity','Curvature','Cumul. displacement','Membrane marker','Actin','ActA',...
%         'Location','southeast')

% legend('ActA','Actin','Cumulative disp.','Curvature','Location','southeast')
%tmpframes = tmpframes - maxt;
xlim([tmpframes(1)-1,tmpframes(end)+1]); 
xticks([tmpframes(1):tickstep:tmpframes(end)])

ylim([0,8])
ax7.Parent.YLabel.String = 'Entropy';
ax7.Parent.YLabel.FontName = 'Arial';
%ax7.Parent.YLabel.FontSize = tlsize;
ax7.Parent.XLabel.String = 'Frames';
ax7.Parent.XLabel.FontName = 'Arial';
%ax7.Parent.XLabel.FontSize = tlsize;
% title(strcat(outputname,{' '},': entropy in boundary quantities'),'Interpreter','none',...
%     'FontSize',tlsize)

ax7.Parent.FontSize = tlsize;
title({outputname;'Entropy in boundary quantities'},'Interpreter','none',...
    'FontSize',tlsize+3)
 
hold off
f1.Position(3:4) = [800,500];

drawnow;
%     if newsavefigs
%         print(f1,strcat(newoutputfolder,guvname,'_',tstamp1,'_fig_globalentropy'),'-dpng')
%     end



end % end plotentropylocalshape


function [f1,Eallguv,Kmed,Sallguv,KmedS,tmpframes] = plotentropyglobalshape(fignum,dEChem,dShapes,pstatalign,fileparams)
v2struct(fileparams);
v2struct(pstatalign);
% v2struct(varshapes); % global shape features
% % recall: varshapes = v2struct(dEccent,dAspratio_ml,dAspratio_ind,dAreas);


nshapes = size(dShapes,2);


tlsize = 14;
dcolors = [0 0.4470 0.7410;  % ActA entropy
        0.8500 0.3250 0.0980; % actin entropy

        0.4660 0.6740 0.1880; % Eccent
        0.4940 0.1840 0.5560; % Aspratio_ml
        0.3010 0.7450 0.9330; % Aspration_ind
        0.6350 0.0780 0.1840]; % Areas


% [maxf,idmaxf] = max(dNframes); %GUV with most frames
% tickstep = df{idmaxf}.dataParams.tickstep;
% [maxt,~] = max(dTrapa); %GUV with latest time after rapa is on
% maxd = max(maxt-dTrapa);
% maxf1 = maxf + maxd;

tmpframes = [1:maxf1]';
tmpframes = tmpframes - maxt;

transp = 0.5; %transparency
nky = size(dEChem,2); % # of kymographs for which entropy is computed

%    f1 = figure(fignum); clf(f1);
f1 = figure(fignum); 
hold on


Eallguv = cell(nky,1);
Kmed = cell(nky,1);
ax8 = gobjects(nky+nshapes,1);

for i=1:nky % entropy for each kymo
    
    entky = dEChem(:,i);%entropy of ith kymo,all GUVs
    tmpEallguv = nan(maxf1,n_analysis);

    for k = 1:n_analysis % each guv


        % align time series of entropy based on largest t_rapa
    
        tmpEky = nan(maxf1,1);
        ell  = maxt-dTrapa(k);


        tmpEky(ell+1:ell+dNframes(k)) = entky{k}; %kth GUV 
        tmpEallguv(:,k) = tmpEky;

        ax7 = plot(tmpframes,tmpEky,'-','LineWidth',0.5,'Color',[dcolors(i,:),transp]);
   

    end
    Eallguv{i} = tmpEallguv;
    Kmed{i} = median(tmpEallguv,2,'omitnan');
    %ax8 = plot(Kmed,'LineWidth',2,'Color',[dcolors(i,:),1]);

end


for i = 1:nky
    ax8(i) = plot(tmpframes,Kmed{i},'LineWidth',3,'Color',[dcolors(i,:),1]);
end
xline(0,'LineWidth',2);

% legend(ax8(1:nky),["ActA";"Actin"],'FontSize',tlsize,'Location','eastoutside')

ylim([0,8])
ax7.Parent.YLabel.String = 'Entropy';
ax7.Parent.YLabel.FontName = 'Arial';

%% now plot shape parameters on right y-axis

yyaxis right

Sallguv = cell(nshapes,1);
KmedS = cell(nshapes,1);

for i = 1:nshapes

    shapeID= dShapes(:,i); % ith shape parameter of all GUVs
    tmpSallguv = nan(maxf1,n_analysis);


    for k = 1:n_analysis
        % align time series of entropy based on largest t_rapa
    
        tmpshape = nan(maxf1,1);
        ell  = maxt-dTrapa(k);


        tmpshape(ell+1:ell+dNframes(k)) = shapeID{k}; %kth GUV 
        tmpSallguv(:,k) = tmpshape;

        ax7 = plot(tmpframes,tmpshape,'-','LineWidth',0.5,'Color',[dcolors(i+nky,:),transp]);
   

    end

    Sallguv{i} = tmpSallguv;
    KmedS{i} = median(tmpSallguv,2,'omitnan');
    %ax8 = plot(Kmed,'LineWidth',2,'Color',[dcolors(i,:),1]);


end

% ax9 = gobjects(nshapes,1);
for i = 1:nshapes
    ax8(i+nky) = plot(tmpframes,KmedS{i},'-','LineWidth',3,'Color',[dcolors(i+nky,:),1]);
end


legend(ax8(1:end),["ActA";"Actin";"Eccent.";"AspRatio-Mat";"AspectRatio-Ind";"Norm. Area"],...
    'FontSize',tlsize,'Location','eastoutside')

% legend(ax8(nky+1:end),["Eccent.";"AspRatio-Mat";"AspectRatio-Ind";"Norm. Area"],...
%     'FontSize',tlsize,'Location','eastoutside')

% if do_absvals
% legend(ax8(1:nky),["ActA";"Actin";"Abs.Val. Cumulative disp.";"Abs.Val. Curvature"],'Location','southeast',...
%     'FontSize',tlsize,'Location','eastoutside')
% else
% legend(ax8(1:nky),["ActA";"Actin";"Cumulative disp.";"Curvature"],'Location','southeast',...
%     'FontSize',tlsize,'Location','eastoutside')
% end
%     legend('Velocity','Curvature','Cumul. displacement','Membrane marker','Actin','ActA',...
%         'Location','southeast')

% legend('ActA','Actin','Cumulative disp.','Curvature','Location','southeast')
%tmpframes = tmpframes - maxt;
xlim([tmpframes(1)-1,tmpframes(end)+1]); 
xticks([tmpframes(1):tickstep:tmpframes(end)])

ylim([0,0.85])
ax7.Parent.YLabel.String = 'Shape quantities';
ax7.Parent.YLabel.FontName = 'Arial';
%ax7.Parent.YLabel.FontSize = tlsize;
ax7.Parent.XLabel.String = 'Frames';
ax7.Parent.XLabel.FontName = 'Arial';
%ax7.Parent.XLabel.FontSize = tlsize;
% title(strcat(outputname,{' '},': entropy in boundary quantities'),'Interpreter','none',...
%     'FontSize',tlsize)

ax7.Parent.FontSize = tlsize;
title({outputname;'Entropy in boundary quantities'},'Interpreter','none',...
    'FontSize',tlsize+3)
 
hold off
f1.Position(3:4) = [800,500];

drawnow;
%     if newsavefigs
%         print(f1,strcat(newoutputfolder,guvname,'_',tstamp1,'_fig_globalentropy'),'-dpng')
%     end



end % end plotentropyglobalshape



function [R1out, C1out,RR1out,CC1out,scp] = spatialcorr(refidx,dP,dK,corrtype,pstatalign,fileparams)
% R1out correlation coefficients
% C1out # of spatial points meeting inclusion criteria
% plot spatial correlations for fixedkymos over time


v2struct(fileparams);
v2struct(pstatalign);
%double quotes for string array!!!
%ch = ["Membrane marker","Actin","ActA"];
ch = ["Membrane marker","ActA","Actin"];

% [maxf,idmaxf] = max(dNframes); %GUV with most frames
% %tickstep = dP{idmaxf}.tickstep; %df{idmaxf}.dataParams.tickstep;
% [maxt,~] = max(dTrapa); %GUV with latest time after rapa is on
% maxd = max(maxt-dTrapa);
% maxf1 = maxf + maxd;



xreflabel = ch(refidx);
othidx = setdiff([1:3],refidx,'sorted'); % channel numbers of other 

if do_absvals
    xothlabel = [ch(othidx(1)),ch(othidx(2)),"velocity","abs.val. cumulative disp.","abs.val. curvature"];
else
    xothlabel = [ch(othidx(1)),ch(othidx(2)),"velocity","cumulative disp.","curvature"];
end
% preallocate outputs
nxoth = length(xothlabel);
R1out = cell(n_analysis,nxoth);
C1out = cell(n_analysis,1); % of spatial points in each GUVkymo meeting criteria for correlation
tstring = strings(1,nxoth);
%nthresh = 100;
mthresh = nan(n_analysis,1);

allpadkymoref = cell(n_analysis,1);
allpadkymommb = cell(n_analysis,1);
allpadkymocomp = cell(n_analysis,nxoth);

RR1out = cell(nxoth,1);


for j = 1:nxoth% comparison variables

    % each comparison variable has a new figure 
    % will show correlation relative to reference variable
%     fignum = nfig0+j;
%     f11 = figure(fignum); hold on
     tstring(j) = strcat(xreflabel,'-',xothlabel(j)," correlation");



    for i = 1:n_analysis % GUV number
%         tmpkflr  = df{i}.dataKymoVar.kymofluor; % kymos for fluorescent channels
%         tmpkvel  = df{i}.dataKymoVar.kymovel; % kymo velocities
%         tmpkdc = df{i}.dataKymoVar.kymodcumul; % kymo cumulative disp
%         tmpkcurv  = df{i}.dataKymoVar.kymocurv; % kymo curvatures

        tmpkflr  = {dK(i).memb,dK(i).ActA,dK(i).actin}; % kymos for fluorescent channels
        tmpkvel  = dK(i).vel; % kymo velocities
        tmpkdc = dK(i).dc; % kymo cumulative disp
        tmpkcurv  = dK(i).curv; % kymo curvatures

        % normalization to ensure features are comparable across guvs


        if do_norm_corr == 1
            a1 = tmpkflr{1}; a2 = tmpkflr{2}; a3 = tmpkflr{3};
            tmpkflr{1} = a1./mean(a1,'omitnan');
            tmpkflr{2} = a2./mean(a2,'omitnan');
            tmpkflr{3} = a3./mean(a3,'omitnan');
        elseif do_norm_corr == 2
            a1 = tmpkflr{1}; a2 = tmpkflr{2}; a3 = tmpkflr{3};
            tmpkflr{1} = a1./vecnorm(a1);
            tmpkflr{2} = a2./vecnorm(a2);
            tmpkflr{3} = a3./vecnorm(a3);
        elseif do_norm_corr == 3
            a1 = tmpkflr{1}; a2 = tmpkflr{2}; a3 = tmpkflr{3};
            tmpkflr{1} = a1./max(a1(:));
            tmpkflr{2} = a2./max(a2(:));
            tmpkflr{3} = a3./max(a3(:));
        elseif do_norm_corr == 4
            a1 = tmpkflr{1}; a2 = tmpkflr{2}; a3 = tmpkflr{3};
            tmpkflr{1} = a1./max(a1,[],1);
            tmpkflr{2} = a2./max(a2,[],1);
            tmpkflr{3} = a3./max(a3,[],1);
        end

%         figure(101); imagesc(tmpkflr{1}); colorbar;
%         figure(102); imagesc(tmpkflr{2});colorbar;
%         figure(103); imagesc(tmpkflr{3});colorbar;
        % make cumulative displacement first frame all nans
        tmpkdc(:,1) = nan;
        
        % compute corr to absolute values - NOT MUCH DIFFERENT
        if do_absvals
            tmpkvel = abs(tmpkvel);
            tmpkdc = abs(tmpkdc);
            tmpkcurv = abs(tmpkcurv);
        end

        % collect all comparison kymos for a given GUV
        tmpcompall = [tmpkflr(othidx(1)),tmpkflr(othidx(2)),tmpkvel,tmpkdc,tmpkcurv];
        % reference kymo for given GUV
        kymoref0 = tmpkflr{refidx}; % ActA is reference for comparision with other signals
            

    %% --------- compute reference-comparison spatial correlation
      
        % jth comparison variable
        kymocomp0 = tmpcompall{j}; % jth comparison variable
    
        tmpnframes = dNframes(i);
        tmpR0 = nan(tmpnframes,1);
        tmpC0 = tmpR0;
        aa0 = ~isnan(kymoref0);
        totalreals = sum(aa0(:));

        mmb0 = tmpkflr{1};
        mmb_min = min(min(mmb0));
        mmb_max = max(max(mmb0));
        msort = sort(mmb0(:),'ascend','MissingPlacement','last');
        %mthresh = mmb_min+pthresh*(mmb_max-mmb_min);
        % threshold for membrane marker for each GUV
        mthresh(i) = msort(round(pthresh*totalreals));

         %ignore nans and ignore regions where membrane marker is high
        real0 = ~isnan(kymoref0);
        low0 = mmb0<=mthresh(i);
        sd0 = logical(real0.*low0); % selected data
        
%        % can be wrong when pthresh < 1
%        % sets high memb marker regions to 0
%         kymoref = kymoref0.*sd0;
%         kymocomp = kymocomp0.*sd0;
%         mmb = mmb0.*sd0;

    kymoref = kymoref0; 
    kymocomp = kymocomp0; 
    mmb = mmb0;
    % set all points not meeting inclusion criteria to nans
    kymoref(~sd0) = nan;
    kymocomp(~sd0) = nan;
    mmb(~sd0) = nan;

        for k = 1:tmpnframes    

%             %ignore nans and ignore regions where membrane marker is high
%             real1 = ~isnan(kymoref(:,k));
%             low1 = mmb(:,k)<=mthresh(i);
%             sd = logical(real1.*low1); % selected data
%            totalsd = sum(sd(:));

            totalsd = sum(sd0(:,k));
            sdk = sd0(:,k);

            if totalsd >= nthresh
                if do_partial_corr
                    %tmpR0(k) = partialcorr(kymoref(sd,k),kymocomp(sd,k),mmb(sd,k),'type',corrtype);
                    tmpR0(k) = partialcorr(kymoref(sdk,k),kymocomp(sdk,k),mmb(sdk,k),'type',corrtype);
                else
                    %tmpR0(k) = corr(kymoref(sd,k),kymocomp(sd,k),'type',corrtype);
                    tmpR0(k) = corr(kymoref(sdk,k),kymocomp(sdk,k),'type',corrtype);
                end
            else
                tmpR0(k) = nan;
            end
            tmpC0(k) = totalsd;
        end
            



        tmpR1 = nan(maxf1,1);
        tmpC1 = nan(maxf1,1);
        ell  = maxt-dTrapa(i);
        tmpR1(ell+1:ell+dNframes(i)) = tmpR0;
        tmpC1(ell+1:ell+dNframes(i)) = tmpC0;


        R1out{i,j} = tmpR1;
        C1out{i} = tmpC1; % counts, # meeting inc criteria, same for all correlation variables

        
        % use the filtered kymo  
        padkcomp = nan(size(kymocomp,1),maxf1);        
        padkcomp(:,ell+1:ell+dNframes(i)) = kymocomp;        
        allpadkymocomp{i,j} = padkcomp;

        if j == 1 % only do this once
        padkref = nan(size(kymoref,1),maxf1);
        padkref(:,ell+1:ell+dNframes(i)) = kymoref;
        allpadkymoref{i} = padkref;

        padkmmb = nan(size(mmb,1),maxf1);
        padkmmb(:,ell+1:ell+dNframes(i)) = mmb;
        allpadkymommb{i} = padkmmb;


        end

            
    end % end of ith GUV analysis



   % these kymos have same number of columns
   % stack the rows
    allkref = cat(1,allpadkymoref{:});
    allkcompj = cat(1,allpadkymocomp{:,j});
    allkmmb = cat(1,allpadkymommb{:});

    allsd = ~isnan(allkref);

    padframes = size(allkref,2);

    RR1tmp = nan(padframes,1);
    CC1tmp = RR1tmp;

    for m = 1:padframes

        sdm = allsd(:,m);
        totalsdm = sum(sdm);
        if totalsdm >= nthresh
            if do_partial_corr
                [RR1tmp(m),pval] = partialcorr(allkref(sdm,m),allkcompj(sdm,m),allkmmb(sdm,m),'type',corrtype);
            else
                [RR1tmp(m),pval] = corr(allkref(sdm,m),allkcompj(sdm,m),'type',corrtype);
            end
        else
            RR1tmp(m) = nan;
        end
        CC1tmp(m) = totalsdm;

    end

%     for m = 50
% 
%         sdm = allsd(:,m);
%         totalsdm = sum(sdm);
%         if totalsdm >= nthresh
%             if do_partial_corr
%                 [RR1tmp(m),pval] = partialcorr(allkref(sdm,m),allkcompj(sdm,m),allkmmb(sdm,m),'type',corrtype);
%             else
%                 [RR1tmp(m),pval] = corr(allkref(sdm,m),allkcompj(sdm,m),'type',corrtype);
%             end
%         else
%             RR1tmp(m) = nan;
%         end
%         CC1tmp(m) = totalsdm;
% 
%     end
%     


    RR1out{j} = RR1tmp;
    CC1out = CC1tmp;
    


end % end of jth comparison variable   

scp = v2struct(maxt,maxf1,tstring,outputname,corrtype,xothlabel,xreflabel,pthresh,mthresh,nthresh);

end


function [figsout,R1med] = plotspatialcorr(nfig0,R1,C1,RR1,CC1,scparams,pstatalign,fileparams)

% recall:   scp = v2struct(maxt,maxf1,tstring,outputname,corrtype,xothlabel,xreflabel,mthresh,pthresh);
v2struct(scparams); 
v2struct(fileparams);
v2struct(pstatalign);

tlsize = 14;
% # of correlation plots (5)
nxoth = length(xothlabel);
figsout = gobjects(nxoth+2,1);
R1med = cell(nxoth,1);
tmpframes = [1:maxf1]';
tmpframes = tmpframes - maxt;

for j = 1:nxoth
    fignum = nfig0+j;

    
    for i = 1:n_analysis

        f11 = figure(fignum); hold on
        
        ax7 = gca;
        p1 = plot(tmpframes,R1{i,j},'LineWidth',1.5);
        % set transparency and color
        p1.Color(1:3) = [0 0.4470 0.7410]; 
        p1.Color(4) = 0.5;

        if n_analysis == 1
            text(0.7*maxf1,0.8,{['memb_cutoff: ',num2str(mthresh(i))],['prop_cutoff: ',num2str(pthresh)]},...
            'Interpreter','none','FontSize',12);
        end

    end

    if n_analysis > 1
        text(0.7*maxf1,0.8,['overall_prop: ',num2str(pthresh)],...
        'Interpreter','none','FontSize',tlsize-2);
    end


   
%     R1cat = cat(2,R1{:,j});
%     R1med = median(R1cat','omitnan');
%     plot(R1med,'k' ,'LineWidth',2.5)



%     hold off

%     title(ax7,strcat(outputname,{' '},tstring{j},' (',corrtype,')'),'Interpreter','none',...
%     'FontSize',tlsize);

    title(ax7,{outputname; strcat(tstring{j},' (',corrtype,')')},'Interpreter','none',...
    'FontSize',tlsize);

       
    R1cat = cat(2,R1{:,j});
    R1med{j} = median(R1cat','omitnan')';
    plot(tmpframes,R1med{j},'k' ,'LineWidth',2.5)

% bulk frame to frame
    plot(tmpframes,RR1{j},'r','LineWidth',2.5)


    
    xline(0,'LineWidth',2); %maxt
    yline(0);

    ax7.YLim=[-1,1];
    ax7.XLim=[0, maxf1];

    ax7.YLabel.String = strcat(corrtype,{' '},'correlation coef.');
    ax7.YLabel.FontName = 'Arial';
    %ax7.YLabel.FontSize = tlsize;
    ax7.XLabel.String = 'Frames';
    ax7.XLabel.FontName = 'Arial';
    %ax7.XLabel.FontSize = tlsize;
    ax7.FontSize = tlsize;

xlim([tmpframes(1)-1,tmpframes(end)+1]); 
xticks([tmpframes(1):tickstep:tmpframes(end)])

    legend(ax7,[dguvnames(:);"median";"bulk";'';''],'Location','eastoutside');

    % using [ ] allows separate lines for comma-separated text

%     title(ax7.Parent,strcat(outputname,{' '},tstring,' (',corrtype,')'),...
%         'FontSize',tlsize,'Interpreter','none');
    hold off
    
    f11.Position(3:4) = [800,500];
    drawnow;
    figsout(j) = f11;

end


%% plot # of points used for analysis
fignum = nfig0+nxoth+1;
f12 = figure(fignum); hold on

for i = 1:n_analysis
    plot(C1{i})
    ylabel('# of points for corr. (below memb cutoff)','FontSize',tlsize)
end
yline(nthresh,'LineWidth',2);
xline(maxt,'LineWidth',2);
% C1cat = cat(2,C1{:});
% C1max = max(max(C1cat));
yl = ylim;

legend(dguvnames(:),'Location','eastoutside');
% title(strcat(outputname,{' '},'# of spatial points used for each GUV at each frame',' (',corrtype,')'),'Interpreter','none',...
%     'FontSize',tlsize);
title({outputname;strcat('# of spatial points used for each GUV at each frame',' (',corrtype,')')}...
    ,'Interpreter','none','FontSize',tlsize);


text(0.7*maxf1,0.8*yl(2),['overall_prop: ',num2str(pthresh)],...
    'Interpreter','none','FontSize',tlsize-2);

hold off
f12.Position(3:4) = [800,500];
figsout(end-1) = f12;


%% plot # of points used for bulk analysis
fignum = nfig0+nxoth+2;
f13 = figure(fignum); hold on

for i = 1:n_analysis
    plot(CC1)
    ylabel('# of points for corr. (below memb cutoff)','FontSize',tlsize)
end
yline(nthresh,'LineWidth',2);
xline(maxt,'LineWidth',2);
% C1cat = cat(2,C1{:});
% C1max = max(max(C1cat));
yl = ylim;

%legend(dguvnames(:),'Location','eastoutside');
% title(strcat(outputname,{' '},'# of spatial points used for each GUV at each frame',' (',corrtype,')'),'Interpreter','none',...
%     'FontSize',tlsize);
title({outputname;strcat('# of spatial points from all GUVs at each frame',' (',corrtype,')')}...
    ,'Interpreter','none','FontSize',tlsize);


text(0.7*maxf1,0.8*yl(2),['overall_prop: ',num2str(pthresh)],...
    'Interpreter','none','FontSize',tlsize-2);

hold off
f13.Position(3:4) = [800,500];
figsout(end) = f13;

end


function [R1out,C1out,pval,refAgg,compAgg,scp] = signalcorr(refidx,dP,dK,corrtypeSig,pstatalign,fileparams)


v2struct(fileparams);
v2struct(pstatalign);
%double quotes for string array!!!
%ch = ["Membrane marker","Actin","ActA"];
ch = ["Membrane marker","ActA","Actin"];

% [maxf,idmaxf] = max(dNframes); %GUV with most frames
% tickstep = dP{idmaxf}.tickstep; %df{idmaxf}.dataParams.tickstep;
% [maxt,~] = max(dTrapa); %GUV with latest time after rapa is on
% maxd = max(maxt-dTrapa);
% maxf1 = maxf + maxd;


xreflabel = ch(refidx);
refstr = strcat(xreflabel," (AU)");

othidx = setdiff([1:3],refidx,'sorted'); % channel numbers of other 
if do_absvals
    xothlabel = [ch(othidx(1)),ch(othidx(2)),"velocity","abs.val. cumulative disp.","abs.val. curvature"];
    
    compstr = [strcat(ch(othidx(1))," (AU)");
               strcat(ch(othidx(2))," (AU)");
               "velocity (px/frame)";
               "abs.val. cumulative disp. (px)";
               "abs.val. curvature (1/px)"];
else
    xothlabel = [ch(othidx(1)),ch(othidx(2)),"velocity","cumulative disp.","curvature"];
    
    compstr = [strcat(ch(othidx(1))," (AU)");
               strcat(ch(othidx(2))," (AU)");
               "velocity (px/frame)";
               "cumulative disp. (px)";
               "curvature (1/px)"];
end



% preallocate outputs
nxoth = length(xothlabel);
refAgg = cell(n_analysis,2); %aggregated kymoref in 2 matrices - before /after Rapa
compAgg = cell(n_analysis,2,nxoth); % aggregated other kymos in 2 matrices - before /after Rapa
mmbAgg = refAgg;

R1out = nan(2,nxoth); % full correlations before/after rapa between ref and comparison variables
C1out = nan(2,nxoth); % # of data points for before/after rapa full corr
C1ref = nan(n_analysis,2); % # of data points in ref variable meeting inc criteria before/after 
pval = nan(2,nxoth);

tstring = strings(1,nxoth);

%nthresh = 100; % min number of points per GUV
mthresh = nan(n_analysis,1);

for j = 1:nxoth
    % each comparison variable has a new figure 
    tstring(j) = strcat(xreflabel,'-',xothlabel(j)," correlation");
end

for i = 1:n_analysis % GUV number

    tmpkflr  = {dK(i).memb,dK(i).ActA,dK(i).actin}; % kymos for fluorescent channels
    tmpkvel  = dK(i).vel; % kymo velocities
    tmpkdc = dK(i).dc; % kymo cumulative disp
    tmpkcurv  = dK(i).curv; % kymo curvatures

    if do_norm_corr == 1
        a1 = tmpkflr{1}; a2 = tmpkflr{2}; a3 = tmpkflr{3};
        tmpkflr{1} = a1./mean(a1,'omitnan');
        tmpkflr{2} = a2./mean(a2,'omitnan');
        tmpkflr{3} = a3./mean(a3,'omitnan');
    elseif do_norm_corr == 2
        a1 = tmpkflr{1}; a2 = tmpkflr{2}; a3 = tmpkflr{3};
        tmpkflr{1} = a1./vecnorm(a1);
        tmpkflr{2} = a2./vecnorm(a2);
        tmpkflr{3} = a3./vecnorm(a3);
    elseif do_norm_corr == 3
        a1 = tmpkflr{1}; a2 = tmpkflr{2}; a3 = tmpkflr{3};
        tmpkflr{1} = a1./max(a1(:));
        tmpkflr{2} = a2./max(a2(:));
        tmpkflr{3} = a3./max(a3(:));
    elseif do_norm_corr == 4
        a1 = tmpkflr{1}; a2 = tmpkflr{2}; a3 = tmpkflr{3};
        tmpkflr{1} = a1./max(a1,[],1);
        tmpkflr{2} = a2./max(a2,[],1);
        tmpkflr{3} = a3./max(a3,[],1);
    end


   
    tmptrapa = dTrapa(i);
    tmpnpts = size(tmpkcurv,1); % # of points

    % make cumulative displacement first frame all nans
    tmpkdc(:,1) = nan;

    
    % compute corr to absolute values - NOT MUCH DIFFERENT
    if do_absvals
        tmpkvel = abs(tmpkvel);
        tmpkdc = abs(tmpkdc);
        tmpkcurv = abs(tmpkcurv);
    end
    
    % collect all comparison kymos for a given GUV
    tmpcompall = [tmpkflr(othidx(1)),tmpkflr(othidx(2)),tmpkvel,tmpkdc,tmpkcurv];
    % reference kymo for given GUV
    kymoref = tmpkflr{refidx}; % reference for comparision with other signals

    
    % setup inclusion criteria for REFERENCE CHANNEL only
    aa0 = ~isnan(kymoref); % ignore nans in ref channel
    totalreals = sum(aa0(:));

    mmb = tmpkflr{1};
    mmb_min = min(min(mmb));
    mmb_max = max(max(mmb));
    msort = sort(mmb(:),'ascend','MissingPlacement','last');
    %mthresh = mmb_min+pthresh*(mmb_max-mmb_min);
    % threshold for membrane marker for each GUV
    mthresh(i) = msort(round(pthresh*totalreals));
%     if pthresh == 1
%         mthresh(i) = max(mmb(:));
% 
%     end
%    

    %ignore nans and ignore regions where membrane marker is high
    aa1 = mmb<=mthresh(i); 
    sdref = logical(aa0.*aa1);


    
    %% set time windows for "pre" and "post" scatter plots
    
    %stimes1 = 3;
    stimes1 = 1 : tmptrapa;
    %stimes1 = tmptrapa - 5: tmptrapa;

    %stimes2 = 15;
    %stimes2 = tmptrapa+1:dNframes(i);
    stimes2 = tmptrapa + 10: tmptrapa + 44;




    %% count number of points meeting inclusion criteria before/after input
    totalsd1 = sum(sdref(:,stimes1),'all');



    totalsd2 = sum(sdref(:,stimes2),'all');

    % for bulk, allow more data
    nthresh2 = nthresh/2;
    if totalsd1 >= nthresh2 && totalsd2 >= nthresh2
        krefadj = kymoref;
        kmmbadj = mmb;
        krefadj(~sdref) = nan; % set all points not meeting inclusion criteria to nans
        kmmbadj(~sdref) = nan; 
    else
        krefadj = nan(size(kymoref));
        kmmbadj = nan(size(mmb));
        disp(strcat('GUV number =',num2str(Gidx_analysis(i))))
        disp('not enough points meeting nthresh')
    end


    % aligned pre and post rapa matrices
    tA0 = nan(tmpnpts,maxt);
    tA1 = nan(tmpnpts,maxf1-maxt);
    ell = maxt - dTrapa(i);


    tA0(:,ell+1:ell+length(stimes1)) = krefadj(:,stimes1);
    %tA1(:,1:dNframes(i)-(tmptrapa+lag1)) = krefadj(:,(tmptrapa+lag1)+1:end);
    tA1(:,1:length(stimes2)) = krefadj(:,stimes2);
    refAgg{i,1} = tA0;
    refAgg{i,2} = tA1;

    mA0 = nan(tmpnpts,maxt);
    mA1 = nan(tmpnpts,maxf1-maxt);

    mA0(:,ell+1:ell+length(stimes1)) = kmmbadj(:,stimes1);
    %mA1(:,1:dNframes(i)-(tmptrapa+lag1)) = kmmbadj(:,(tmptrapa+lag1)+1:end);
    mA1(:,1:length(stimes2)) = kmmbadj(:,stimes2);
    mmbAgg{i,1} = mA0;
    mmbAgg{i,2} = mA1;

    for k = 1:nxoth
        tC0 = nan(tmpnpts,maxt);
        tC1 = nan(tmpnpts,maxf1-maxt);
        % ignore nans and high memb marker regions
        if totalsd1 >= nthresh2 && totalsd2 >= nthresh2
            kcmpadj = tmpcompall{k};
            kcmpadj(~sdref) = nan; % set all points not meeting inclusion criteria to nans
        else
            kcmpadj = nan(size(kymoref));
        end      


        tC0(:,ell+1:ell+length(stimes1)) = kcmpadj(:,stimes1);
        %tC1(:,1:dNframes(i)-(tmptrapa+lag1)) = kcmpadj(:,(tmptrapa+lag1)+1:end);
        tC1(:,1:length(stimes2)) = kcmpadj(:,stimes2);

        compAgg{i,1,k} = tC0;
        compAgg{i,2,k} = tC1;            

    end



    % # of points meeting inclusion crit. in ref variable, for ith guv
    C1ref(i,:) = [totalsd1,totalsd2];
end




% now compute correlations between ref variable and comparison variables over all frames
% meeting inclusion criteria
for k = 1:nxoth
    % these are already filtered for inclusion criteria
    kA1 = refAgg(:,1); % ref kymo for all GUVs, pre-rapa
    kA1 = cat(1,kA1{:}); % stack ref kymos vertically
    kC1 = compAgg(:,1,k); % kth comparison kymos for all GUVs, pre-rapa,  
    kC1 = cat(1,kC1{:}); % stack kth comp kymos vertically

    kM1 = mmbAgg(:,1); % memb kymo for all GUVs, pre-rapa
    kM1 = cat(1,kM1{:}); % stack mmb kymos vertically

    kA2 = refAgg(:,2); % ref kymo for all GUVs, post-rapa
    kA2 = cat(1,kA2{:}); % stack ref kymos vertically
    kC2 = compAgg(:,2,k); % kth comparison kymos for all GUVs, post-rapa,  
    kC2 = cat(1,kC2{:}); % stack kth comp kymos vertically

    kM2 = mmbAgg(:,2); % memb kymo for all GUVs, pre-rapa
    kM2 = cat(1,kM2{:}); % stack mmb kymos vertically

    % find points with real values in both ref and comp variable
    % there's a mismatch only for velocity and cumuldips bc of all nans in first frame
    rA1 = ~isnan(kA1);
    rC1 = ~isnan(kC1);
    r1idx = logical(rA1.*rC1); % find points with real values in both ref and comp variable

    rA2 = ~isnan(kA2);
    rC2 = ~isnan(kC2);
    r2idx = logical(rA2.*rC2);
    
    s1idx = sum(r1idx(:));
    s2idx = sum(r2idx(:));

    %vector of reals
    kA1f = kA1(r1idx);
    kC1f = kC1(r1idx);
    kM1f = kM1(r1idx);

    kA2f = kA2(r2idx);
    kC2f = kC2(r2idx);
    kM2f = kM2(r2idx);

    
    if do_partial_corr
        [R1out(1,k),pval(1,k)] = partialcorr(kA1f,kC1f,kM1f,'type',corrtypeSig);
   
        [R1out(2,k),pval(2,k)] = partialcorr(kA2f,kC2f,kM2f,'type',corrtypeSig);

    else

        if ~isempty(kA1f)

        [R1out(1,k),pval(1,k)] = corr(kA1f,kC1f,'type',corrtypeSig);
   
        [R1out(2,k),pval(2,k)] = corr(kA2f,kC2f,'type',corrtypeSig);

        else


        end

    end

    C1out(:,k) = [s1idx;s2idx];

end



scp = v2struct(maxt,maxf1,tstring,xothlabel,xreflabel,mthresh,nthresh,...
    refidx,othidx,corrtypeSig,refstr,compstr);



% 
% 
% 
% guvname = dParams.guvname;
% tstamp1 = dParams.tstamp1;
% nbins = dParams.nbins;
% t_rapa = dParams.t_rapa;
% 
% kymovel  = dKtmp.vel; % kymo velocities
% kymocurv  = dKtmp.curv; % kymo curvatures
% kymodcumul = dKtmp.dc; % kymo cumulative disp
% kymofluor  = {dKtmp.memb,dKtmp.ActA,dKtmp.actin}; % kymos for fluorescent channels
% 
% 
% % if guvname(1) == "L"
% %     local = patches{1};
% %     distant = patches{2};
% % else
% %     local = 1:nbins;
% %     distant = 1:nbins;
% % 
% % end
% 
% region = patch;
% %allpts = 1:nbins;
% 
% kymoref = kymofluor{3}; %1 memb, 2 ActA, 3 actin
% xplabel = 'Actin';
% xstring = [xplabel,{' '},'intensity (AU)'];
% 
% %% plot actin-memb
% 
% 
% fignum = 40;
% ystring = 'Membrane intensity (AU)';
% tstring = 'actin-membrane correlation';
% scp = v2struct(region,t_rapa,ystring,tstring,guvname);
% kymocomp = kymofluor{1};
% figsout(1) = plotsignalcorr(fignum,kymoref,kymocomp,scp);
% 
% %% plot actin-ActA
% 
% 
% fignum = 41;
% ystring = 'ActA intensity (AU)';
% tstring = 'actin-ActA correlation';
% scp = v2struct(region,t_rapa,ystring,tstring,guvname);
% kymocomp = kymofluor{2};
% figsout(2) = plotsignalcorr(fignum,kymoref,kymocomp,scp);
% 
% %% ----- plot actin-velocity correlation 
% 
% fignum = 42;
% ystring = 'Membrane velocity (pixels/frame)';
% tstring = [xplabel,'-velocity correlation'];
% scp = v2struct(region,t_rapa,xstring,ystring,tstring,guvname);
% kymocomp = kymovel;
% figsout(3) = plotsignalcorr(fignum,kymoref,kymocomp,scp);
% 
% %
% %% plot actin- cumulative displacement
% 
% fignum = 43;
% ystring = 'Cumulative displacement (pixels)';
% tstring = 'actin-cumulative displacement correlation';
% scp = v2struct(region,t_rapa,ystring,tstring,guvname);
% kymocomp = kymodcumul;
% figsout(4) = plotsignalcorr(fignum,kymoref,kymocomp,scp);
% 
% %% --------- plot actin-curvature full correlation
% 
% 
% fignum = 44;
% ystring = 'Membrane curvature (1/pixels)';
% tstring = 'actin-curvature correlation';
% scp = v2struct(region,t_rapa,ystring,tstring,guvname);
% kymocomp = kymocurv;
% figsout(5) = plotsignalcorr(fignum,kymoref,kymocomp,scp);
% 
%  

end



function figsout = plotsignalcorr(nfig0,R1,C1,refAgg,compAgg,scparams,pstatalign,fileparams)

% recall 
% scp = v2struct(maxt,maxf1,tstring,xothlabel,xreflabel,mthresh,nthresh,...
%    refidx,othidx,corrtypeSig,refstr,compstr);

v2struct(scparams); 
v2struct(fileparams);
v2struct(pstatalign);

tlsize = 14;
% # of correlation plots (5)
nxoth = length(xothlabel);
figsout = gobjects(nxoth+1,1);
plims = struct2cell(pstatalign.pranges);
cmpidx = [othidx,4:6];

% recall
% refAgg = cell(n_analysis,2); %aggregated kymoref in 2 matrices - before /after Rapa
% compAgg = cell(n_analysis,2,nxoth); % aggregated other kymos in 2 matrices - before /after Rapa
% R1out = nan(2,nxoth); % full correlations before/after rapa between ref and comparison variables
% C1out = nan(2,nxoth); % # of data points for before/after rapa full corr

for k = 1:nxoth
    fignum = nfig0+k;
    
    kA1i = refAgg(:,1); % ref kymo for all GUVs, pre-rapa
    kA1 = cat(1,kA1i{:}); % stack ref kymos vertically
    kC1i = compAgg(:,1,k); % kth comparison kymos for all GUVs, pre-rapa,  
    kC1 = cat(1,kC1i{:}); % stack kth comp kymos vertically

    kA2i = refAgg(:,2); % ref kymo for all GUVs, post-rapa
    kA2 = cat(1,kA2i{:}); % stack ref kymos vertically
    kC2i = compAgg(:,2,k); % kth comparison kymos for all GUVs, post-rapa,  
    kC2 = cat(1,kC2i{:}); % stack kth comp kymos vertically


    f14 = figure(fignum); 
    t14 = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
    
    f14a = nexttile(1,[2 2]);   

    ax7 = binscatter(kA1(:),kC1(:),200);
    colormap(gca,'winter');
%     legend('test')
 
    % axis limits for ref and comp variable

    xmin = min([kA1(:);kA2(:)]); xmax = max([kA1(:);kA2(:)]);
    ymin = min([kC1(:);kC2(:)]); ymax = max([kC1(:);kC2(:)]);
    

    ax7.Parent.YLim=[ymin, ymin + 1.2*(ymax-ymin)];
    ax7.Parent.XLim=[xmin, xmin + 1.2*(xmax-xmin)];
   

    ax7.Parent.YLabel.String = compstr(k);
    ax7.Parent.YLabel.FontName = 'Arial';
    ax7.Parent.YLabel.FontSize = tlsize;
    ax7.Parent.XLabel.String = refstr;
    ax7.Parent.XLabel.FontName = 'Arial';
    ax7.Parent.XLabel.FontSize = tlsize;


    xl=xmin+1*(xmax-xmin);
    yl=ymin+1.1*(ymax-ymin);
    text(xl,yl,sprintf(strcat(corrtypeSig,...
        ' Corr.: %5.2f'),R1(1,k)),'FontSize',tlsize-2,...
    'HorizontalAlignment','right')

    xlb=xmin+1*(xmax-xmin);
    ylb=ymin+1.15*(ymax-ymin);
    text(xlb,ylb,['overall_prop: ',num2str(pthresh)],...
        'Interpreter','none','FontSize',tlsize-2,'HorizontalAlignment','right')

    title(f14a,strcat(tstring(k),': Pre-rapamycin'),'FontSize',tlsize)




    f14b = nexttile(3,[2 2]);
   
    ax7 = binscatter(kA2(:),kC2(:),200);
    colormap(gca,'spring');
 

    ax7.Parent.YLim=[ymin, ymin + 1.2*(ymax-ymin)];
    ax7.Parent.XLim=[xmin, xmin + 1.2*(xmax-xmin)];

    ax7.Parent.YLabel.String = compstr(k);
    ax7.Parent.YLabel.FontName = 'Arial';
    ax7.Parent.YLabel.FontSize = tlsize;
    ax7.Parent.XLabel.String = refstr;
    ax7.Parent.XLabel.FontName = 'Arial';
    ax7.Parent.XLabel.FontSize = tlsize;


    xl=xmin+1*(xmax-xmin);
    yl=ymin+1.1*(ymax-ymin);
    text(xl,yl,sprintf(strcat(corrtypeSig,...
        ' Corr.: %5.2f'),R1(2,k)),'FontSize',tlsize-2,...
    'HorizontalAlignment','right')

    xlb=xmin+1*(xmax-xmin);
    ylb=ymin+1.15*(ymax-ymin);
    text(xlb,ylb,['overall_prop: ',num2str(pthresh)],...
        'Interpreter','none','FontSize',tlsize-2,'HorizontalAlignment','right')

    title(f14b,strcat(tstring(k),': Post-rapamycin'),'FontSize',tlsize)

    hold off
    
    f14.Position(3:4) = [800,400];
    drawnow;
    figsout(k) = f14;

end



%% plot # of points used for analysis
fignum = nfig0+nxoth+1;
f15 = figure(fignum); hold on

bar(C1)
ylabel('# of points for corr. (below memb cutoff)','FontSize',tlsize)

yline(nthresh,'LineWidth',2);
% C1cat = cat(2,C1{:});
% C1max = max(max(C1cat));
yl = ylim;

legend(tstring(:),'Location','west','FontSize',tlsize);
title({string(outputname);strcat("# of spatial points used for each correlation",...
    " (",corrtypeSig,")")},'Interpreter','none',...
    'FontSize',tlsize);
text(1.5,0.8*yl(2),['overall_prop: ',num2str(pthresh)],...
    'Interpreter','none','FontSize',tlsize-2,'HorizontalAlignment','right');

hold off
%f15.Position(3:4) = [800,500];
figsout(end) = f15;

end



function plotspatiotempcorr(allData,patches,corrtype,params)

v2struct(params);

cent   = allData.allCentroids;
areas  = allData.allAreas;
cols   = allData.cols;
rows   = allData.rows;
nframes = allData.nframes;
listFrames = allData.listFrames;
fnames  = allData.fnames; % name of channels
kymovel  = allData.kymovel; % kymo velocities
kymocurv  = allData.kymocurv; % kymo curvatures
kymodcumul = allData.kymodcumul; % kymo cumulative disp
kymofluor  = allData.kymofluor; % kymos for fluorescent channels
membch = allData.membch;
fluorch = allData.fluorch;
guvname = allData.guvname;
savgolay_smooth = allData.savgolay_smooth;
computecurv = allData.computecurv;
im0    = allData.im0;
imL    = allData.imL;
allImages  = allData.allImages;
allIrgb = allData.allIrgb;
savefigs = allData.savefigs;
tstamp1 = allData.tstamp1;
outputfolder = allData.outputfolder;
alignbdy = allData.alignbdy; %allBoundaries;
alignvel = allData.alignvel; %allVelocities;
aligncurv = allData.aligncurv; 
alignfluor = allData.alignfluor; %allIntensities;
nbins = allData.nbins;
deltaf = allData.deltaf;
thetacenter = allData.thetacenter;
t_rapa = allData.t_rapa;
vflag = any(~isnan(kymovel(:))); % check if velocity is available
kymofixedtheta = allData.kymofixedtheta;
kymofixedpts = allData.kymofixedpts;

allpatches = cat(1,patches{1},patches{2});

    
%% --------- plot actin-curvature full correlation

tlist = cell(1,3);
tlist{1} = 1:t_rapa;
tlist{2} = t_rapa+1:nframes;
%tlist{3} = 1:nframes;
for i = 1:2
    f20 = figure(20+i);


    k3 = kymocurv(allpatches,tlist{i})'; % curvature from -pi to pi columns
    k2 = kymofluor{2}(allpatches,tlist{i})'; % actin channel, compare over all positions, same frames as vel
    %s = k2>0; % useful to avoid nans, compare only when normalized actin signal is > 0 
%     k3 = k3(s); 
%     k2 = k2(s);
    k2k3data = [k2, k3];
    n1 = length(patches{1});
    n2 = length(patches{2});
    b = [0,n1,n1+n2,2*n1+n2];
    bl = b+0.5;
    c = floor(nbins/4); 
    markticks = b+c; %sort([b,b+c],'ascend'); % b + c
    % matrix of correlation coefficient between pairs of columns of input matrix
    R = corr(k2k3data,'type',corrtype);

    ax1=imagesc(R); axis equal; hold on
    xline([bl(2),bl(4)],'LineWidth',1); xline(bl(3),'LineWidth',2.5)
    yline([bl(2),bl(4)],'LineWidth',1); yline(bl(3),'LineWidth',2.5)
    ax1.Parent.CLim = [-1,1]; colorbar;
    set(gca,'YDir','normal')
    ax1.Parent.YLim=[0.5, size(k2k3data,2)+0.5];
    ax1.Parent.XLim=[0.5, size(k2k3data,2)+0.5];
    xticks(markticks); yticks(markticks);
    if kymofixedtheta 
        xticklabels({'Actin: \theta=0','Actin: \theta=\pi','Curv: \theta=0'...
            ,'Curv: \theta=\pi'})
        yticklabels({'Actin: \theta=0','Actin: \theta=\pi','Curv: \theta=0'...
            ,'Curv: \theta=\pi'})
        if guvname(1) == 'L'
            xlabel({'','\theta = angular position relative to micropipette'})
        else
            thetac = num2str(wrapToPi(deg2rad(thetacenter)/pi));
            xlabel({'',['\theta = angular position relative to reference angle (',thetac,'\pi)']});
        end
    elseif kymofixedpts && guvname(1) == 'L'
        xticklabels({'Actin: Bins_{near}','Actin: Bins_{far}','Curv: Bins_{near}','Curv: Bins_{far}'})
        yticklabels({'Actin: Bins_{near}','Actin: Bins_{far}','Curv: Bins_{near}','Curv: Bins_{far}'})
    elseif kymofixedpts && guvname(1) == 'G'
        xticklabels({'Actin: Bins_{center}','Actin: Bins_{last}','Curv: Bins_{center}','Curv: Bins_{last}'})
        yticklabels({'Actin: Bins_{center}','Actin: Bins_{last}','Curv: Bins_{center}','Curv: Bins_{last}'})
    end

    %ax1.Parent.XLabel.FontName = 'Arial';

    ax1.Parent.XAxis.FontSize = 12;
    ax1.Parent.YAxis.FontSize = 12;
    ax1.Parent.XLabel.FontSize = 11; % must come after XAxis, bc its a subset of those

    if i == 1 
        titlelist = [guvname,': Pre-Rapamacyin'];
    elseif i ==2
        titlelist = [guvname,': Post-Rapamacyin'];
%     elseif i == 3
%         titlelist = [guvname,': All frames'];
    end

%     title({guvname,'Actin-curvature correlation coefficients'},...
%         {'computed for time-series of different angular positions'})
    % title ( main title with multiple lines, subtitle )
    title({titlelist,['Actin-curvature ',corrtype,' corr.']},...
        {'computed for time-series of different angular positions'})

    hold off
    
    f20.Position(3:4) = [700,600];
    %shading  flat; 
    drawnow;




    if newsavefigs && i == 1
        print(f20,strcat(newoutputfolder,guvname,'_',tstamp1,'_fig_actincurv_preRapcorr_',corrtype(1:5)),'-dpng')

    elseif newsavefigs && i == 2
        print(f20,strcat(newoutputfolder,guvname,'_',tstamp1,'_fig_actincurv_postRapcorr_',corrtype(1:5)),'-dpng')

%     elseif newsavefigs && i == 3
% 
%         print(f20,strcat(newoutputfolder,guvname,'_',tstamp1,'_fig_actincurv_fullcorr'),'-dpng')

    end

end


%% --------- plot actin-membrane marker full correlation

tlist = cell(1,3);
tlist{1} = 1:t_rapa;
tlist{2} = t_rapa+1:nframes;
%tlist{3} = 1:nframes;
for i = 1:2
    f30 = figure(30+i);


    k3 = kymofluor{1}(allpatches,tlist{i})'; % curvature from -pi to pi columns
    k2 = kymofluor{2}(allpatches,tlist{i})'; % actin channel, compare over all positions, same frames as vel
    %s = k2>0; % useful to avoid nans, compare only when normalized actin signal is > 0 
%     k3 = k3(s); 
%     k2 = k2(s);
    k2k3data = [k2, k3];
    n1 = length(patches{1});
    n2 = length(patches{2});
    b = [0,n1,n1+n2,2*n1+n2];
    bl = b+0.5;
    c = floor(nbins/4); 
    markticks = b+c; %sort([b,b+c],'ascend'); % b + c
    % matrix of correlation coefficient between pairs of columns of input matrix
    R = corr(k2k3data,'type',corrtype);

    ax1=imagesc(R); axis equal; hold on
    xline([bl(2),bl(4)],'LineWidth',1); xline(bl(3),'LineWidth',2.5)
    yline([bl(2),bl(4)],'LineWidth',1); yline(bl(3),'LineWidth',2.5)
    ax1.Parent.CLim = [-1,1]; colorbar;
    set(gca,'YDir','normal')
    ax1.Parent.YLim=[0.5, size(k2k3data,2)+0.5];
    ax1.Parent.XLim=[0.5, size(k2k3data,2)+0.5];
    xticks(markticks); yticks(markticks);
    if kymofixedtheta 
        xticklabels({'Actin: \theta=0','Actin: \theta=\pi','Memb: \theta=0'...
            ,'Memb: \theta=\pi'})
        yticklabels({'Actin: \theta=0','Actin: \theta=\pi','Memb: \theta=0'...
            ,'Memb: \theta=\pi'})
        if guvname(1) == 'L'
            xlabel({'','\theta = angular position relative to micropipette'})
        else
            thetac = num2str(wrapToPi(deg2rad(thetacenter)/pi));
            xlabel({'',['\theta = angular position relative to reference angle (',thetac,'\pi)']});
        end
    elseif kymofixedpts && guvname(1) == 'L'
        xticklabels({'Actin: Bins_{near}','Actin: Bins_{far}','Memb: Bins_{near}','Memb: Bins_{far}'})
        yticklabels({'Actin: Bins_{near}','Actin: Bins_{far}','Memb: Bins_{near}','Memb: Bins_{far}'})
    elseif kymofixedpts && guvname(1) == 'G'
        xticklabels({'Actin: Bins_{center}','Actin: Bins_{last}','Memb: Bins_{center}','Memb: Bins_{last}'})
        yticklabels({'Actin: Bins_{center}','Actin: Bins_{last}','Memb: Bins_{center}','Memb: Bins_{last}'})
    end

    %ax1.Parent.XLabel.FontName = 'Arial';

    ax1.Parent.XAxis.FontSize = 12;
    ax1.Parent.YAxis.FontSize = 12;
    ax1.Parent.XLabel.FontSize = 11; % must come after XAxis, bc its a subset of those

    if i == 1 
        titlelist = [guvname,': Pre-Rapamacyin'];
    elseif i ==2
        titlelist = [guvname,': Post-Rapamacyin'];
%     elseif i == 3
%         titlelist = [guvname,': All frames'];
    end

%     title({guvname,'Actin-curvature correlation coefficients'},...
%         {'computed for time-series of different angular positions'})
    % title ( main title with multiple lines, subtitle )
    title({titlelist,['Actin-membrane marker ',corrtype,' corr.']},...
        {'computed for time-series of different angular positions'})

    hold off
    
    f30.Position(3:4) = [700,600];
    %shading  flat; 
    drawnow;




    if newsavefigs && i == 1
        print(f30,strcat(newoutputfolder,guvname,'_',tstamp1,'_fig_actinmemb_preRapcorr_',corrtype(1:5)),'-dpng')

    elseif newsavefigs && i == 2
        print(f30,strcat(newoutputfolder,guvname,'_',tstamp1,'_fig_actinmemb_postRapcorr_',corrtype(1:5)),'-dpng')

%     elseif newsavefigs && i == 3
% 
%         print(f20,strcat(newoutputfolder,guvname,'_',tstamp1,'_fig_actincurv_fullcorr'),'-dpng')

    end

end




end



function [Uproj,Usvd,percentvar,allsvalues] = pcabases(Z,pcadim,varargin)

if ~isempty(varargin)
    varname = varargin{1};
end

% here use all frames, save none for testing

% if Z is > 2 dimensional matrix, 
% create a data matrix with vectorized data
if length(size(Z))>2
    Xu = Z(:);
    X = reshape(Xu,size(Z,1)*size(Z,2),size(Z,3));
else
    X = Z;
end

% find mean image in training images
meanframe = mean(X,2);

% figure
% imagesc(reshape(meanframe,size(U,1),size(U,2))); axis equal; colorbar
% if ~isempty(varname)
%     title(strcat('mean frame: ',varname))
% else
%     title('mean frame')
% end

%center all images 
cX = X - meanframe;

% SVD to compute eigendecomposition
% use "thin" SVD to get only eigenvectors with non-zero eigenvalues
[Usvd, Ssvd, Vsvd] = svd(cX,'econ');

%eigenvalues of covariance matrix = square of singular values of X
allsvalues = diag(Ssvd);

figure();
plot(allsvalues)
hold on
plot(pcadim,allsvalues(pcadim),'r*')
title('List of all non-zero singular values')
xlabel('index of singular values in decreasing order')
ylabel('singular value')

ax= gca; ax.FontSize = 14;
ax.Title.FontSize = 17;
ax.TitleFontWeight = "normal";

%percent variance explained is based on cumulative sum of eigenvalues (square of singluar values)
percentvar = 100*cumsum(allsvalues.^2)/sum(allsvalues.^2);

% PCA dimensionality (# of eigenvectors used to construct subspace)
N = [1:size(Ssvd,1)]';

% variance explained for each PCA dimensionality
varN = percentvar(N);

figure
plot(N,varN,'-o')
hold on
plot(pcadim,varN(pcadim),'r*')
title('Percent variance explained')
xlabel('# of eigenvectors (PCA dimensionality)')
ylabel('percent variance captured')

ax= gca; ax.FontSize = 14;
ax.Title.FontSize = 17;
ax.TitleFontWeight = "normal";

% compute projection matrix for each dimensionality
% use top N = pcadim eigenvectors as basis
% pcadim = 20; % based on percent variance explained 
Un_svd = Usvd(:,1:pcadim);
PM = Un_svd*Un_svd';
% project each frame 
[ptrain, Xtraincoord] = projectdata(cX,PM,Un_svd);

ptrain = ptrain + meanframe;
Uproj = reshape(ptrain,size(Z,1),size(Z,2),size(Z,3));


end

function [ptrain, Xtraincoord] = projectdata(cX,PM,Un)
% INPUT:
% cX: 4608 pixel x 150 training images/frames, centered
% PM: 4608 x 4608 projection matrix to PCA subspace
%     note: PM = Un*pinv(Un) = Un*Un' where Un is 4608 x N 
% Un: 4608 x N (N columns for top N eigenvectors)
% OUTPUT:
% ptrain: 4608 x 150, new frames after projection to PCA subspace

% Xtraincd: N x 150, new coordinates for projected training images


    Xtrain = cX; %4608 x 150
  
    n_train = size(cX,2);%150
    N = size(Un,2);
    Xtraincoord = zeros(N,n_train);
   

    % compute projections and coordinates of all training images
    ptrain = zeros(size(Xtrain)); 
    for i = 1:n_train 
        ptrain(:,i) = PM*Xtrain(:,i);
        Xtraincoord(:,i) = Un'*Xtrain(:,i);
    end

end


function varargout = v2struct(varargin)

%% Description
%    v2struct has dual functionality in packing & unpacking variables into structures and
%    vice versa, according to the syntax and inputs.
%
%    Function features:
%      * Pack variables to structure with enhanced field naming
%      * Pack and update variables in existing structure
%      * Unpack variables from structure with enhanced variable naming
%      * Unpack only specific fields in a structure to variables
%      * Unpack without over writing existing variables in work space
%
%    In addition to the obvious usage, this function could by highly useful for example in
%    working with a function with multiple inputs. Packing variables before the call to
%    the function, and unpacking it in the beginning of the function will make the function
%    call shorter, more readable, and you would not have to worry about arguments order any
%    more. Moreover you could leave the function as it is and you could pass same inputs to
%    multiple functions, each of which will use its designated arguments placed in the
%    structure.
%
%- Syntax
%    Pack
%      S = v2struct
%      S = v2struct(x,y,z,...)
%      S = v2struct(fieldNames)
%      S = v2struct(A,B,C,..., fieldNames)
%      S = v2struct(x,..., nameOfStruct2Update, fieldNames)
%      v2struct
%      v2struct(x,y,z,...)
%      v2struct(fieldNames)
%      v2struct(A,B,C,..., fieldNames)
%      v2struct(x,..., nameOfStruct2Update, fieldNames)
%
%    Unpack
%      v2struct(S)
%      [a,b,c,...] = v2struct(S)
%      v2struct(S,fieldNames)
%      [a,b,c,...] = v2struct(S,fieldNames)
%
%- Inputs & Outputs
%    Pack - inputs
%      x,y,z,...           - any variable to pack. can be replaced by fieldNames below.
%      nameOfStruct2Update - optional, name of structure to update if desired.
%      fieldNames          - optional, cell array of strings, which must include a cell
%                            with the string 'fieldNames' and must be the last input.
%    Pack - outputs 
%      S - the packed structure. If there is no output argument then a structure named
%          Sv2struct would be created in the caller workspace.
%
%    Unpack - inputs
%      S          - name of structure to be unpacked.
%      fieldNames - optional, cell array of strings, which must include a cell with the
%                   string 'fieldNames' and must be the last input.
%    Unpack - outputs          
%      a,b,c,... - variables upacked from the structure.
%                  if there are no output arguments then variables would be created in
%                  the caller workspace with naming according to name of inputs.
%
%- Examples
%  % see 'Usage example' section below for convenient presentation of these examples.
%
%    % NOTE: whenever using filedNames cell array please note the following
%    %  1. fieldNames cell array must include a cell with the string 'fieldNames'
%    %  2. fieldNames cell array input must be the last input.
%
%  % Pack
%      x = zeros(3); x2 = ones(3); y = 'Testing123'; z = cell(2,3);
%      fieldNames1 = {'fieldNames', 'x', 'y', 'z'};
%      fieldNames2 = {'fieldNames', 'a', 'b', 'c'};
%      fieldNames3 = {'fieldNames', 'x'};
%      nameOfStruct2Update = 'S';
%
%       % The four examples below return structure S with same values however the
%       % structure's field names are defined differently in every syntax. 
%      % Example 1.
%      % structure field names defined by variables names.
%       S = v2struct(x,y,z) 
%      % Example 2.
%      % structure field names defined according to the cell array fieldNames. 
%       % NOTE: variables with the names in fieldNames1 must exist in the caller workspace.
%       S = v2struct(fieldNames1) 
%      % Example 3.
%      % same as #1. but arguments are passed explicitly
%       S = v2struct(zeros(3), 'Testing123', cell(2,3), fieldNames1) 
%      % Example 4.
%      % field names defined by content of fieldNames2 while
%      % the values are set according to the passed arguments. In this case the structure
%      % S returned would be: S.a=x, S.b=y, S.c=z
%       S = v2struct(x,y,z, fieldNames2) 
%
%      % Example 5.
%      % update structure S. The fields that would be updated are according to content
%      % of fieldNames3. Note that you must pass a variable with the name
%      % 'nameOfStruct2Update' placed before 'fieldNames3'. This variable should contain
%      % the name of the structure you want to update as a string. Also note that if you
%      % set an output structure name which is different than the one stated in
%      % nameOfStruct2Update a new structure would be created and the structure that was
%      % meant to be updated would not get updated.
%       S.oldField = 'field to be saved for future use'
%       S = v2struct(x2, nameOfStruct2Update, fieldNames3)
%
%      % Example 6.
%      % pack all variables in caller workspace. Call without input arguments.
%        S = v2struct
%
%      % The following examples return the same results as the examples above but the
%      % structure would be returned with the default name 'Sv2struct'. Be cautious as
%      % this might lead to overriding of arguments.
%      % Example 7.
%       v2struct(x,y,z)
%      % Example 8.
%       v2struct(fieldNames1)
%      % Example 9.
%       v2struct(zeros(3), 'Testing123', cell(2,3), fieldNames1)
%      % Example 10.
%       v2struct(x,y,z, fieldNames2)
%      % Example 11.
%       S.oldField = 'field to be saved for future use'
%       v2struct(x2, nameOfStruct2Update, fieldNames3)
%      % Example 12.
%       v2struct
%
%  % Unpack
%      clear S x x2 y z fieldNames1 fieldNames2 fieldNames3 nameOfStruct2Update
%      S.x = zeros(3); S.y = 'Testing123'; S.z = cell(2,3);
%      fieldNames3 = {'fieldNames','x','z'};
%
%      % Example 1.
%      % This example creates or overwrites variables x, y, z in the caller with the
%      % contents of the corresponding named fields.
%       v2struct(S)
%
%      % Example 2.
%      % This example assigns the contents of the fields of the scalar structure
%      % S to the variables a,b,c rather than overwriting variables in the caller. If
%      % there are fewer output variables than there are fields in S, the remaining fields
%      % are not extracted.
%       [a,b,c] = v2struct(S)
%
%      % Example 3.
%      % This example creates or overwrites variables x and z in the caller with the
%      % contents of the corresponding named fields.
%       v2struct(S, fieldNames3)
%
%      % Example 4.
%      % This example assigns the contents of the fields 'x' and 'z' defined by
%      % fieldNames3 of the scalar structure S to the variables a and b rather than
%      % overwriting variables in the caller. If there are fewer output variables than
%      % there are fields in S, the remaining fields are not extracted.
%       [a,b] = v2struct(S, fieldNames3)
%
%       % This example unpacks variables 'y' and 'z' only without overwriting variable 'x'. 
%       % NOTE the addition of the field named 'avoidOverWrite' to the structure to be
%       % unpacked. This is mandatory in order to make this functionality work. The
%       % contents of this field can be anything, it does not matter. 
%      S.avoidOverWrite = '';
%      x = 'do not overwrite me';
%      v2struct(S)
%
%- Usage example (includes sub-functions)
%    1. run attached v2structDemo1.m file for on screen presentation of examples.
%    2. run attached v2structDemo2.m file and read comments in file for a suggestion of
%       how to use v2struct in managing input to other functions with improved usability.
%
%- Revision history
%    2011-05-19, Adi N., Creation
%    2011-05-29, Adi N., Update structure added, some documentation and demo function changes
%    2011-06-02, Adi N., Fixed updating structure functionality
%    2011-06-05, Adi N., Added functionality: avoid overwritring existing variables, added
%                        unpacking examples to demo1 .m file.
%    2011-06-30, Adi N., fieldNames usage corrected, now must include a specific string to
%                        be triggered. Documentation enhanced. Code tweaked.
%    2011-07-14, Adi N., Fixed bug in packing with variables only
%    2011-08-14, Adi N., Clarified warning and error when packing/unpacking with
%                        fieldNames.
%    2011-09-12, Adi N., Added easy packing of all variables in caller workspace (thanks 
%                        to Vesa Lehtinen for the suggestion), fixed bug in warning
%                        handling in packing case, edited comments.
%
%    Inspired by the function: mmv2truct - D.C. Hanselman, University of Maine, Orono, ME
%    04469 4/28/99, 9/29/99, renamed 10/19/99 Mastering MATLAB 5, Prentice Hall,
%    ISBN 0-13-858366-8


% parse input for field names
if isempty(varargin)
   gotCellArrayOfStrings = false;
   toUnpackRegular = false;
   toUnpackFieldNames = false;
   gotFieldNames = false;
else
   gotCellArrayOfStrings = iscellstr(varargin{end});
   
   toUnpackRegular = (nargin == 1) && isstruct(varargin{1});
   if toUnpackRegular
      fieldNames = fieldnames(varargin{1})';
      nFields = length(fieldNames);
   end
   
   gotFieldNames = gotCellArrayOfStrings & any(strcmpi(varargin{end},'fieldNames'));
   if gotFieldNames
      fieldNamesRaw = varargin{end};
      % indices of cells with actual field names, leaving out the index to 'fieldNames' cell.
      indFieldNames = ~strcmpi(fieldNamesRaw,'fieldNames');
      fieldNames = fieldNamesRaw(indFieldNames);
      nFields = length(fieldNames);
   end
   toUnpackFieldNames = (nargin == 2) && isstruct(varargin{1}) && gotFieldNames;
end


% Unpack
if toUnpackRegular || toUnpackFieldNames 
   
   struct = varargin{1};
   assert(isequal(length(struct),1) , 'Single input nust be a scalar structure.');
   CallerWS = evalin('caller','whos'); % arguments in caller work space
   
   % update fieldNames according to 'avoidOverWrite' flag field.
   if isfield(struct,'avoidOverWrite')
      indFieldNames = ~ismember(fieldNames,{CallerWS(:).name,'avoidOverWrite'});
      fieldNames = fieldNames(indFieldNames);
      nFields = length(fieldNames);
   end
   
   if toUnpackRegular % Unpack with regular fields order
      if nargout == 0 % assign in caller
         for iField = 1:nFields
            assignin('caller',fieldNames{iField},struct.(fieldNames{iField}));
         end
      else % dump into variables
         for iField = 1:nargout
            varargout{iField} = struct.(fieldNames{iField});
         end
      end
      
   elseif toUnpackFieldNames % Unpack with fields according to fieldNames
      if nargout == 0 % assign in caller, by comparing fields to fieldNames
         for iField = 1:nFields
            assignin('caller',fieldNames{iField},struct.(fieldNames{iField}));
         end
      else % dump into variables
         assert( isequal(nFields, nargout) , ['Number of output arguments',...
            ' does not match number of field names in cell array']);
         for iField = 1:nFields
            varargout{iField} = struct.(fieldNames{iField});
         end
      end
   end
   
% Pack   
else
   % build cell array of input names
   CallerWS = evalin('caller','whos');
   inputNames = cell(1,nargin);
   for iArgin = 1:nargin
      inputNames{iArgin} = inputname(iArgin);
   end
   nInputs = length(inputNames);
   
      % look for 'nameOfStruct2Update' variable and get the structure name
   if ~any(strcmpi(inputNames,'nameOfStruct2Update')) % no nameOfStruct2Update
      nameStructArgFound = false;
      validVarargin = varargin;
   else % nameOfStruct2Update found
      nameStructArgFound = true;
      nameStructArgLoc = strcmp(inputNames,'nameOfStruct2Update');
      nameOfStruct2Update = varargin{nameStructArgLoc};
      % valid varargin with just the inputs to pack and fieldNames if exists
      validVarargin = varargin(~strcmpi(inputNames,'nameOfStruct2Update'));
      % valid inputNames with just the inputs name to pack and fieldNames if exists
      inputNames = inputNames(~strcmpi(inputNames,'nameOfStruct2Update'));
      nInputs = length(inputNames);
      % copy structure from caller workspace to enable its updating
      if ismember(nameOfStruct2Update,{CallerWS(:).name}) % verify existance
         S = evalin('caller',nameOfStruct2Update);
      else
         error(['Bad input. Structure named ''',nameOfStruct2Update,...
            ''' was not found in workspace'])
      end
   end
   
   % when there is no input or the input is only variables and perhaps
   % also nameOfStruct2Update   
   if ~gotFieldNames
      % no input, pack all of variables in caller workspace
      if isequal(nInputs, 0)
         for iVar = 1:length(CallerWS)
            S.(CallerWS(iVar).name) = evalin('caller',CallerWS(iVar).name);
         end
         % got input, check input names and pack
      else
         for iInput = 1:nInputs
            if gotCellArrayOfStrings % called with a cell array of strings
               errMsg = sprintf(['Bad input in cell array of strings.'...
                           '\nIf you want to pack (or unpack) using this cell array as'...
                           ' designated names'...
                           '\nof the structure''s fields, add a cell with the string'...
                           ' ''fieldNames'' to it.']);
            else
               errMsg = sprintf(['Bad input in argument no. ', int2str(iArgin),...
                                 ' - explicit argument.\n'...
                          'Explicit arguments can only be called along with a matching'...
                          '\n''fieldNames'' cell array of strings.']);
            end
            assert( ~isempty(inputNames{iInput}), errMsg);
            S.(inputNames{iInput}) = validVarargin{iInput};
         end
         
         % issue warning for possible wrong usage when packing with an input of cell array of
         % strings as the last input without it containing the string 'fieldNames'.
         if gotCellArrayOfStrings
            name = inputNames{end};
            % input contains structure and a cell array of strings
            if (nargin == 2) && isstruct(varargin{1})
               msgStr = [inputNames{1},''' and ''',inputNames{2},''' were'];
               % input contains any arguments with an implicit cell array of strings
            else
               msgStr = [name, ''' was'];
            end
            warnMsg = ['V2STRUCT - ''%s packed in the structure.'...
               '\nTo avoid this warning do not put ''%s'' as last v2struct input.'...
               '\nIf you want to pack (or unpack) using ''%s'' as designated names'...
               ' of the'...
               '\nstructure''s fields, add a cell with the string ''fieldNames'' to'...
               ' ''%s''.'];
            fprintf('\n')
            warning('MATLAB:V2STRUCT:cellArrayOfStringNotFieldNames',warnMsg,msgStr,...
                     name,name,name)
         end
      end
   % fieldNames cell array exists in input
   elseif gotFieldNames
      nVarToPack = length(varargin)-1-double(nameStructArgFound);
      if nVarToPack == 0 % no variables to pack
         for iField = 1:nFields
            S.(fieldNames{iField}) = evalin('caller',fieldNames{iField});
         end
         
         % else - variables to pack exist
         % check for correct number of fields vs. variables to pack
      elseif ~isequal(nFields,nVarToPack)
         error(['Bad input. Number of strings in fieldNames does not match',...
                'number of input arguments for packing.'])
      else
         for iField = 1:nFields
            S.(fieldNames{iField}) = validVarargin{iField};
         end
      end
      
   end % if ~gotFieldNames

if nargout == 0
   assignin( 'caller', 'Sv2struct',S );
else
   varargout{1} = S;
end

end % if nargin
end
