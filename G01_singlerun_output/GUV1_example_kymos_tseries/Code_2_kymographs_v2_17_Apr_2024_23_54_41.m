format compact
clearvars, clc
close all
warning('off','all')
set(0,'DefaultFigureWindowStyle','normal') 
set(0,'DefaultAxesTitleFontWeight','normal');
set(0,'DefaultAxesFontName','Arial')

%% IMPORTANT:
% check all sections in the code marked by "**SETUP**" (search the code)


%% **SETUP** type of analysis
runtype = 1; % 1: run kymograph analysis on a single target; can use just a few frames for quick checks
             % 2: run full kymograph analysis on 1 or more targets from a table
             % 3: find and load previous data and just adjust plots
             % 4: manually select file and adjust

savefigs  = 2;   % keep as 1 or 2 to save figures and movies
             % 0 = no saving; 
             % 1 = save mfiles,kymos,figs ; 
             % 2 = same as 1 + save dataRaw(images,other boundary quants etc)

 % part of the name of output folder - ALWAYS HAVE BACKSLASH
 datatype = 'GUV1_example_kymos_tseries/'; 

%%
if runtype == 1 
    %% **SETUP** Option 1: Adjust these for analysis on single target  
        targname = 'G01';
        listFrames = 1:20; % frames to use for kymograph analysis
        t_rapa = 7; % last frame before external input
        nsize = 1; % controls inner/outer masks that determine membrane region of a target;
                    % decrease for smaller targets, increase for larger; 
                    % check the time-series of inner/outer masks using variable: do_plotMaskseries
        thetacenter = 0; %  angle (degrees rel. to x-axis); a reference point to align center of kymograph e.g. for needle experiments

    % NOTE: once the above settings are good, can save them in a table 
    % see target_settings.csv which has 4 columns corresponding to above properties
    % see runtype == 2 (next section) for running analysis from a table

    runall(runtype,datatype,targname,t_rapa,thetacenter,listFrames,nsize,savefigs);

elseif runtype == 2 
    %% **SETUP** For analysis on multipe targets (cells, guvs, etc), use this table defining the key settings for each object

    % read the table of target settings (.csv file)
    Targ_set = readtable('target_settings.csv');
        
    % choose which targets to analyze
    % these correspond to rows in the table
    targ_idx_analysis = [1];%1,5,9,12,13,15:17
    n_analysis = length(targ_idx_analysis); % # of targets to be analyzed
    
    for i = 1:n_analysis 
        close all;
        ii = targ_idx_analysis(i);
        targname = Targ_set.targname{ii};
        lastframe = Targ_set.lastframe(ii);
        listFrames = 1:lastframe;
        t_rapa = Targ_set.t_rapa(ii);
        nsize = Targ_set.nsize(ii);
        thetacenter = Targ_set.thetacenter(ii);
        runall(runtype,datatype,targname,t_rapa,thetacenter,listFrames,nsize,savefigs);
    end

elseif runtype == 3 % load previously saved data

    % list of ALL folder
    Sfolds = dir('G*output'); % list of Global output folders

    df_list = cell(n_analysis,1); % list of .mat files for SELECT target
    %df = cell(n_analysis,1); % variables in each .mat file 
    tstampnew = datestr(now,'dd_mmm_yyyy_HH_MM_SS');

    for i = 1:n_analysis
        close all;
        tmpfolder = fullfile(Sfolds(targ_idx_analysis(i)).name,datatype);
        tmpfile = dir([tmpfolder,'*.mat']);
        df_list{i} = fullfile(tmpfolder,tmpfile.name);
        load(df_list{i},'dataRaw','dataParams','dataKymoBin','dataKymoVar');

        % adjust this to ensure internal saving of png within function goes to correct folder

        dataParams.outputfolder = tmpfolder;
               
        kf = dataParams.kymofixed;

        %% add mask timeseries 
        % input variables
        tstampnewtext =  strcat(tstampnew,'_MASKTS');
        dataParams.tstamp1 = tstampnewtext;
        tmpKymo = dataKymoBin;
        makekymo = dataParams.makekymo;
        % plot movies and save only based on the first kymograph set (fixedpts or fixedtheta usually)
        plottseries(dataParams,dataRaw,makekymo,tmpKymo,kf(1));

        % ADJUST FONT SIZE
%         tstampnewtext =  strcat(tstampnew,'_FONT14');
%         dataParams.tstamp1 = tstampnewtext;
% 
%         count = 0; % png saving step depends on count
%         adjfix(1) = plotphyskymos(dataParams,dataKymoBin,kf(1),count);
%         adjfix(2) = plotbiochemkymos(dataParams,dataKymoBin,kf(1),count);
%         count = 5;
%         adjvar(1) = plotphyskymos(dataParams,dataKymoVar,kf(2),count);
%         adjvar(2) = plotbiochemkymos(dataParams,dataKymoVar,kf(2),count);
%         if savefigs > 0
%             %savefig(f01,strcat(dataParams.outputfolder,dataParams.cellname,'_',tstampnewtext,'_NEW_macrosummary.fig'),'compact');
%             savefig(adjfix,strcat(dataParams.outputfolder,dataParams.cellname,'_',tstampnewtext,'_figsoutFix.fig'),'compact');
%             savefig(adjvar,strcat(dataParams.outputfolder,dataParams.cellname,'_',tstampnewtext,'_figsoutVar.fig'),'compact');
%         end


        % ADJUST BINARIZATION FOR PATTERN VISUALIZATION
%         % adjust time stamp to include adjustment
%         tstampnewtext =  strcat(tstampnew,'_PATTERNS');
%         dataParams.tstamp1 = tstampnewtext;
% 
%         count = 10;
% 
%         kmemb = dataKymoBin.kymofluor{1};
%         kacta = dataKymoBin.kymofluor{3}; % acta
%         kactin = dataKymoBin.kymofluor{2}; % actin
%         kcurv = dataKymoBin.kymocurv;
% 
%         actinthresh = median(kactin(:));
% 
%         scrv = sort(abs(kcurv(:)),'descend','MissingPlacement','last');
%         ncrv = length(scrv);
%         crvthresh = scrv(round((0.10)*ncrv)); % top 20%
% 
%         kactin1 = (kactin>=actinthresh); % threshold value different for each target
%         kcurv1 = (abs(kcurv)>=crvthresh);
%         %figure; imagesc(kactin1);
% 
%         %dataKymoBin.kymofluor{2} = kactin1;
%         dataKymoBin.kymocurv = kcurv1;
%         
%         
%         figsoutFix(1) = plotphyskymos(dataParams,dataKymoBin,kf(1),count);
% %         clrs = [0 0.4470 0.7410; 0.9290 0.6940 0.1250];
% %         ax = gca; colorbar off
% %         ax.CLim = [0,1];
% %         colormap(ax,clrs);
% %         colorbar
% 
% %        figsoutFix(2) = plotbiochemkymos(dataParams,dataKymoBin,kf(1),count);    
% %         ax = gca; colorbar off
% %         ax.CLim = [0,1];
% %         colormap(ax,clrs);
% %         colorbar
% 
%         count = 20;
% 
%         kactinvar = dataKymoVar.kymofluor{2}; % actin
%         kcurvvar = dataKymoVar.kymocurv;
% 
%         kactinvar1 = (kactinvar>=actinthresh)&(~isnan(kactinvar)); % threshold value different for each target
%         kcurvvar1 = (abs(kcurvvar)>=crvthresh)&(~isnan(kcurvvar));
%         %figure; imagesc(kactin1);
% 
%         kactinVnull = isnan(kactinvar);
%         kactinVfin = double(kactinvar1);
%         kactinVfin(kactinVnull) = nan;
%         
%         kcurvVnull = isnan(kcurvvar);
%         kcurvVfin = double(kcurvvar1);
%         kcurvVfin(kcurvVnull) = nan;
% 
% 
% %        dataKymoVar.kymofluor{2} = kactinVfin;
%         dataKymoVar.kymocurv = kcurvVfin;
% 
%         figsoutVar(1) = plotphyskymos(dataParams,dataKymoVar,kf(2),count);
% %         ax = gca; colorbar off
% %         ax.CLim = [0,1];
% %         colormap(ax,clrs);
% %         colorbar
% %        figsoutVar(2) = plotbiochemkymos(dataParams,dataKymoVar,kf(2),count);
% %         ax = gca; colorbar off
% %         ax.CLim = [0,1];
% %         colormap(ax,clrs);
% %         colorbar
%         if savefigs > 0
%             %savefig(f01,strcat(dataParams.outputfolder,dataParams.targname,'_',tstampnewtext,'_NEW_macrosummary.fig'),'compact');
%             savefig(figsoutFix,strcat(dataParams.outputfolder,dataParams.targname,'_',tstampnewtext,'_figsoutFix.fig'),'compact');
%             savefig(figsoutVar,strcat(dataParams.outputfolder,dataParams.targname,'_',tstampnewtext,'_figsoutVar.fig'),'compact');
%         end


    end

elseif runtype == 4 % chose manually
%     close all;
%     tstampnew = datestr(now,'dd_mmm_yyyy_HH_MM_SS');
% 
%     [filelist,filepath]=uigetfile('*.mat','MultiSelect','on',...
%      'Select LIST to plot'); pause(0.5); 
%     
%     tmpfullfile = fullfile(filepath,filelist);
%     load(tmpfullfile,'dataRaw','dataParams');
%     
%     dataParams.do_fixedbox = 1;
% 
%     % adjust this to ensure internal saving of png within function goes to correct folder
%     dataParams.outputfolder = filepath;
%     dataParams.tstamp1 = strcat(tstampnew,'_NEW');
% 
% 
%     f01 = plotmovement(dataParams,dataRaw);
% 
%     if savefigs == 2
%         savefig(f01,strcat(dataParams.outputfolder,dataParams.targname,'_',tstampnew,'_NEW_macrosummary.fig'),'compact');
%     elseif savefigs == 1
%         savefig(f01,strcat(dataParams.outputfolder,dataParams.targname,'_',tstampnew,'_NEW_macrosummary.fig'),'compact');
%     end



end



function runall(runtype,datatype,targname,t_rapa,thetacenter,listFrames,nsize,savefigs)
%--------------------------------------------------------------------------  
% SET ITYPE IN INTERPBOUNDARY! 'pchip' or 'makima'
% SET HOW FLUOR INTENSITIES ARE COMPUTED!
% SET kymograph smoothing window
% targname         = 'G09';
% t_rapa          = 7;       % last frame before rapamycin added, approx = time 0 for recruitment
% thetacenter    = 0;      % angle (degrees rel. to x-axis) reference point to align center of kymograph
% listFrames      = [1:74];   % choose frames (not necessarily consecutive) to include in analysis, in ascending order
% %IMPORTANT FOR SETTING SIZE OF INNER OUTER MASKS!!!
% if strcmp(targname(1),"G"); nsize = 0.7; % 1 for G01, 
% else; nsize = 3; % 3 for L48, check tseries inner/outer masks 
% end

%% **SETUP** most settings can stay constant; 
% adjust the PLOTTING OPTIONS section to control number of plots

% PLOTTING OPTIONS: 0/false or 1/true
do_plotAlignseries   = 1;  % display boundary, centroid, kymo center point
do_plotTimeseries    = 1;    % display MAIN time series : FLUORESC
do_plotMaskseries    = 1;    % plot inner, outer mask time series
do_plotBoundaryOverlay = 0; % separate plot of boundary overlay over time 
do_fixedbox          = 0;  % in summary fig, plot with fixed box size for centroid and boundary evolution
do_plotVelseries_pdist2     = 0;   % keep 0; slower and less accurate; displacements with velocity in color (pdist2 based)
do_plotVelseries_motiontrack = 1; % plotting displacements with velocity in color (more accurate, motion-track based)
do_plotEllipseries   = 1;  % plot bounding ellipse ; if = 1, do_aspratio_ind must = 1
do_aspratio_ind     = 1;  % 1 = do independent asp ratio computation, 
                          % 0 = use matlab region props maj/min axis to compute aspect ratio

% more detailed settings below:                      
basystem = 'mac';
kymofixed = [2,0];  
    % options for the two kymographs generated; first number is for the 1st kymo, etc
    % for first kymograph
    % 1: fixedpts:  use fixed # of boundary/perimeter points for 1st kymograph (= # pts in longest perim, post-interp)
    % 2: fixedtheta:  % use fixed # of angular bins for 1st kymograph
    % for second kymograph - always use 0 (varying number of boundary points)
    % 0: kymograph with varying number of boundary points - useful for checking cells/guvs that change size;

nbins   = 360; % # of bins to use for 1st kymo ( kymofixedpts or kymofixedtheta)
pvary   = 3; % multiple of perimeter for final interpolation of varying kymo

m0      = 3; % boundary points per pixel; will multiply perimeter, used in pre-computation interp
m1      = 6; % muliple of perimeter, post-computation interp - needed for small targets to get enough spatial bins

% see getIntensity function for how intensity at the membrane is computed
% an inner and outer boundary is determined for the boundary of the target
plperp     = 0.6; % 60% outward from line perpendicular to boundary (from inner towards outer boundary) for intensity computation
do_top3px  = 0; % (set this to 0 to use the mean; or 1 to use only top 3 px for fluor intensity)
do_lumennorm = 0; % 0/false, or 1/true: divide biochem signals by lumen signal (use 0 usually)
do_membnorm = 0; % 0/false, or 1/true: divide biochem signals by membrane marker signal (use 0 usually)

veltrack1  = 1; % compute velocity using motion tracking (1) or pdist2 (0) for fixed kymo
veltrack2  = 1; % compute velocity using motion tracking (1) or pdist2 (0) for variable kymo
    % alignvel will be based on whatever is chosen here. bdyvel is still pdist2 based

% ---Select curvature calculation method
computecurv = 2; % 1 for circumcenter; 2 for fitcircle

% ----Select smoothing
savgolay_smooth = true;  % do savitzky-golay smoothing of boundary for curvature estimates
darc      = 15;  % distance (pixels) of arc segment used for savgolay 
darc2     = 15;  % distance (pixels) of arc segment used for curvature computation
% parc           = 0.05;  % prop of perimeter: arc segment used for savgolay 
% parc2          = 0.05;  % prop of perimeter: arc segment used for curvature computation
do_movavg  = 0; % 0/false, or 1/true: do moving average
s_movavg   = 2*m1; % if do_movavg, use this for # of points right/left 

%------------------------------------------
minfont = 20;
sat = 1; % factor for caxis upper limit fluor signals
satlow = 0; % factor for caxis lower limit of fluor signal


%% **SETUP** Fluorescent marker labels and Folder settings 
%--------------------------------
% make sure this matches where the tiff files for the smoothed images and
% masks are located
inputfolder     = strcat(targname,'_input/',targname,'_acs_v0/');
%---------------------
if runtype == 1
makekymo     = 1;    % create kymograph
outputfolder    = strcat(targname,'_singlerun_output/',datatype); % name of output folder
elseif runtype == 2
makekymo     = 1;    % Create kymograph
outputfolder    = strcat(targname,'_tablerun_output/',datatype); % name of output folder    
else % only re-runs
makekymo    = 0;    % No kymograph
outputfolder    = strcat(targname,'_output/',dataype); 
end

% % name of files in the input folder - label fluor channels 
fname1          = strcat(targname,'_mCh_memb'); % membrane marker file name
fname2          = strcat(targname,'_YFP_actin'); % actin marker file name
fname3          = strcat(targname,'_CFP_ActA');
fnames          = {fname1;fname2;fname3};
% if strcmp(targname(1),'L')
%     fname1          = strcat(targname,'_mCh_memb'); % membrane marker file name
%     fname2          = strcat(targname,'_YFP_actin'); % actin marker file name
%     fname3          = strcat(targname,'_CFP_ActA');
%     fname4          = strcat(targname,'_BF');
%     fname5          = strcat(targname,'_647_dye');
%     fnames          = {fname1;fname2;fname3;fname4;fname5};
%     needlech            = 4 ; % BF channel for needle
% elseif strcmp(targname(1),'G')
%     fname1          = strcat(targname,'_mCh_memb'); % membrane marker file name
%     fname2          = strcat(targname,'_YFP_actin'); % actin marker file name
%     fname3          = strcat(targname,'_CFP_ActA');
%     fnames          = {fname1;fname2;fname3};
% end

maskfname       = strcat(targname,'_mask_inclusion');
membch       = 1 ; % channel of membrane marker
fluorch         = [1:3]; % membrane 1, actin 2, ActA 3;
nfluor          = length(fluorch);
nchannels        = length(fnames);
fullFileNames   = cell(nchannels,1);
for i = 1:nchannels
    fullFileNames{i}  = strcat(inputfolder,fnames{i},'_final.tif');
%     if ~isfile(fullFileNames{i})
%         fullFileNames{i} = strcat(inputfolder,fnames{i},'.tif');
%     end
end
tstamp1 = datestr(now,'dd_mmm_yyyy_HH_MM_SS');
finfo        = imfinfo(fullFileNames{1});
nallframes      = length(finfo);
cols         = finfo.Width;
rows         = finfo.Height;
bits         = finfo(1).BitDepth;
im0          = listFrames(1);        
imL          = min(nallframes,listFrames(end));   
deltaf     = imL - im0 + 1; % useful when skipping frames during testing
nframes      = min(nallframes,length(listFrames));
dt           = 1; % keep 1 for pixels/frame; use seconds per frame for pixels/sec 
% dx           = 1/finfo.Xresolution; %1/2.2764 = 0.44 microns/pixel
if exist(outputfolder,'dir') == 0
    mkdir(outputfolder);
end
%---------------------------
%ticks for kymographs
if nframes < 10
    tickstep = 1;
else
    tickstep = floor(nframes/10);
end
%--------------------------------------------------------------------------   

maskFileName = strcat(inputfolder,maskfname,'_final.tif');
if ~isfile(maskFileName)
    maskFileName = strcat(inputfolder,maskfname,'.tif');
end
maskInfo   = imfinfo(maskFileName);
if nallframes~=length(maskInfo)
    disp('Error: Number of frames in the image and mask files must be the same')
end
% save initial settings
dataParams = v2struct;


%----------------------------------------------------------
% cell arrays that are not channel specific
allCurvatures  = cell(nframes,1);                    % Curvatures, all frames
allpdist2Velocities  = cell(nframes,1);                    % Boundary velocity, all fr.
allpdist2Dcumuls  = cell(nframes,1);
allBoundaries  = cell(nframes,1);                    % Boundary positions, all fr.

% channel specific
allImages = cell(nframes,nchannels);   % original images, frame i, channel j
Igs     = cell(nframes, nchannels);
Iinit   = cell(nframes, nchannels);
allIrgb = cell(nframes, nchannels);
allIntensities = cell(nframes,nfluor);    % boundary intensities for all frames and fluorescent channels
alliBW  = cell(nframes,nfluor);
alloBW  = cell(nframes,nfluor);
allinner = cell(nframes,nfluor);
allouter = cell(nframes,nfluor);

% numeric arrays
allCentroids   = zeros(nframes,2);                   % Centroid, all frames
allAspRatio_ml  = zeros(nframes,1);   % compute using matlab ratio of major/minor axis length from region props
allAspRatio_ind      = zeros(nframes,1);     % independent computation of Aspect ratio (ratio of major/minor axis of minimum bounding ellipse)
allA0          = cell(nframes,1);  % independent computation allows for plotting the best fit ellipse
allA0c         = cell(nframes,1);  % independent computation allows for plotting the best fit ellipse
allEccents     = zeros(nframes,1);    %Eccentricity
allMajAxL      = zeros(nframes,1);                   % Major axis length
allLengths     = zeros(nframes,1);                   % Length of boundaries
allPerims      = zeros(nframes,1);
allAreas       = zeros(nframes,1);                   % target areas, all frames
allBW          = zeros(rows,cols,nframes,'logical'); % Segmentation all frames - only based on 1 channel


%% compute boundary points and smooth
%----------------------------------------------------
% This section does velocity/boundary computations for membrane channel
% and intensity computations for every channel
    
for i = 1:nframes

    movieFrame=listFrames(i); %im0:imL;  % TIFF indices   
    frame = i; %movieFrame-im0+1; % indices for saving data
    fprintf('Frame number: %4d. ',movieFrame);

    for k = 1:nchannels
        %---------------------------------------------
        % gaussian smoothing for membrane (I) and actin (I2)
        [Xgs,Xorig] = readframe(fullFileNames{k},movieFrame); %uint16 output (not rgb)    
        allImages{frame,k} = Xorig;
        Igs{frame,k} = Xgs;


        % define initial images for computations
        % segment based on I (gaussian smoothed membrane maker) 
        % or segment based on Iorig (original membrane marker)
        Xinit = Xgs;
        Iinit{frame,k} = Xinit; 

        % make boundary marker into rgb image for plotting
        % use imadjust to increase contrast
        Xrgb = repmat(imadjust(Xinit),1,1,3);
        allIrgb{frame,k} = Xrgb;
    end 

    BW = logical(imread(maskFileName,movieFrame));
    [B,L,ncells]     = bwboundaries(BW,'noholes');
    stats     = regionprops(L,'Area','Centroid','Perimeter','Eccentricity','MajorAxisLength','MinorAxisLength');

    % identify boundary points for single object (exclude new ones) 
    if ncells==1
         % BOUNDARY FROM MASK
        bdyraw      = B{1}; % (row,col) coordinates (y,x)
        area     = stats.Area;
        centroid = stats.Centroid;
        perimeter = stats.Perimeter;
        eccent   = stats.Eccentricity;
        mjaxl   = stats.MajorAxisLength;
        aspratio_ml = stats.MajorAxisLength/stats.MinorAxisLength;

        % interpolate boundary to have spacing of 1 pixel
        bdy = interpboundary(bdyraw,m0*perimeter+1);
    elseif ncells > 1 && frame > 1
        tmpcenters = zeros(ncells,2);
        tmpdist    = zeros(ncells,1);
        for k=1:ncells
            tmpcenters(k,:) = stats(k).Centroid;
            tmpdist(k)      = norm(tmpcenters(k,:)-allCentroids(frame-1,:));
        end
        [~,kstar] = min(tmpdist);
        % BOUNDARY FROM MASK
        bdyraw       = B{kstar}; % (row,col) coordinates (y,x)
        area      = stats(kstar).Area;
        centroid  = stats(kstar).Centroid;
        perimeter = stats(kstar).Perimeter;
        eccent = stats(kstar).Eccentricity;
        mjaxl   = stats(kstar).MajorAxisLength;
        aspratio_ml = stats(kstar).MajorAxisLength/stats(kstar).MinorAxisLength;
        % interpolate boundary to have spacing of 1 pixel; 
        bdy = interpboundary(bdyraw,m0*perimeter+1); % #o unique points + 1
    end


    %---------------------------------------------

    % optional savitzky-golay smoothing of boundary
    % Now smooth with a Savitzky-Golay sliding polynomial filter
    %deltaarc1 = parc*perimeter*m0; % #points in arc = perimeter fraction * points per pixel (m0)
    deltaarc1 = darc*m0; % #points in arc = arc distance (pixels) * points per pixel (m0)
    windowWidth = 2*(round(0.5*deltaarc1))+1; % must be odd; adjust based on length of interpBdy
    polynomialOrder = 2; % keep this as 2 for quadratic smoothing of boundary
    % use membrane channel rgb to display
    bdy_sg = sgsmoothboundary(allIrgb{frame,membch},bdy,windowWidth,polynomialOrder);


    %--------------------------------------
    if movieFrame > im0
        oBorigs = Borigs; % save old boundary for velocity before updating
    end
    %----------------------------------------

    % ============ SELECT boundary to use
    if savgolay_smooth
        Borigs = bdy_sg; % savitzky-golay smoothed boundary
        % sgolay smoothed boundary will only be approximately circular 
        % (not closed loop usually due to edge effects)
    else
        Borigs = bdy; % unsmoothed boundary
    end


    %----------------------------------------
    %% compute aspect ratio
    %prm0 = v2struct(do_plotEllipseries,movieFrame,minfont);
    if do_aspratio_ind
        [aspratio_ind,A0,A0c,~] = aspectratio(Borigs);
    end

    %% needle position check
    

    %% compute velocities (pixels/frame)
    

    if movieFrame > im0  %&& ~veltrack

        % Get displacement between old boundary and new boundary over 1 frame
        % (row,col) coordinates (y,x)
        % dt = 1 so it's same as displacements
        pVel = v2struct(do_plotVelseries_pdist2,savefigs,movieFrame,im0,outputfolder,...
            targname,tstamp1,minfont);
        [vel,idx] = velocity_pdist2(Borigs,oBorigs,dt,pVel,allIrgb{frame,membch});    
        dcumul = dcumul(idx) + dt*vel;
        %velold = vel;
    else
        vel = nan(length(Borigs),1);
        %velold = zeros(length(Borigs),1);
        dcumul = zeros(length(Borigs),1);
    end

    %% compute intensity from all channels
    % At the ith frame
    Iinit_curr = Iinit(frame,:); % cell array of all channels at current frame

    pFluor = v2struct(centroid,movieFrame,im0,imL,do_plotTimeseries,do_lumennorm,do_membnorm,do_top3px,savefigs,outputfolder,...
        tstamp1,fnames,fluorch,nfluor,nchannels,nsize,plperp,bits);
    [fluor,iBW,oBW,inner,outer] = getIntensity(Borigs,Iinit_curr,BW,pFluor); % get actin intensity +/- boundary 


    % ---------------


    %% compute k-circumcenter

%     % =========== CIRCUMCENTER
%     % circumcenter approach to estimate curvature 
% 
%     % define spaced boundary points for curvature estimation 
%     %defBk = 'distance';   % define distance (in pixels) separating boundary points
%     defBk = 'number';     % define number of boundary points
%     %defBk= 'original';   % use original boundary points
%     
%     Bk_dist = parc2*perimeter*m0;    %#points spaced out in arc = perimeter fraction * points per pixel (m0)
%     Bk_num = round(perimeter/deltaarc1);  % number of unique (non-circular) spaced points around boundary
%     curvatureThresh = 1;     %0.06 the maximum allowed value of the curvature measure
%     loopclose = 1;              % 0 - if open boundaries | 1 - if closed boundaries
%     iflag = 2;    % two curvature interpolation options for circumcenter approach 
%     % needed b/c approach uses spaced boundary points
%     % 1 use scatteredinterpolant to match coordinates of bdy ; 
%     % 2 use interpboundary to match length of bdy
% 
%     % spaced boundary points for estimating k (curvature)
%     % Bkpoints is circular
%     Bkpoints = defBkpoints(I1rgb,Borigs,defBk,Bk_dist,Bk_num);
% 
%     Pcurv1 = v2struct(curvatureThresh,loopclose,perimeter);
%     [Bkcurv0] = curvature1_v5(Bkpoints,Pcurv1);
%     % Bkcurv0 is NOT circular;  Borigs,Bkpoints are circular
%     % make  Bkcurv circular
%     Bkcurv = [Bkcurv0;Bkcurv0(1)];
%    
%     % interpolate curvature to match length or coordinates of bdy (iflag = 1 or 2) 
%     % and plot overlay as fig 21 
%     curv1 = plotcurvature1_v6(I1rgb,Borigs,Bkcurv,Bkpoints,iflag);
%     % -------------
%% compute k-fitcircle
    % ====================== FITCIRCLE
    % fitcircle approach to estimate curvature
    %deltaarc2 = parc2*perimeter*m0; % 0.1, %#points in arc = perimeter fraction * points per pixel (m0)
    deltaarc2 = darc2*m0; % 0.1, %#points in arc = arc distance (pixels) * points per pixel (m0)
    deltaarc2_half=round(0.5*deltaarc2); % 2*delta+1 is used for symmetry
    c2 = curve_fitcircle(Borigs,deltaarc2_half); % use savgolay smoothed or raw boundary
    % Borigs is the selected boundary for circumcenter/fitcircle = bdy or bdy_sg
    curv2 = plotcurvature2_v1(allIrgb{frame,membch},Borigs,c2); % plot on original boundary

    % test 
    %c2 = curve_fitcircle(Bkpoints,1); % use savgolay smoothed or raw boundary
    %curv2 = plotcurvature2_v1(I1rgb,bdy,c2,Bkpoints,iflag);
    
    
    % Select which curvature approach to use
    if computecurv == 1
        curv = curv1; % use circumcenter approach
        % show boundary spacing, smoothing used for circumcenter 
        % f1 = displayboundary(I1rgb,bdy,Borigs,Bkpoints);
    else
        curv = curv2; % use fitcircle approach
    end


    
    %% post computation interpolation of all boundary measurements
    

    n   = round(m1*perimeter)+1;%  set n to be: (number of new UNIQUE boundary points) + 1
    % if want spacing of 0.5 pixel unit, set # of unique boundary points = 2*perimeter
    % n = 2*perimeter + 1


    % gfp, vel only approximately circular; bdy, curv are circ

    iFluor = cell(1,nfluor); % pre-allocation is required
    % interpolate bdy quantities
    % pass fluor as comma separated list, not as a cell array
    [iBdy,iVel,iCurv,iDcumul,iFluor{:}] = interpboundary(Borigs,n,vel,curv,dcumul,fluor{:});
    % interpBdy and interpCurv are circular; others just very close
    % make circular all
    iVel(end) = iVel(1);
    iCurv(end) = iCurv(1);
    iDcumul(end)= iDcumul(1);
    matiFluor = cell2mat(iFluor);
    matiFluor(end,:) = matiFluor(1,:);
    iFluor = mat2cell(matiFluor,length(iFluor{1}),ones(1,nfluor));

    if savefigs > 0 % save m file
        if movieFrame == imL
            % save mfile if saving figs 
            currentfile = strcat(mfilename,'.m');
            destinationfile = strcat(outputfolder,mfilename,'_',tstamp1,'.m');
            copyfile(currentfile,destinationfile);   
        end
    end

    if frame>1
        allpdist2Velocities{frame} = iVel; % BSA
    else
        allpdist2Velocities{frame} = nan(length(iBdy),1); % BSA 
    end
    
    allPerims(frame)    = perimeter;
    allCurvatures{frame}  = iCurv; 
    allpdist2Dcumuls{frame} = iDcumul;
    allBoundaries{frame}  = iBdy; 
    allIntensities(frame,:) = iFluor; % row = frame, cols = nfluor
    allCentroids(frame,:) = centroid;  
    allEccents(frame) = eccent;
    allMajAxL(frame) = mjaxl;
    allAspRatio_ml(frame) = aspratio_ml;
    if do_aspratio_ind
        allAspRatio_ind(frame) = aspratio_ind;
        allA0{frame}     = A0;
        allA0c{frame}     = A0c; 
    end
    allLengths(frame)   = length(iBdy);
    allAreas(frame)     = area;
    allBW(:,:,frame)      = BW;
    alliBW(frame,:)       = iBW;
    alloBW(frame,:)       = oBW;
    allinner(frame,:)     = inner;
    allouter(frame,:)     = outer;


   
    if frame>1
        deltaArea             = 100*(allAreas(frame,:)/allAreas(frame-1,:)-1);
        deltaVelocity         = norm(allCentroids(frame,:)-allCentroids(frame-1,:));
        fprintf('Change area: %6.1f%%, Centroid displacement %6.2f pixels\n',deltaArea,deltaVelocity);
    else
        fprintf('\n');
    end
end % end of nframes


%% make kymographs and save data
%----------------------------------------------------
dataRaw = v2struct; % pack all variables into single structure array 
if makekymo

    kf = kymofixed; % e.g. [2 , 0] where 2 means fixedtheta, and 0 means varying 

    for i = 1:length(kf)
        
        binData = cell(1,4);
        if i==1
        [allkymos, alignall, binData1,binData2] = makekymos(dataRaw,kf(i));
        else
        [allkymos, alignall] = makekymos(dataRaw,kf(i));
        end
        tmpKymo.alignangles = binData1(:,1);
        tmpKymo.alignDtheta = binData1(:,2);
        tmpKymo.alignbinidx = binData1(:,3);
        tmpKymo.edges = binData2;
        
        % these are the kymos with aligned quantities centered for display
        tmpKymo.kymovel = allkymos{1};
        tmpKymo.kymocurv = allkymos{2};
        tmpKymo.kymodcumul = allkymos{3};
        tmpKymo.kymofluor = allkymos(4:end); % mch-memb, yfp-actin, cfp-actA
        
        % this is never binned
        tmpKymo.alignbdy = alignall(:,1); 
    
        % these are the aligned quantities just as lists
        % these may be spatially binned if using kymofixedtheta or kymofixedpts 
        tmpKymo.alignvel = alignall(:,2);
        tmpKymo.aligncurv = alignall(:,3);
        tmpKymo.aligndcumul = alignall(:,4);
        tmpKymo.alignfluor = alignall(:,5:end);
    
  
        

        if kf(i) > 0   % kymofixedpts or fixedtheta, vs i == 1
            % save kymooutput
            dataKymoBin = tmpKymo;
            % plot movies and save only based on the first kymograph set (fixedpts or fixedtheta usually)
            plottseries(dataParams,dataRaw,makekymo,tmpKymo,kf(i));
            % plot figures and output them
            pause(0.1);
            figsout0 = plotmovement(dataParams,dataRaw);
            count = 10;
            figsoutFix(1) = plotphyskymos(dataParams,tmpKymo,kf(i),count);
            figsoutFix(2) = plotbiochemkymos(dataParams,tmpKymo,kf(i),count);
            

            %figsout0 = v2struct(f3);
            %figsoutFix = v2struct(f12,f13);
        else % kymovarying
            % save kymooutput
            dataKymoVar = tmpKymo;
            % plot figures and output them
            count = 20;
            figsoutVar(1) = plotphyskymos(dataParams,tmpKymo,kf(i),count);
            figsoutVar(2) = plotbiochemkymos(dataParams,tmpKymo,kf(i),count);
            %figsoutVar = v2struct(f22,f23);
        end
        
    end

    if savefigs == 2
        save(strcat(outputfolder,targname,'_',tstamp1,'_data'),'dataRaw','dataParams',...
            'dataKymoBin','dataKymoVar')
        savefig(figsout0,strcat(outputfolder,targname,'_',tstamp1,'_figsout0.fig'),'compact');
        savefig(figsoutFix,strcat(outputfolder,targname,'_',tstamp1,'_figsoutFix.fig'),'compact');
        savefig(figsoutVar,strcat(outputfolder,targname,'_',tstamp1,'_figsoutVar.fig'),'compact');
    elseif savefigs == 1
        save(strcat(outputfolder,targname,'_',tstamp1,'_data'),'dataParams',...
            'dataKymoBin','dataKymoVar');
        savefig(figsout0,strcat(outputfolder,targname,'_',tstamp1,'_figsout0.fig'),'compact');
        savefig(figsoutFix,strcat(outputfolder,targname,'_',tstamp1,'_figsoutFix.fig'),'compact');
        savefig(figsoutVar,strcat(outputfolder,targname,'_',tstamp1,'_figsoutVar.fig'),'compact');
    end

else
    % plot tseries that are not based on kymos (like bounding ellipse and alignment and masks)
    plottseries(dataParams,dataRaw,makekymo);
    figsout0 = plotmovement(dataParams,dataRaw);
    if savefigs == 2
        save(strcat(outputfolder,targname,'_',tstamp1,'_data'),'dataRaw','dataParams');
        savefig(figsout0,strcat(outputfolder,targname,'_',tstamp1,'_figsout0.fig'),'compact');
    elseif savefigs == 1
        save(strcat(outputfolder,targname,'_',tstamp1,'_data'),'dataParams'); 
        savefig(figsout0,strcat(outputfolder,targname,'_',tstamp1,'_figsout0.fig'),'compact');
    end
end

end % end of runall

%% ------------ auxiliary functions

function plottseries(dataParams,dataRaw,makekymo,varargin)

if ~isempty(varargin)
    tmpKymo = varargin{1};
    kf = varargin{2};
end
cols   = dataParams.cols;
rows   = dataParams.rows;
nframes = dataParams.nframes;
listFrames = dataParams.listFrames;
fnames  = dataParams.fnames; % name of channels
membch = dataParams.membch;
fluorch = dataParams.fluorch;
targname = dataParams.targname;
savgolay_smooth = dataParams.savgolay_smooth;
computecurv = dataParams.computecurv;
im0    = dataParams.im0;
imL    = dataParams.imL;
do_plotEllipseries = dataParams.do_plotEllipseries;
do_plotTimeseries = dataParams.do_plotTimeseries;
do_plotAlignseries = dataParams.do_plotAlignseries;
do_plotMaskseries = dataParams.do_plotMaskseries;
do_aspratio_ind = dataParams.do_aspratio_ind;
savefigs = dataParams.savefigs;
tstamp1 = dataParams.tstamp1;
outputfolder = dataParams.outputfolder;
sat = dataParams.sat;
satlow = dataParams.satlow;
deltaf = dataParams.deltaf;
nfluor = dataParams.nfluor;
tickstep = dataParams.tickstep;
minfont = dataParams.minfont;
basystem = dataParams.basystem;

if makekymo
% these are the kymos
kymovel  = tmpKymo.kymovel; % aligned velocities
kymodcumul = tmpKymo.kymodcumul;
kymocurv  = tmpKymo.kymocurv; % aligned curvatures
kymofluor  = tmpKymo.kymofluor; % aligned kymos for fluorescent channels

% aligned data, these can be spatially binned, but not centered for display like kymos
alignbdy = tmpKymo.alignbdy; % this will never be binned - longer than alignvel/curv/fluor
alignvel = tmpKymo.alignvel; 
aligndcumul = tmpKymo.aligndcumul;
aligncurv = tmpKymo.aligncurv; 
alignfluor = tmpKymo.alignfluor; 
% % full aligned data, with  varying # points (m1*perim) or fixed long perimeter (max m1*perim)
% longbdy = dataKymos.longbdy;
% longvel = dataKymos.longvel;
% longcurv = dataKymos.longcurv;
% longdcumul = dataKymos.longdcumul;
% longfluor = dataKymos.longfluor;
end

cent   = dataRaw.allCentroids;
areas  = dataRaw.allAreas;
allImages  = dataRaw.allImages;
allIrgb = dataRaw.allIrgb;
alloBW = dataRaw.alloBW;
alliBW = dataRaw.alliBW;
allinner = dataRaw.allinner;
allouter = dataRaw.allouter;
if do_aspratio_ind % these are aspect ratio quantities
allA0 = dataRaw.allA0;
allA0c = dataRaw.allA0c;
end
allBoundaries = dataRaw.allBoundaries;

%vflag = any(~isnan(kymovel(:))); % check if velocity is available

%% plot min bounding ellipse

if do_plotEllipseries 
    if ~do_aspratio_ind % must have do_aspratio_ind OFF to compute aspect ratio quantities
        disp('No bounding ellipse available: turn do_aspratio_ind ON for independent aspect ratios with A0 and A0c')
    else
        pause(0.1)
    %    fname = strcat(outputfolder,targname,'_',tstamp1,'_bounding_ellipse');
    %     if strcmp('linux',basystem)
    %         v = VideoWriter(fname,'Motion JPEG AVI');
    %     else
    %         v = VideoWriter(fname,'MPEG-4');
    %     end
    %     
    %     v.FrameRate = 5; %playback framerate
    %     open(v);
    
        f34 = figure(34);
        tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        f34.Position(3:4) = [300, 300];
        for j = 1:nframes
            movieFrame=listFrames(j); %im0:imL;  % TIFF indices
            tmpIrgb = allIrgb{j,membch};
            bdyraw = allBoundaries{j};
            tmpA0 = allA0{j};
            tmpA0c = allA0c{j};
    
            byraw = bdyraw(:,1);
            bxraw = bdyraw(:,2);
            xpts = [bxraw,byraw];
    

            nexttile(1,[2,2]);
            imshow(tmpIrgb); axis equal; hold on;
            plot(xpts(:,1),xpts(:,2),'r','LineWidth',3); 
            FontSize = minfont; %min(minfont,floor(min(cols,rows)/10));
            text(6,6,num2str(movieFrame),'Color','white','FontSize',FontSize)
            title('Bounding ellipse used for eccentricity','FontSize',FontSize)
        
            %set(gca,'YDir','reverse');
              
            Ellipse_plot(tmpA0,tmpA0c);
    
            hold off;
            drawnow
            
            % save as tiff series
            if savefigs > 0
                %pause(0.05)
                F34 = getframe(f34);
                if movieFrame == im0
                    imwrite(F34.cdata,strcat(outputfolder,targname,'_',tstamp1,'_bounding_ellipse.tif'));
                else
                    imwrite(F34.cdata,strcat(outputfolder,targname,'_',tstamp1,'_bounding_ellipse.tif'),...
                    'WriteMode','append');
                end
            end
                
            %writeVideo(v,F4.cdata);
    
        end
        %close(v);
        %close (f34);
    end
end



    %% plot boundary alignment
if do_plotAlignseries 
    for j = 1:nframes         % indices for arrays
        movieFrame=listFrames(j); %im0:imL;  % TIFF indices
        f60 = figure(60); 
        tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        f60.Position(3:4) = [400, 400];
        nexttile(1,[2,2]);
        tmpbdy = alignbdy{j};
        tmpI1rgb = allIrgb{j,membch};
        x0 = cent(j,1); y0 = cent(j,2);
        imshow(tmpI1rgb); axis equal; 
        hold on; 
        ax60 = plot(x0,y0,'m*','MarkerSize',12,'LineWidth',3);
        ax60.Parent.Title.String = strcat(targname,{' '},'tracking ref. angle and center)');
        FontSize = minfont; %min(minfont,floor(min(cols,rows)/10));
        text(6,6,num2str(movieFrame),'Color','white','FontSize',FontSize)
        %plot(tmpbdy(:,2),tmpbdy(:,1),'ro','MarkerSize',10); 
        plot(tmpbdy(:,2),tmpbdy(:,1),'r-','LineWidth',3); 


        midpoint = round(length(tmpbdy)/2);
        plot(tmpbdy(midpoint,2),tmpbdy(midpoint,1),'c*','MarkerSize',12,'LineWidth',3);  

        
        hold off;
        drawnow


        % save as tiff series
        if savefigs > 0
            %pause(0.05)
            F60 = getframe(f60);
            if movieFrame == im0
                imwrite(F60.cdata,strcat(outputfolder,targname,'_',tstamp1,'_tseries_align.tif'));
            else
                imwrite(F60.cdata,strcat(outputfolder,targname,'_',tstamp1,'_tseries_align.tif'),...
                'WriteMode','append');
            end
        end
    end
    %close(f60);
end % end of plot alignseries

    %% timeseries tiff main boundary measurements
if do_plotTimeseries && makekymo
    pause(1);
    
    % find caxis lims for fluor quants
    fluorcax = cell(1,nfluor);
    for k = 1:nfluor
        if kf > 0 %kymofixedpts || kymofixedtheta
            a3 = cat(2,alignfluor{:,k}); % these are same dimensions so can concatenize
        else
            a3 = cat(2,kymofluor{k}); % can only concat if same dimensions, so use kymo which has padding
        end
        fluormax = sat*max(a3(:));
        fluormin = satlow*min(a3(:));
        fluorcax{k} = [fluormin, fluormax];
    end
    % find caxis lims for bdy quants
    avel = cat(2,alignvel{:});
    adc = cat(2,aligndcumul{:});
    acurv = cat(2,aligncurv{:});  

    velmax = max(abs(avel(:)));
    dcmax = max(abs(adc(:)));
    curvmax = max(abs(acurv(:)));
    velcax = [-velmax,velmax];
    dccax = [-dcmax,dcmax];
    curvcax = [-curvmax,curvmax];


    %biophysical quantitites panel
    bio_panel = {'membrane-marker';'actin';'ActA'};

    phys_panel = {'velocity';'cumulative-disp.';'curvature'};

    % fields: name of panel, signal values, color axes, background image
    % each row is a different figure 
    multifig = struct('names',{bio_panel(:),phys_panel(:)},...
        'vals',{alignfluor,[alignvel,aligndcumul,aligncurv]},...
        'caxes',{fluorcax(:),{velcax;dccax;curvcax}},...
        'images',{allIrgb,[allIrgb(:,membch),allIrgb(:,membch),allIrgb(:,membch)]});

    
    % now plot multi-panel time series together
    for k = 1:length(multifig)       



        for j = 1:nframes         % indices for arrays
            movieFrame=listFrames(j); %im0:imL;  % TIFF indices
            %tmpI1rgb = allIrgb{j,:};
            tmpI1pan = multifig(k); % kth element of structure array is kth figure
            %tmpbdy = alignbdy{j};
            fignum = k+70; % figure numbers must not overlap
            currFrame = j; % current Frame

            % tmpvel = alignvel{j};
            % tmpdc = aligndcumul{j};
            % tmpcurv = aligncurv{j};
            % 
            % tmpfluor = alignfluor{j,k};
            % tmpI2rgb = allIrgb{j,k};


            f2vars = v2struct(savgolay_smooth,computecurv,minfont,membch);
            f2 = displaymultipanel(fignum,f2vars,alignbdy,tmpI1pan,currFrame,movieFrame);

       
            if savefigs > 0 % create tiff series
                %pause(0.05)
                F2 = getframe(f2);
    
                if movieFrame == im0
                    imwrite(F2.cdata,strcat(outputfolder,targname,'_',tstamp1,'_tseries_main_',cat(2,multifig(k).names{:}),'.tif'));
                else
                    imwrite(F2.cdata,strcat(outputfolder,targname,'_',tstamp1,'_tseries_main_',cat(2,multifig(k).names{:}),'.tif'),...
                        'WriteMode','append');
                end
            end

        end % end of nframes
        %close(f2)

    end % end of nfluor
end %end of plot main time series

%% plot mask time series

if do_plotMaskseries

    for k = 1:nfluor % include membrane channel here
        f8 = figure(8); clf(f8); 
        f8.Position(3:4) = [400, 400];
        tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

        for j = 1:nframes 
            movieFrame = listFrames(j);    
          
            nexttile(1,[2,2]);
            Iflrgb = allIrgb{j,k};
            iBW = alliBW{j,k};
            oBW = alloBW{j,k};
            inner = allinner{j,k};
            outer = allouter{j,k};
            if makekymo
                tmpbdy = alignbdy{j};
            else
                tmpbdy = allBoundaries{j};
            end
           
            f8ax = imshow(Iflrgb);     axis equal; hold on; 
            green = cat(3, zeros(rows,cols), ones(rows,cols), zeros(rows,cols)); 
            h8 = imshow(green); 
            % Use mask as the AlphaData for the solid green image. 
            set(h8, 'AlphaData', 0.5*oBW) ; 
            red = cat(3, ones(rows,cols), zeros(rows,cols), zeros(rows,cols)); 
            h9 = imshow(red); 
            % Use mask as the AlphaData for the solid green image. 
            set(h9, 'AlphaData', 0.5*iBW) ; 
            % plot boundary points
            h10=plot(inner(:,1),inner(:,2),'r--'); 
            plot(outer(:,1),outer(:,2),'c--');
            plot(tmpbdy(:,2),tmpbdy(:,1),'w--');
            hold off; 
            FontSize = minfont; %min(minfont,floor(min(cols,rows)/10));
            text(6,6,num2str(movieFrame),'Color','white','FontSize',FontSize)
            h10.Parent.YLim=[1 rows];
            h10.Parent.XLim=[1 cols];
            title(['Normalized',' ',fnames{k}(5:end)],'Interpreter','none')
            
    %         f8.Position(3:4) = [400, 400];
            drawnow
            
            
            if savefigs>0
                %pause(0.05)
                F8 = getframe(f8); 
                if movieFrame == im0
                    imwrite(F8.cdata,strcat(outputfolder,targname,'_',tstamp1,'_tseries_mask_',fnames{k}(5:end),'.tif'));
                else
                    imwrite(F8.cdata,strcat(outputfolder,targname,'_',tstamp1,'_tseries_mask_',fnames{k}(5:end),'.tif'),...
                        'WriteMode','append');
                end
            end



        end

    end

end %if do_plotmaskseries
   
end 

function figsout = plotphyskymos(dataParams,tmpKymo,kf,count)

cols   = dataParams.cols;
rows   = dataParams.rows;
nframes = dataParams.nframes;
listFrames = dataParams.listFrames;
fnames  = dataParams.fnames; % name of channels
membch = dataParams.membch;
fluorch = dataParams.fluorch;
targname = dataParams.targname;
savgolay_smooth = dataParams.savgolay_smooth;
computecurv = dataParams.computecurv;
im0    = dataParams.im0;
imL    = dataParams.imL;
savefigs = dataParams.savefigs;
tstamp1 = dataParams.tstamp1;
outputfolder = dataParams.outputfolder;
t_rapa = dataParams.t_rapa;
tickstep = dataParams.tickstep;
dt = dataParams.dt;
deltaf = dataParams.deltaf;
nbins = dataParams.nbins;

kymovel  = tmpKymo.kymovel; % aligned velocities
kymocurv  = tmpKymo.kymocurv; % aligned curvatures
kymodcumul = tmpKymo.kymodcumul; % aligned cumulative displacements
kymofluor  = tmpKymo.kymofluor; % aligned kymos for fluorescent channels
alignbdy = tmpKymo.alignbdy; %allBoundaries;
alignvel = tmpKymo.alignvel; %allVelocities;
aligncurv = tmpKymo.aligncurv; 
alignfluor = tmpKymo.alignfluor; %allIntensities;
edges = tmpKymo.edges;
% 
% cent   = dataRaw.allCentroids;
% areas  = dataRaw.allAreas;
% allImages  = dataRaw.allImages;
% allIrgb = dataRaw.allIrgb;

vflag = any(~isnan(kymovel(:))); % check if velocity is available

tfs = 17; %title font size
afs = 14; % axes labels font size
xlinec = [0.1 0.1 0.1];

 %% ---------------- plot kymographs
    f12 = figure(2+count); 
    t2 = tiledlayout('flow');
    %t2=tiledlayout(3,3);
    t2.Padding     = 'compact';
    t2.TileSpacing = 'compact';

    nexttile(1,[1 3])
    if vflag
    ax4 = imagesc(kymovel,'AlphaData',~isnan(kymovel));
    cmax          = max(abs(kymovel(:)));
    colormap(ax4.Parent,'parula');
    ax4.Parent.CLim = [-cmax,cmax]; 
    colorbar;
%     xticks(1:tickstep:size(kymovel,2)); xticklabels({im0:tickstep:imL});
    xticks(tickstep:tickstep:size(kymovel,2)+1); xticklabels({im0-1+tickstep:tickstep:imL});

%   recall:  
%   edges = - pi : 2*pi/nbins : pi ; 
%   edges = 1 : (n-1)/nbins : n;

    hold on
    xline(t_rapa+0.5,'LineWidth',3,'Color',xlinec,'LineStyle','--');
    if targname(1) == 'L'
        if kf == 2 %fixedtheta
            plot(im0-0.5,0,'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        elseif kf == 1 %fixedpoints
            plot(im0-0.5,round(nbins/2),'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        end
    end
    hold off

    if kf == 2 %kymofixedtheta
        edgeticks = 0.5:nbins/6:nbins+0.5;
        % edgelabs =  pi : -2*pi/6 : -pi;
        edgelabs = {'\pi';'2\pi/3';'\pi/3';'0';'-\pi/3';'-2\pi/3';'-\pi'};
        yticks(edgeticks); yticklabels(edgelabs);
%     elseif kymofixedpts
%         edgeticks = 0.5:nbins/6:nbins+0.5;
%         edgelabs = 0:nbins/6:nbins;
%         yticks(edgeticks); yticklabels({edgelabs});
    end



    ax4.Parent.XLabel.String = 'Frame Number';
    %ax4.Parent.XLabel.FontSize = afs;
    if kf == 2 %kymofixedtheta
    ax4.Parent.YLabel.String = 'Angular position (binned)';
    elseif kf == 1 %kymofixedpts
    ax4.Parent.YLabel.String = 'Perimeter (binned points)';
    else
    ax4.Parent.YLabel.String = 'Perimeter (points)';
    end

    ax4.Parent.YLabel.FontName = 'Arial';
    %ax4.Parent.YLabel.FontSize = afs;
    shading  flat; 
    end

    ax4.Parent.FontSize = afs;
    title('Velocity (pixels/frame) along perimeter','FontSize',tfs);



    %------------------------------------
    nexttile(4,[1 3])

    ax5           = imagesc(kymodcumul,'AlphaData',~isnan(kymodcumul));
    shading flat
%     cm            = zeros(256,3,'uint8');
%     cm(1:256,2)   = linspace(0,255,256);
%     cm(1:256,3)   = linspace(0,255,256);
    cmax = max(abs(kymodcumul(:)));
    %cmin = min(kymo2(kymo2>0));
    colormap(ax5.Parent,'parula');
    %ax5.Parent.Colormap = cm;
    ax5.Parent.CLim = [-cmax,cmax];
    colorbar;
%     xticks(1:tickstep:size(kymofluor,2)); xticklabels({im0:tickstep:imL});
    xticks(tickstep:tickstep:size(kymodcumul,2)); xticklabels({im0-1+tickstep:tickstep:imL});
    
    hold on
    xline(t_rapa+0.5,'LineWidth',3,'Color',xlinec,'LineStyle','--');
    if targname(1) == 'L'
        if kf == 2 %fixedtheta
            plot(im0-0.5,0,'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        elseif kf == 1 %fixedpoints
            plot(im0-0.5,round(nbins/2),'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        end
    end
    hold off
    
    if kf == 2 %kymofixedtheta
        edgeticks = 0.5:nbins/6:nbins+0.5;
        % edgelabs =  pi : -2*pi/6 : -pi;
        edgelabs = {'\pi';'2\pi/3';'\pi/3';'0';'-\pi/3';'-2\pi/3';'-\pi'};
        yticks(edgeticks); yticklabels(edgelabs);
    end
    ax5.Parent.XLabel.String = 'Frame Number';
    ax5.Parent.XLabel.FontName = 'Arial';
    %ax5.Parent.XLabel.FontSize = afs;
    if kf ==2 %kymofixedtheta
    ax5.Parent.YLabel.String = 'Angular position (binned)';
    elseif kf == 1 %kymofixedpts
    ax5.Parent.YLabel.String = 'Perimeter (binned points)';
    else
    ax5.Parent.YLabel.String = 'Perimeter (points)';
    end
    %ax5.Parent.YLabel.FontName = 'Arial';
    ax5.Parent.YLabel.FontSize = afs;

    ax5.Parent.FontSize = afs;
    title('Cumulative displacement (pixels) along perimeter','FontSize',tfs);




    
    %----------------------------------------------
        nexttile(7,[1 3])



    ax6           = imagesc(kymocurv,'AlphaData',~isnan(kymocurv));
    shading flat
%     cmax = max(abs(kymocurv(:)));
%     colormap(ax6.Parent,'parula');
%     ax6.Parent.CLim = [-cmax,cmax];
%     colorbar;

% BINARIZATION 
%     clrs = [0.2422 0.1504 0.6603; 0.9290 0.6940 0.1250];
%     ax6.Parent.CLim = [0,1];
%     colormap(ax6.Parent,clrs);
%     colorbar;

% ORIGINAL COLOR SCHEME
    cmax = max(abs(kymocurv(:)));
    colormap(ax6.Parent,'parula');
    ax6.Parent.CLim = [-cmax,cmax];
    colorbar;

%     xticks(1:tickstep:size(kymocurv,2)); xticklabels({im0:tickstep:imL});
    xticks(tickstep:tickstep:size(kymocurv,2)); xticklabels({im0-1+tickstep:tickstep:imL});
    
    hold on
    xline(t_rapa+0.5,'LineWidth',3,'Color',xlinec,'LineStyle','--');
    if targname(1) == 'L'
        if kf == 2 %fixedtheta
            plot(im0-0.5,0,'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        elseif kf == 1 %fixedpoints
            plot(im0-0.5,round(nbins/2),'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        end
    end
    hold off

    if kf == 2 %kymofixedtheta
        edgeticks = 0.5:nbins/6:nbins+0.5;
        %edgelabs =  pi : -2*pi/6 : -pi;
        edgelabs = {'\pi';'2\pi/3';'\pi/3';'0';'-\pi/3';'-2\pi/3';'-\pi'};
        yticks(edgeticks); yticklabels(edgelabs);
    end
    ax6.Parent.XLabel.String = 'Frame Number';
    ax6.Parent.XLabel.FontName = 'Arial';
    %ax6.Parent.XLabel.FontSize = afs;
    if kf == 2 %kymofixedtheta
    ax6.Parent.YLabel.String = 'Angular position (binned)';
    elseif kf == 1 %kymofixedpts
    ax6.Parent.YLabel.String = 'Perimeter (binned points)';
    else
    ax6.Parent.YLabel.String = 'Perimeter (points)';
    end
    ax6.Parent.YLabel.FontName = 'Arial';
    %ax6.Parent.YLabel.FontSize = afs;

    ax6.Parent.FontSize = afs;
    if computecurv == 1
        title('Curvature (circumcenter) along perimeter','FontSize',tfs);
    else
       title('Curvature (fitcircle) along perimeter','FontSize',tfs);
    end
    
    %


    title(t2,[targname,': physical kymographs']);

        % adjust the figure positions 
    %[left bottom width height]
    f12.Position(3:4) = [700 700];
    drawnow;

%     %% adjust font size in all figures
% fh = findall(0,'Type','Figure');
% set( findall(fh, '-property', 'fontsize'), 'fontsize', 14)

    %figsout = v2struct(f12);
    figsout = f12;
    if savefigs>0 && kf>0
    print(f12,strcat(outputfolder,targname,'_',tstamp1,'_fig_physkymos_Fix'),'-dpng')
    elseif savefigs>0 && kf==0
    print(f12,strcat(outputfolder,targname,'_',tstamp1,'_fig_physkymos_Var'),'-dpng')

    end


end

function figsout = plotbiochemkymos(dataParams,tmpKymo,kf,count)

cols   = dataParams.cols;
rows   = dataParams.rows;
nframes = dataParams.nframes;
listFrames = dataParams.listFrames;
fnames  = dataParams.fnames; % name of channels
membch = dataParams.membch;
fluorch = dataParams.fluorch;
targname = dataParams.targname;
savgolay_smooth = dataParams.savgolay_smooth;
computecurv = dataParams.computecurv;
im0    = dataParams.im0;
imL    = dataParams.imL;
savefigs = dataParams.savefigs;
tstamp1 = dataParams.tstamp1;
outputfolder = dataParams.outputfolder;
t_rapa = dataParams.t_rapa;
tickstep = dataParams.tickstep;
dt = dataParams.dt;
deltaf = dataParams.deltaf;
nbins = dataParams.nbins;

sat = dataParams.sat;
satlow = dataParams.satlow;

kymovel  = tmpKymo.kymovel; % aligned velocities
kymocurv  = tmpKymo.kymocurv; % aligned curvatures
kymodcumul = tmpKymo.kymodcumul; % aligned cumulative displacements
kymofluor  = tmpKymo.kymofluor; % aligned kymos for fluorescent channels
alignbdy = tmpKymo.alignbdy; %allBoundaries;
alignvel = tmpKymo.alignvel; %allVelocities;
aligncurv = tmpKymo.aligncurv; 
alignfluor = tmpKymo.alignfluor; %allIntensities;
edges = tmpKymo.edges;

vflag = any(~isnan(kymovel(:))); % check if velocity is available

tfs = 17; % title font size
afs = 14; % axes font size
xlinec = [0.1 0.1 0.1];

 %% ---------------- plot kymographs
    f13 = figure(3+count); 
    t2 = tiledlayout('flow');

    %t2=tiledlayout(3,3);
    t2.Padding     = 'compact';
    t2.TileSpacing = 'compact';

    nexttile(1,[1 3])
    if vflag
    kymfl = kymofluor{fluorch(1)};
    ax4 = imagesc(kymfl,'AlphaData',~isnan(kymfl));
    cmax = max(kymfl(:));
    cmin = min(kymfl(:));
    colormap(ax4.Parent,'parula');

    ax4.Parent.CLim = [satlow*cmin,sat*cmax]; 

    colorbar;
%     xticks(1:tickstep:size(kymovel,2)); xticklabels({im0:tickstep:imL});
    xticks(tickstep:tickstep:size(kymfl,2)+1); xticklabels({im0-1+tickstep:tickstep:imL});
    
    hold on
    xline(t_rapa+0.5,'LineWidth',3,'Color',xlinec,'LineStyle','--');
    if targname(1) == 'L'
        if kf == 2 %fixedtheta
            plot(im0-0.5,0,'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        elseif kf == 1 %fixedpoints
            plot(im0-0.5,round(nbins/2),'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        end
    end
    hold off
    
    if kf == 2 %kymofixedtheta
        edgeticks = 0.5:nbins/6:nbins+0.5;
        % edgelabs =  pi : -2*pi/6 : -pi;
        edgelabs = {'\pi';'2\pi/3';'\pi/3';'0';'-\pi/3';'-2\pi/3';'-\pi'};
        yticks(edgeticks); yticklabels(edgelabs);
    end
    ax4.Parent.XLabel.String = 'Frame Number';
    %ax4.Parent.XLabel.FontSize = afs;
    if kf == 2 %kymofixedtheta
    ax4.Parent.YLabel.String = 'Angular position (binned)';
    elseif kf == 1 % kymofixedpts
    ax4.Parent.YLabel.String = 'Perimeter (binned points)';
    else
    ax4.Parent.YLabel.String = 'Perimeter (points)';
    end
    ax4.Parent.YLabel.FontName = 'Arial';
    %ax4.Parent.YLabel.FontSize = afs;
    shading  flat; 
    end
    ax4.Parent.FontSize = afs;
    title('Membrane marker intensity along perimeter','FontSize',tfs);




    % ------------------------
    nexttile(4,[1 3])
    kymfl = kymofluor{fluorch(3)};
    ax6           = imagesc(kymfl,'AlphaData',~isnan(kymfl));
    shading flat
    cmax = max(kymfl(:));
    cmin = min(kymfl(:));


    colormap(ax6.Parent,'parula');
    ax6.Parent.CLim = [satlow*cmin,sat*cmax];

    colorbar;
%     xticks(1:tickstep:size(kymocurv,2)); xticklabels({im0:tickstep:imL});
    xticks(tickstep:tickstep:size(kymfl,2)); xticklabels({im0-1+tickstep:tickstep:imL});
    
    hold on
    xline(t_rapa+0.5,'LineWidth',3,'Color',xlinec,'LineStyle','--');
    if targname(1) == 'L'
        if kf == 2 %fixedtheta
            plot(im0-0.5,0,'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        elseif kf == 1 %fixedpoints
            plot(im0-0.5,round(nbins/2),'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        end
    end
    hold off
    
    if kf == 2 %kymofixedtheta
        edgeticks = 0.5:nbins/6:nbins+0.5;
        % edgelabs =  pi : -2*pi/6 : -pi;
        edgelabs = {'\pi';'2\pi/3';'\pi/3';'0';'-\pi/3';'-2\pi/3';'-\pi'};
        yticks(edgeticks); yticklabels(edgelabs);
    end
    ax6.Parent.XLabel.String = 'Frame Number';
    ax6.Parent.XLabel.FontName = 'Arial';
    %ax6.Parent.XLabel.FontSize = afs;
    if kf == 2 %kymofixedtheta
        ax6.Parent.YLabel.String = 'Angular position (binned)';
    elseif kf == 1 %kymofixedpts
        ax6.Parent.YLabel.String = 'Perimeter (binned points)';
    else
        ax6.Parent.YLabel.String = 'Perimeter (points)';
    end
    ax6.Parent.YLabel.FontName = 'Arial';
    %ax6.Parent.YLabel.FontSize = afs;
    ax6.Parent.FontSize = afs;
    title('ActA intensity along perimeter','FontSize',tfs);



    %
    %------------------------------
    nexttile(7,[1 3])


    kymfl = kymofluor{fluorch(2)};
    ax5           = imagesc(kymfl,'AlphaData',~isnan(kymfl));
    shading flat

%     cmax = max(kymfl(:));
%     cmin = min(kymfl(:));
%     colormap(ax5.Parent,'parula');
%     ax5.Parent.CLim = [satlow*cmin,sat*cmax];
%     colorbar;

    % BINARIZATION
%     clrs = [0.2422 0.1504 0.6603; 0.9290 0.6940 0.1250];
%     ax5.Parent.CLim = [0,1];
%     colormap(ax5.Parent,clrs);
%     colorbar;

    % ORIGINAL COLOR SCHEME
    cmax = max(kymfl(:));
    cmin = min(kymfl(:));
    colormap(ax5.Parent,'parula');
    ax5.Parent.CLim = [satlow*cmin,sat*cmax];
    colorbar;



%     xticks(1:tickstep:size(kymofluor,2)); xticklabels({im0:tickstep:imL});
    xticks(tickstep:tickstep:size(kymfl,2)); xticklabels({im0-1+tickstep:tickstep:imL});
    
    hold on
    xline(t_rapa+0.5,'LineWidth',3,'Color',xlinec,'LineStyle','--');
    if targname(1) == 'L'
        if kf == 2 %fixedtheta
            plot(im0-0.5,0,'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        elseif kf == 1 %fixedpoints
            plot(im0-0.5,round(nbins/2),'r>','MarkerSize',8,'MarkerFaceColor',[1 0.5 0]);
        end
    end
    hold off
    
    if kf == 2 %kymofixedtheta
        edgeticks = 0.5:nbins/6:nbins+0.5;
        %edgelabs =  pi : -2*pi/6 : -pi;
        edgelabs = {'\pi';'2\pi/3';'\pi/3';'0';'-\pi/3';'-2\pi/3';'-\pi'};
        yticks(edgeticks); yticklabels(edgelabs);
    end
    ax5.Parent.XLabel.String = 'Frame Number';
    ax5.Parent.XLabel.FontName = 'Arial';
    %ax5.Parent.XLabel.FontSize = afs;
    if kf == 2 %kymofixedtheta
    ax5.Parent.YLabel.String = 'Angular position (binned)';
    elseif kf == 1 %kymofixedpts
    ax5.Parent.YLabel.String = 'Perimeter (binned points)';
    else
    ax5.Parent.YLabel.String = 'Perimeter (points)';
    end
    ax5.Parent.YLabel.FontName = 'Arial';
    %ax5.Parent.YLabel.FontSize = afs;
    ax5.Parent.FontSize = afs;
    title('Actin intensity along perimeter','FontSize',tfs);


%     
%     title(t2,[targname,': biochemical kymographs (top ',num2str(100*(1-sat)),...
%         '% sat.)']);
    title(t2,[targname,': biochemical kymographs']);
    drawnow

%     %% adjust font size in all figures
% fh = findall(0,'Type','Figure');
% set( findall(fh, '-property', 'fontsize'), 'fontsize', 14)

        % adjust the figure positions 
    %[left bottom width height]
    f13.Position(3:4) = [700 700];

    %figsout = v2struct(f13);
    figsout = f13;
    if savefigs>0 && kf>0
    print(f13,strcat(outputfolder,targname,'_',tstamp1,'_fig_biochemkymos_Fix'),'-dpng')
    elseif savefigs>0 && kf==0
    print(f13,strcat(outputfolder,targname,'_',tstamp1,'_fig_biochemkymos_Var'),'-dpng')
    end



end

function figsout = plotmovement(dataParams,dataRaw)

cols   = dataParams.cols;
rows   = dataParams.rows;
nframes = dataParams.nframes;
listFrames = dataParams.listFrames;
targname = dataParams.targname;
im0    = dataParams.im0;
imL    = dataParams.imL;
savefigs = dataParams.savefigs;
tstamp1 = dataParams.tstamp1;
outputfolder = dataParams.outputfolder;
tickstep = dataParams.tickstep;
t_rapa = dataParams.t_rapa;
do_aspratio_ind = dataParams.do_aspratio_ind;
do_fixedbox = dataParams.do_fixedbox;
do_plotBoundaryOverlay = dataParams.do_plotBoundaryOverlay;

cent   = dataRaw.allCentroids;
areas  = dataRaw.allAreas;
aspratio_ml  = dataRaw.allAspRatio_ml;
aspratio_ind = dataRaw.allAspRatio_ind;
eccent   = dataRaw.allEccents;
bdy = dataRaw.allBoundaries;
allMajAxL = dataRaw.allMajAxL;



afs = 14; % axes font size
tfs = afs + 3; % title font size
%% ------------- plot centroid and velocity
    f3 = figure(3); 
    t1=tiledlayout(4,3); %(5,3) if including centroid velocity
    fdelta = imL-im0+1;
    cm=colormap(parula(fdelta)); % use as many colors as there frames
%     minRow = rows;
%     maxRow = 1;
%     minCol = cols;
%     maxCol = 1;
    minRow = 1;
    maxRow = rows;
    minCol = 1;
    maxCol = cols;
    velocityCentroid = [nan;vecnorm(diff(cent)')'];

    if do_fixedbox
        % FOR FIXED BOX SIZE
        boxt1 = [1,180,1,180]; % xmin, xmax, ymin, ymax
        maxrc = boxt1([2,4]);
        dispadj = round(maxrc/2); % center x, center y
        dxadj = cent(1,:)-fliplr(dispadj); % x, y adjust
    else
        boxt1 = [minCol maxCol minRow maxRow];
        dxadj   = [0 0];
    end

    for k = 1:nframes         % indices for arrays

        movieFrame=listFrames(k); %movieframe, usually same as im0:imL;  % TIFF indices
        %b = mod(a,m) returns the remainder after division of a by m, 
        %thiscm = cm(mod(k-1,255)+1,:); 
        thiscm = cm(movieFrame-im0+1,:);
        thiscm2 = [0 0.4470 0.7410];
        backcm = 0.3*[1 1 1];
        tmpbdy=bdy{k};

        if do_fixedbox
            boxsize = boxt1;
        else
            % adjust overall min, max after each frame
            minCol = min(minCol,min(tmpbdy(:,2)-dxadj(1)));
            maxCol = max(maxCol,max(tmpbdy(:,2)-dxadj(1)));
            minRow = min(minRow,min(tmpbdy(:,1)-dxadj(2)));
            maxRow = max(maxRow,max(tmpbdy(:,1)-dxadj(2)));
            boxsize = [minCol maxCol minRow maxRow];
        end
        
        nexttile(1)
        % plot centroid: already in x, y coords (not row, col coords)
        %ax1=plot(cent(k,1),rows+1-cent(k,2),'.','Color',thiscm,'MarkerSize',30);
        ax1=plot(cent(k,1)-dxadj(1),cent(k,2)-dxadj(2),'.','Color',thiscm,'MarkerSize',30);
        axis equal, hold on
        set(ax1.Parent,'Color',backcm)
        set(ax1.Parent,'Ydir','reverse')
        axis(boxsize)

        nexttile(2)  

        % plot boundary: convert row index to y coord
        %ax2=plot(tmpbdy(:,2),rows+1-tmpbdy(:,1),'-','Color',thiscm,'LineWidth',2);
        ax2=plot(tmpbdy(:,2)-dxadj(1),tmpbdy(:,1)-dxadj(2),'-','Color',thiscm,'LineWidth',2);
        axis equal, hold on
        set(ax2.Parent,'Color',backcm)
        set(ax2.Parent,'Ydir','reverse')        
        axis(boxsize)
        ax2.Parent.CLim = [im0 imL+1]; 


        
    end

    % add colorbar to ax2

    %cc = colorbar('Ticks',im0+0.5:tickstep:imL+0.5,'TickLabels',{im0:tickstep:imL});  
    cc = colorbar('Ticks',im0-0.5+tickstep:tickstep:imL+0.5,'TickLabels',{im0-1+tickstep:tickstep:imL});  
    cc.Label.String = 'Frames';
    cc.Label.FontSize = afs;
    cc.FontSize = afs;

    % next  plots are stair plots with only 1 color, FAST
    % stair eccentricity, major axis length/asp ratio, area, centroid velocity plot with color variation 

    nexttile(4,[1 3])

    % plot eccentricity
    tmpeccent = nan(fdelta,1);
    tmpeccent(listFrames)=eccent;
    ax3a = stairs([tmpeccent;tmpeccent(end)],'Color',thiscm2,'LineWidth',2);
    ax3a.Parent.CLim = [im0 imL+1]; 
    %set(gca,'Color',backcm)
    hold on



    
    nexttile(7,[1 3])

%     % plot major axis length
%     tmpmaxl = nan(fdelta,1);
%     tmpmaxl(listFrames)=allMajAxL;
%     ax3b=stairs([tmpmaxl;tmpmaxl(end)],'.-','Color',thiscm2,'LineWidth',2);
%     ax3b.Parent.CLim = [im0 imL+1]; 
%     hold on

    % plot aspect ratio
    tmpaspr = nan(fdelta,1);
%     if do_aspratio_ind
%         tmpaspr(listFrames)=aspratio_ind;
%     else % use matlab regionprops
%         tmpaspr(listFrames)=aspratio_ml;
%     end
    tmpaspr(listFrames)=aspratio_ml;
    ax3b = stairs([tmpaspr;tmpaspr(end)],'Color',thiscm2,'LineWidth',2);
    ax3b.Parent.CLim = [im0 imL+1]; 
    %set(gca,'Color',backcm)
    hold on

    

    nexttile(10,[1 3])

    % plot areas
    tmpareas = nan(fdelta,1);
    tmpareas(listFrames)=areas;
    ax3l=stairs([tmpareas;tmpareas(end)],'.-','Color',thiscm2,'LineWidth',2);
    ax3l.Parent.CLim = [im0 imL+1]; 
    %set(gca,'Color',backcm)
    hold on
    

%     nexttile(13,[1 3])

%     % stair centroid velocity plot with color variation 
%     tmpcvel = nan(fdelta,1);
%     tmpcvel(listFrames)=velocityCentroid;
%     ax3r=stairs([tmpcvel;tmpcvel(end)],'.-','Color',thiscm2,'LineWidth',2);
%     ax3r.Parent.CLim = [im0 imL+1];
%     %set(gca,'Color',backcm)
%     hold on


    
    if do_fixedbox
        boxsize = boxt1;
    else
        %extra padding when adjusting box for every target
        minCol = max(   1,minCol-5);
        maxCol = min(cols,maxCol+5);
        minRow = max(   1,minRow-5);
        maxRow = min(rows,maxRow+5);
        boxsize = [minCol maxCol minRow maxRow];
    end

    %
    nexttile(1)
    axis(boxsize), axis equal
    hold off

    ax1.Parent.FontName = 'Arial';
    ax1.Parent.FontSize = afs;
    ax1.Parent.XLim = boxsize(1:2);
    ax1.Parent.YLim = boxsize(3:4);
    title('Centroid Position','FontSize',tfs)
    %ax1.Parent.YTickLabel = ax1.Parent.YTickLabel(end:-1:1);
    %
    nexttile(2)
    axis(boxsize), axis equal
    hold off
    ax2.Parent.FontName = 'Arial';
    ax2.Parent.FontSize = afs;
    ax2.Parent.XLim = boxsize(1:2);
    ax2.Parent.YLim = boxsize(3:4);
    title('Boundary','FontSize',tfs)
    %ax2.Parent.YTickLabel = ax2.Parent.YTickLabel(end:-1:1);

    %
    nexttile(4,[1 3])
    xline(t_rapa+1,'LineWidth',2)
    hold off
    ax3a.Parent.YLabel.String = 'Eccentricity';
    ax3a.Parent.YLim=[0 1];
    ax3a.Parent.YLabel.FontName = 'Arial';
    %ax3a.Parent.YLabel.FontSize = afs;
    ax3a.Parent.XLim=[im0 imL+1];
    ax3a.Parent.XTick=im0-0.5+tickstep:tickstep:imL+0.5;
    ax3a.Parent.XTickLabel=im0-1+tickstep:tickstep:imL;
    ax3a.Parent.XLabel.String='Frames';
    %ax3a.Parent.XLabel.FontSize=afs;
    ax3a.Parent.FontSize = afs;
%
    nexttile(7,[1 3])
    xline(t_rapa+1,'LineWidth',2)
    hold off
    ax3b.Parent.YLabel.String = 'Aspect Ratio';
    ax3b.Parent.YLim(1) = 1; % lower limit of aspect ratios
    %ax3b.Parent.YLim=[1 3]; % range of aspect ratios
    %ax3b.Parent.YLabel.String = 'Major Axis Length (pixels)';
    ax3b.Parent.YLabel.FontName = 'Arial';
    %ax3b.Parent.YLabel.FontSize = afs;
    ax3b.Parent.XLim=[im0 imL+1];
    ax3b.Parent.XTick=im0-0.5+tickstep:tickstep:imL+0.5;
    ax3b.Parent.XTickLabel=im0-1+tickstep:tickstep:imL;
    ax3b.Parent.XLabel.String='Frames';
    %ax3b.Parent.XLabel.FontSize=afs;
    ax3b.Parent.FontSize = afs;
    % 
    nexttile(10,[1 3])
    xline(t_rapa+1,'LineWidth',2)
    hold off
    ax3l.Parent.YLabel.String = 'Area (pixels)';
    ax3l.Parent.YLabel.FontName = 'Arial';
    %ax3l.Parent.YLabel.FontSize = afs;
    ax3l.Parent.XLim=[im0 imL+1];
    ax3l.Parent.XTick=im0-0.5+tickstep:tickstep:imL+0.5;
    ax3l.Parent.XTickLabel=im0-1+tickstep:tickstep:imL;
    ax3l.Parent.XLabel.String='Frames';
    %ax3l.Parent.XLabel.FontSize=afs;
    ax3l.Parent.FontSize = afs;
    % 
%     nexttile(13,[1 3])
%     xline(t_rapa+1,'LineWidth',2)
%     hold off
%     ax3r.Parent.YLabel.String = 'Centroid Velocity (pixels/frame)';
%     ax3r.Parent.YLabel.FontName = 'Arial';
%     %ax3r.Parent.YLabel.FontSize = afs;
%     ax3r.Parent.XLim=[im0 imL+1];
%     ax3r.Parent.XTick=im0-0.5+tickstep:tickstep:imL+0.5;
%     ax3r.Parent.XTickLabel=im0-1+tickstep:tickstep:imL;
%     ax3r.Parent.XLabel.String='Frames';
%     %ax3r.Parent.XLabel.FontSize=afs;
%     ax3r.Parent.FontSize = afs;
    %
    t1.Padding     = 'compact';
    t1.TileSpacing = 'compact';
    t1.Title.String = targname;
    t1.Title.FontSize = 12;
    t1.Title.FontName = 'Arial';
    t1.Title.Color = [0 0 0];
    
    % adjust the figure positions 
    %[left bottom width height]
    f3.Position(3:4) = [800, 900];
    drawnow
    
    set(f3, 'InvertHardcopy', 'off')

if do_plotBoundaryOverlay
    f4 = figure(4); hold on
    f4t = tiledlayout('flow','TileSpacing','compact','Padding','compact');
    boxsize = boxt1;
    for k = 1:nframes
        tmpbdy=bdy{k};
        movieFrame=listFrames(k); %movieframe, usually same as im0:imL;  % TIFF indices

        if ~do_fixedbox
            % adjust overall min, max after each frame
            minCol = min(minCol,min(tmpbdy(:,2)-dxadj(1)));
            maxCol = max(maxCol,max(tmpbdy(:,2)-dxadj(1)));
            minRow = min(minRow,min(tmpbdy(:,1)-dxadj(2)));
            maxRow = max(maxRow,max(tmpbdy(:,1)-dxadj(2)));
            boxsize = [minCol maxCol minRow maxRow];
        end
         
        cm4=colormap(parula(fdelta)); % use as many colors as there frames
        thiscm4 = cm4(movieFrame-im0+1,:);

        nexttile(1)
        ax4=plot(tmpbdy(:,2)-dxadj(1),tmpbdy(:,1)-dxadj(2),'-','Color',thiscm4,'LineWidth',6);
        axis equal, hold on
        set(ax4.Parent,'Color',backcm)
        set(ax4.Parent,'Ydir','reverse')        
        axis(boxsize)
        title(strcat('Boundary'),'FontSize',afs+12);
        ax4.Parent.CLim = [im0 imL+1]; 
        %cc = colorbar('Ticks',im0+0.5:tickstep:imL+0.5,'TickLabels',{im0:tickstep:imL});  
        cc = colorbar('Ticks',im0-0.5+tickstep:tickstep:imL+0.5,'TickLabels',{im0-1+tickstep:tickstep:imL});  
        cc.Label.String = 'Frames';

    end


    if ~do_fixedbox
        %extra padding when adjusting box for every target
        minCol = max(   1,minCol-5);
        maxCol = min(cols,maxCol+5);
        minRow = max(   1,minRow-5);
        maxRow = min(rows,maxRow+5);
        boxsize = [minCol maxCol minRow maxRow];
    end
    hold off
    title(f4t,targname,'FontSize',afs+14)
    ax4.Parent.FontName = 'Arial';
    ax4.Parent.FontSize = afs+4;
%     ax4.Parent.XLim = boxsize(1:2);
%     ax4.Parent.YLim = boxsize(3:4);
    axis(boxsize)

    % width height
    f4.Position(3:4) = [500,425];
    drawnow
    set(f4,'InvertHardcopy','off')
    
    % OUTPUTS
    
    figsout(1) = f3;
    figsout(2) = f4;

    if savefigs>0
    print(f3,strcat(outputfolder,targname,'_',tstamp1,'_fig_macrosummary'),'-dpng')
    print(f4,strcat(outputfolder,targname,'_',tstamp1,'_fig_boundary'),'-dpng')
    end

else

    % OUTPUTS
    
    figsout(1) = f3;

    if savefigs>0
    print(f3,strcat(outputfolder,targname,'_',tstamp1,'_fig_macrosummary'),'-dpng')
    end

    
end





end
function figout = displayboundary(I1rgb,bdy,Borigs,Bkpoints)
    % this displays boundary 
    f1= figure(1); clf(f1);
    tf1 = tiledlayout(2,4);
    tf1.Padding     = 'compact';
    tf1.TileSpacing = 'compact';

    nexttile(1,[2 2])
    imshow(I1rgb);
    hold on;
    plot(bdy(1,2),bdy(1,1),'go','LineWidth',10)
    plot(bdy(:,2),bdy(:,1),'g*','LineWidth',2)
    plot(Borigs(1,2),Borigs(1,1),'bo','LineWidth',10)
    plot(Borigs(:,2),Borigs(:,1),'b*','LineWidth',2)
    title('Raw boundary (blue)')
    hold off

    nexttile(3,[2 2])
    imshow(I1rgb);
    hold on;
    plot(bdy(1,2),bdy(1,1),'og','LineWidth',10)
    plot(bdy(:,2),bdy(:,1),'*g','LineWidth',2)
    plot(Bkpoints(1,2),Bkpoints(1,1),'ob','LineWidth',10)
    plot(Bkpoints(:,2),Bkpoints(:,1),'*b','LineWidth',2)
    title('Spaced and smoothed boundary(blue)')
    hold off
   
    sgtitle('Adjustments to raw boundary (green) for circumcenter-based curvature computation')
    drawnow    
    figout = f1;
end

function figout = displaymultipanel(fignum,fvars,I1bdy,I1pan,currFrame,movieFrame)
% recall: f2vars = v2struct(savgolay_smooth,computecurv,movieFrame,tmpbdy,minfont,membch);
v2struct(fvars);
dmark = 45;
FontSize = minfont; %min(minfont,floor(min(cols,rows)/10));
interpBdy = I1bdy{currFrame};
nsubpanel = length(I1pan(1).names); % # of subpanels

f2= figure(fignum); clf(f2);
tf2 = tiledlayout(2,2*nsubpanel);
tf2.Padding     = 'compact';
tf2.TileSpacing = 'compact';
f2.Position(3:4) = [1500, 500];

Ilabels = I1pan.names; % names of boundary quantities
Icax = I1pan.caxes; % color axis for boundary quantity
Irgb = I1pan.images(currFrame,:); % background image for each sub-panel
Ivals = I1pan.vals(currFrame,:); % boundary quantities to plot for each sub-panel



for j = 1:nsubpanel

    nexttile(2*j-1,[2 2])
    imshow(Irgb{j}); 

    hold on;
    f2ax1 = scatter3(interpBdy(:,2),interpBdy(:,1),zeros(length(interpBdy),1)...
        ,[],Ivals{j},'filled','SizeData',dmark);
    title(Ilabels{j})
    colorbar;  colormap(parula);
    f2ax1.Parent.CLim = Icax{j}; %[-0.03,0.03]; %caxis([-0.11 0.11]); 
    text(6,6,num2str(movieFrame),'Color','white','FontSize',FontSize)
    hold off

end


    hold off
 
    drawnow
    % adjust position
    %f2.Position(3:4) = [1500, 500];

    figout = f2;

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


function varargout = aligninitboundary(I1rgb,bdy,z0,thetatop,varargin)
    % !!! make bdy quantities non-circular for alignment
    tmpbdy = bdy(1:end-1,:);
    %figure(71); imshow(I1rgb); hold on; axis equal; plot(tmpbdy(1,2),tmpbdy(1,1),'r*'); hold off
    tmpvarsin = cell(length(varargin),1);
    tmpvarsout = cell(length(varargin),1);
    for i=1:length(varargin)
        tmpvarsin{i} = varargin{i}(1:end-1);
    end
    x0 = z0(1);
    y0 = z0(2);
    theta = atan2(y0-tmpbdy(:,1), tmpbdy(:,2)-x0); % deltay = y0 - bdyrow; delta x = bdycol - x0;
        % theta here will be always between -pi to pi

    thetatop; % this is the angle pi away from thetacenter
        % thetatop will always be between -pi to pi

    %[~,idx]=(min(abs(theta-theta0)));
    [~,idx]=(min(abs(theta-thetatop))); % find index closest to thetapi
    % this will give index of point to align top of kymograph

    % m = 1-idx = negative number --> CW rotation of position 1 label;
    % shifts first m entries to end
    % theta0 marks top row of kymograph
    tmpbdyalign = circshift(tmpbdy, 1-idx); 

    %figure(72); imshow(I1rgb); hold on; axis equal; plot(tmpbdyalign(1,2),tmpbdyalign(1,1),'r*'); hold off

    % make circular at end
    varargout{1} = [tmpbdyalign; tmpbdyalign(1,:)]; 
    for i=1:length(varargin)
        tmpvarsout{i} = circshift(tmpvarsin{i},1-idx);
        % make circular at end
        varargout{i+1} = [tmpvarsout{i}; tmpvarsout{i}(1) ];
    end
end


function varargout = alignboundary(newbdy,oldbdy,varargin)
% newbdy and oldbdy are coordinates (y,x) relative to centroid

% The analysis assumes that the new boundary is longer; if not, switch them
% [alignedbdy,mindst,aligneddata1,...] = alignboundary(newbdy,oldbdy,data1,...)


    % !!! make bdy quantities non-circular for alignment
    tmpvarsin = cell(length(varargin),1);
    tmpvarsout = cell(length(varargin),1);
    for i=1:length(varargin)
        tmpvarsin{i} = varargin{i}(1:end-1);
    end

    sw = logical(length(newbdy)>length(oldbdy));
    % make b1, b2 non-circular
    if sw %new bdy is longer- set this to b1
       b1 = newbdy(1:end-1,:);
       b2 = oldbdy(1:end-1,:);
    else % old bdy is longer- set this to b1
       b1 = oldbdy(1:end-1,:);
       b2 = newbdy(1:end-1,:);
    end
    n1 = length(b1); % number of points in longer bdy 
    n2 = length(b2); % number of points in shorter bdy
    
    d=zeros(n2,1); % based on shorter bdy (b2)
    idx=[];
    % add unique indices (between 1 and n1) to idx
    % until we get enough = n1-n2
    while length(idx)<n1-n2  % difference in number of points
        deltan = n1-n2-length(idx);
        nidx=ceil(n1*rand(deltan,1)); % nidx = deltan random numbers between 1 and n1
        idx=unique(sort([idx;nidx])); % idx = collect unique numbers in idx,nidx in ascending order
    end
    % find indices in 1:n1 but not in idx and sort
    % kidx has n1-(n1-n2) = n2 # of points
    % kidx collects n2 unique random indices between 1:n1, 
    kidx = sort(setdiff(1:n1,idx))';
    % b1s: list of n2 unique random boundary points from b1
    % b1s has same length as b2 (shorter bdy)
    b1s = b1(kidx,:);

    dst = zeros(length(b1s),1);
    
    for i=0:length(b1s)-1
        % find distance b/w x,y coords of b2 and b1s (both same # points)
        % cycle through all possible rotations 
        % b1,b2 are listed CW initially 
        d=b2-circshift(b1s,i);% rotate initial point in b1s CCW i positions
        dst(i+1)=mean(sqrt(sum(d'.^2)));
    end
    %figure(72)
    %plot(dst)
    [mindst,minidx]=min(dst);
    %minidx - 1 = how much b1s initial point should be rotated CCW (+circshift)

    if sw % new bdy longer, then newbdy = b1, so do same CCW direction as b1s
        newbdyalign = circshift(b1,minidx-1);
        varargout{1} = [newbdyalign; newbdyalign(1,:)]; % make circular
        for i=1:length(varargin)
            tmpvarsout{i} = circshift(tmpvarsin{i},minidx-1);
            varargout{i+1} = [tmpvarsout{i}; tmpvarsout{i}(1)];
        end
    else % new bdy shorter, then newbdy = b2, so do opp CW but same amount as b1s 
        % CCW direction for b1s = (minidx-1)
        % CW direction for b2 = - (minidx-1)   
        newbdyalign = circshift(b2,-(minidx-1)); 
        varargout{1} = [newbdyalign; newbdyalign(1,:)]; % make circular
        for i=1:length(varargin)
            tmpvarsout{i}  = circshift(tmpvarsin{i},-(minidx-1));
            varargout{i+1} = [tmpvarsout{i}; tmpvarsout{i}(1)];
        end

    end
end


function varargout = alignboundary_ba(newbdy,oldbdy,varargin)
% inputs:
% newbdy, oldbdy, varargin: circular lists (repetition of first point at end of list)
% varargin: vel, curv, fluor only 

% function determines alignment 
% create equal length boundary vectors (without repetition of first point)
% rotate newbdy and compute sum of displacement map between the two ordered lists (newbdy and oldbdy)
% aligned newbdy = rotation of newbdy that minimizes the sum of displacements 
% computes a motion tracking, with 1:1 correspondence

% outputs: all aligned, circular, and interpolated bac to original newbdy length if newbdy was adjusted
% varargout{1} = newbdy
% varargout{2} = curvature for newbdy
% varargout{3} = velocity (px/frame) computed based on pdist2 (NOT motion tracking)
% varargout{4} = fluor for newbdy


tmpvarsin = cell(1,length(varargin));
tmpvarsinu = cell(1,length(varargin));
tmpvarsout = cell(1,length(varargin));
tmpvarsoutC = cell(1,length(varargin));
varargout = cell(1,length(varargin)+1);

% newbdy and oldbdy are circular lists = # unique points + 1
nnew = length(newbdy)-1; 
nold = length(oldbdy)-1; 


if nnew > nold 
    % make oldbdy longer
    nlong = nnew;
    b1 = interpboundary(oldbdy,nlong+1); 
    b2 = newbdy;
    % store additional arguments (keep circular)
    for i = 1:length(varargin)
        tmpvarsin{i} = varargin{i};
    end
elseif nnew < nold 
    % make newbdy longer, and all other bdy quantities longer using interpolation
    nlong = nold;
    b1 = oldbdy;
    b2 = interpboundary(newbdy,nlong+1); 
    for i = 1:length(varargin)
        [~,tmpvarsin{i}] = interpboundary(newbdy,nlong+1,varargin{i});
    end
else
    nlong = nnew; % same as nold bc lengths are equal in this case
    b1 = oldbdy;
    b2 = newbdy;
    % store additional arguments (keep circular)
    for i = 1:length(varargin)
        tmpvarsin{i} = varargin{i};
    end
end

dst = zeros(nlong,1);

% make bdy quantities non-circular (unique points only for interpolation and alignment)
b1u = b1(1:end-1,:);
b2u = b2(1:end-1,:);
for i = 1:length(varargin)
    tmpvarsinu{i} = tmpvarsin{i}(1:end-1);
end

for i=0:nlong-1
    % find distance b/w x,y coords of b2u and b1u (both same # points)
    % cycle through all possible rotations 
    % b1u,b2u are listed CW initially 
    d=b1u-circshift(b2u,i);% rotate initial point in b2u CCW i positions
    dst(i+1)=mean(sqrt(sum(d'.^2)));
end
% figure(72)
% plot(dst)
[~,minidx]=min(dst);
%minidx - 1 = how much b2 initial point should be rotated CCW (+circshift)

if nnew < nold
    % interpolation needed to downsample b2 back to original bdynew length
    b2align = circshift(b2u,minidx-1); % align b2 
    for i=1:length(varargin)
        tmpvarsout{i} = circshift(tmpvarsinu{i},minidx-1); %align other quantities
    end   
    % make all quantities circular again before interpolation
    b2alignC = [b2align; b2align(1,:)]; 
    for i=1:length(varargin)
        tmpvarsoutC{i} = [tmpvarsout{i}; tmpvarsout{i}(1)];
    end 
    % shorten b2alignC and other quants back to original newbdy length, circular output
    tmpvarsoutC_orig = cell(size(tmpvarsoutC));
    [b2alignC_orig, tmpvarsoutC_orig{:}] = interpboundary(b2alignC,nnew+1,tmpvarsoutC{:}); 

    varargout{1} = b2alignC_orig;
    for i = 1:length(varargin)
        varargout{1+i} = tmpvarsoutC_orig{i};
    end

else % nnew = nold OR nnew > nold
    % no interpolation needed bc b2align is already same length as bdynew
    b2align = circshift(b2u,minidx-1);
    varargout{1} = [b2align; b2align(1,:)]; % make circular
    for i=1:length(varargin)
        tmpvarsout{i} = circshift(tmpvarsinu{i},minidx-1);
        varargout{i+1} = [tmpvarsout{i}; tmpvarsout{i}(1)];
    end
end

end


function newvelT = velocity_motiontrack(bdy1,c1,bdy0,c0,dt,I1rgb,pveltrack)
% inputs:
% newbdy, oldbdy: [row col], circular lists (repetition of first point at end of list)
% newc, oldc: centroid positions [y0, x0]
% pveltrack: structure array

% computes velocity based on aligned boundaries (motion tracking)  
% create aligned/CENTERED equal length boundary vectors (without repetition of first point)
% rotate newbdy and compute sum of displacement map between the two ordered lists (newbdy and oldbdy)
% aligned newbdy = rotation of newbdy that minimizes the sum of displacements 

% computes net displacement = aligned disp + centroid disp
% vel = net disp / dt

% output: all aligned, circular, and interpolated back to original newbdy length if newbdy was adjusted
% varargout{2} = velocity based on motion tracking

v2struct(pveltrack);

newbdy = bdy1 - c1;
oldbdy = bdy0 - c0;

% newbdy and oldbdy are circular lists = # unique points + 1
nnew = length(newbdy)-1; 
nold = length(oldbdy)-1; 


if nnew > nold 
    % make oldbdy longer
    nlong = nnew;
    b1 = interpboundary(oldbdy,nlong+1); 
    b2 = newbdy;

elseif nnew < nold 
    % make newbdy longer, and all other bdy quantities longer using interpolation
    nlong = nold;
    b1 = oldbdy;
    b2 = interpboundary(newbdy,nlong+1); 
else
    nlong = nnew; % same as nold bc lengths are equal in this case
    b1 = oldbdy;
    b2 = newbdy;
end

dst = zeros(nlong,1);

% make bdy quantities non-circular (unique points only for interpolation and alignment)
b1u = b1(1:end-1,:);
b2u = b2(1:end-1,:);

for i=0:nlong-1
    % find distance b/w x,y coords of b2u and b1u (both same # points)
    % cycle through all possible rotations 
    % b1u,b2u are listed CW initially 
    d=b1u-circshift(b2u,i);% rotate initial point in b2u CCW i positions
    dst(i+1)=mean(sqrt(sum(d'.^2)));
end
% figure(72)
% plot(dst)
[~,minidx]=min(dst);
%minidx - 1 = how much b2 initial point should be rotated CCW (+circshift)

if nnew < nold
    % interpolation needed to downsample b2 back to original bdynew length

    b2ualign = circshift(b2u,minidx-1); % align b2 

    b2uorig = b2ualign + c1; % original coords (to account for centroid disp)
    b1uorig = b1u + c0; 

    d0 = (b2uorig - b1uorig) ; % NET displacements between two non-circular lists (new - old bdy same length)
    disp = sqrt(sum(d0'.^2))';

    %[in] = inpolygon(xq,yq,xv,yv) returns [in] indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    % 1 if xq, yq is INSIDE or ON the polygon defined by the query points
    [in3] = inpolygon(b2uorig(:,2),b2uorig(:,1),b1uorig(:,2),b1uorig(:,1));
    inout3 = 1-2*in3;
    svel3 = (1/dt)*disp.*inout3;
    
    % make all quantities circular again before interpolation
    b2uorigC = [b2uorig; b2uorig(1,:)]; 
    b1uorigC = [b1uorig; b1uorig(1,:)];
    velC = [svel3; svel3(1)];
    
    % shorten b2alignC and dispC back to original newbdy length, circular output
    [b2_norig,velC_norig] = interpboundary(b2uorigC,nnew+1,velC); 
    newvelT = velC_norig;
else    % nnew = nold OR nnew > nold
    % no interpolation needed bc b2align is already same length as bdynew
    b2ualign = circshift(b2u,minidx-1);

    b2uorig = b2ualign + c1; % original coords (to account for centroid disp)
    b1uorig = b1u + c0; 


    d0 = (b2uorig - b1uorig) ; % NET displacements between two non-circular lists (new - old bdy same length)
    disp = sqrt(sum(d0'.^2))';

    %[in] = inpolygon(xq,yq,xv,yv) returns [in] indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    % 1 if xq, yq is INSIDE or ON the polygon defined by the query points
    [in3] = inpolygon(b2uorig(:,2),b2uorig(:,1),b1uorig(:,2),b1uorig(:,1));
    inout3 = 1-2*in3;
    svel3 = (1/dt)*disp.*inout3;

    % make all quantities circular again and NO interpolation
    b2uorigC = [b2uorig; b2uorig(1,:)]; 
    b1uorigC = [b1uorig; b1uorig(1,:)];
    b2_norig = b2uorigC;
    velC = [svel3; svel3(1)];
    newvelT = velC;
end


if do_plotVelseries_motiontrack

    cols = size(I1rgb,2); rows = size(I1rgb,1);
    npts = length(b2uorigC);
    idx = 1:round(npts/50):npts;

    % use list of boundary points that are equal and before any last step interpolation 
    bdy = b2uorigC;
    obdy = b1uorigC;

   
    obdym = fliplr(obdy)'; % xcoords row1, ycoords row2
    bdym = fliplr(bdy)'; % xcoords row1, ycoords row2
    x12 = [obdym(1,idx); bdym(1,idx)];
    y12 = [obdym(2,idx);bdym(2,idx)];
    velCidx = velC(idx);
    
    f5 = figure(5);clf(f5);
        f5.Position(3:4) = [500,500];

    tf5 = tiledlayout(2,2); tf5.TileSpacing = 'compact'; tf5.Padding='compact';
    nexttile(1,[2,2]);
    imshow(0.1*I1rgb); axis equal; hold on;
    
    plot(bdy(:,2),bdy(:,1),'Color',[1 1 1],'LineWidth',5);
    plot(obdy(:,2),obdy(:,1),'Color',0.5*[1 1 1],'LineWidth',2,'LineStyle','--');


    for i = 1:length(idx)
        %figure(5);

        f5ax1 = surf([x12(:,i) x12(:,i)],...
            [y12(:,i) y12(:,i)], ...
            zeros(2,2),velCidx(i)*ones(2,2),'FaceColor','none','EdgeColor','flat',...
            'LineWidth',3);

       


    end
        
    colorbar; colormap(parula);
    f5ax1.Parent.CLim = [-5 5]; %caxis([-9, 9]);

    title('Boundary velocity (px/frame)')
    FontSize = minfont; %min(minfont,floor(min(cols,rows)/10));
    text(6,6,num2str(movieFrame),'Color','white','FontSize',FontSize);
        

    hold off; 
    drawnow


    %f5.Position(3:4) = [600,600];

    if savefigs>0 % create tiff series
        %pause(0.05)

        F5 = getframe(f5);
        if movieFrame == im0
            imwrite(F5.cdata,strcat(outputfolder,targname,'_',tstamp1,'_',filetag,'_tseries_vel_motiontrack.tif'));
        else
            imwrite(F5.cdata,strcat(outputfolder,targname,'_',tstamp1,'_',filetag,'_tseries_vel_motiontrack.tif'),...
                'WriteMode','append');
        end
    end


end % end of do_plotVelseries_motiontrack

end

function [I,Iorig]=readframe(fullFileName,frame)
    if ischar(fullFileName) % can read tif file
        Iorig = imread(fullFileName,frame);
    else % can input image directly
        Iorig = fullFileName(:,:,frame);
    end
    if ndims(Iorig)>2
        Iorig = Iorig(:,:,2);
    end
    I = imgaussfilt(Iorig,1);
    I = medfilt2(I,[1 1]);
end

function [svel3,idx] = velocity_pdist2(bdy,obdy,dt,pVel,varargin) % BSA 

    v2struct(pVel);
[dist, idx]  = pdist2(obdy,bdy,'euclidean','Smallest',1);
% dist contains K = 1 smallest pairwise distance to a point in obdy for each observation in bdy. 
% For each observation in Y, pdist2(X,Y,'smallest'/'largest',K) finds the K smallest or largest distances 
% by computing and comparing the distance values to all the observations in X. 

%idx contains: (row) index of point in X that was mapped to each element of Y.

dist  = dist';

% this option requires figure being opened
%oROI1  = drawfreehand('Position', fliplr(obdy),'Visible','off');
%inout1 = 1-2*inROI(oROI1,bdy(:,2),bdy(:,1));
%oROI2 = images.roi.Freehand('Position',fliplr(obdy));
%inout2 = 1-2*inROI(oROI2,bdy(:,2),bdy(:,1));
%inROI: points cannot lie on the roi, only inside/outside. (see 'classify pixels that partially enclosed')
%     hold off
%     drawnow
% vel = 1/dt * displacement (+ for outward, - for inward) 
% svel1 = (1/dt)*dist.*inout1; 

%[in] = inpolygon(xq,yq,xv,yv) returns [in] indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
% 1 if xq, yq is INSIDE or ON the polygon defined by the query points
[in3] = inpolygon(bdy(:,2),bdy(:,1),obdy(:,2),obdy(:,1));
inout3 = 1-2*in3;
svel3 = (1/dt)*dist.*inout3;

%     figure(1)

%     imshow(I1rgb); hold on; axis equal;
%     plot(bdy(:,2),bdy(:,1),'r--');
%     plot(obdy(:,2),obdy(:,1),'g--');
%     plot(bdy(100,2),bdy(100,1),'ro');
%     plot(obdy(100,2),obdy(100,1),'go');
%     hold off; 
%     drawnow

%% plot displacements

if do_plotVelseries_pdist2

    I1rgb = varargin{1};
    cols = size(I1rgb,2); rows = size(I1rgb,1);

    f5 = figure(5);clf(f5);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    f5.Position(3:4) = [500, 500];
    nexttile(1,[2,2]);

    obdym = fliplr(obdy(idx,:))'; % xcoords row1, ycoords row2
    bdym = fliplr(bdy)'; % xcoords row1, ycoords row2
    x12 = [obdym(1,:); bdym(1,:)];
    y12 = [obdym(2,:);bdym(2,:)];

    imshow(0.1*I1rgb); hold on; axis equal;
    plot(bdy(:,2),bdy(:,1),'Color',[1 1 1],'LineWidth',5);
    plot(obdy(:,2),obdy(:,1),'Color',0.5*[1 1 1],'LineWidth',2,'LineStyle','--');


    for i = 1:length(bdy)
       % figure(5);

        f5ax = surf([x12(:,i) x12(:,i)],...
            [y12(:,i) y12(:,i)], ...
            zeros(2,2),svel3(i)*ones(2,2),'FaceColor','none','EdgeColor','flat',...
            'LineWidth',3);


    end
    colorbar;    colormap(parula);
    f5ax.Parent.CLim = [-5 5]; %caxis([-9, 9]);
    title('Boundary velocity (px/frame)')
    FontSize = min(minfont,floor(min(cols,rows)/10));
    text(6,6,num2str(movieFrame),'Color','white','FontSize',FontSize);

    hold off; 
    drawnow
    %f5.Position(3:4) = [600,600];

    if savefigs > 0 % create tiff series
        %pause(0.05)
        F5 = getframe(f5);
        if movieFrame == im0
            imwrite(F5.cdata,strcat(outputfolder,targname,'_',tstamp1,'_',filetag,'_tseries_vel_pdist2.tif'));
        else
            imwrite(F5.cdata,strcat(outputfolder,targname,'_',tstamp1,'_',filetag,'_tseries_vel_pdist2.tif'),...
                'WriteMode','append');
        end
   end


end % end of do_plotVelseries_pdist2

end


function [getFluor,getiBW,getoBW,getinner,getouter] = getIntensity(bdy,Iinit,BW,pfluor)
% I will be a cell array based on fluorescent channels


v2struct(pfluor); % extract parameters

% if bits == 16
%     maxfl = 65535;
% end

getFluor = cell(1,nfluor); %intensities
getiBW = cell(1,nfluor);
getoBW  = cell(1,nfluor);
getinner = cell(1,nfluor);
getouter = cell(1,nfluor);

count = 0;
for k = 1:nchannels
    if ismember(k,fluorch)
    count = count + 1;

    %% **SETUP** adjust if needed
    del   = 4;% 4 originally # of points left or right used to find perpendicular line segment 
    lperp = nsize*8;% px outer and px inner for perpendicular line segment
    pin = nsize*8; % px to erode for luminal mask
    pout1 = nsize*8;  % px to dilate for 1st outer ring of external mask
    pout2 = nsize*16; % px to dilate for 2nd outer ring of external mask 
    ptop = [1]; % prop of top gfp intensities to use along lperp 
    % KEEP ptop = 1 for most targets; can adjust lower as needed; 
    % USE PLPERP  to adjust how far outward to compute top fluorescent intensities
    npts  = length(bdy); 
    fluor   = zeros(npts,1); % for ptop1
    linner = zeros(npts,2);
    louter = zeros(npts,2);
    iBW  = imerode(BW,strel('disk',round(pin),4));% 6,4
    oBW1 = imdilate(BW,strel('disk',round(pout1),4));% 6,4
    oBW2 = imdilate(BW,strel('disk',round(pout2),4));% 16,4
    oBW  = immultiply(imcomplement(oBW1),oBW2);
    Ifl    = Iinit{k}; % single channel fluorescent signal
    Iflrgb = repmat(imadjust(Ifl),1,1,3);
    cols = size(Ifl,2); rows = size(Ifl,1);  
    cent = centroid;

    for pt=1:npts
        x0  = bdy(pt,2);
        y0  = bdy(pt,1);
        cb  = circshift(bdy,del+1-pt);
        seg = cb(1:2*del+1,[2 1]);
        [z, r, ~] = fitcircle(seg); % z center, r radius
        if(isnan(r))
            r=rold;
            z=zold;
        end
        xc = z(1); yc = z(2);
        % find normal to curve and extend line segment through in polar coords
        theta = atan2(yc-y0,x0-xc); % dy = yc - row# ; dx = col# - xc
        d = -lperp:0.5:lperp;% 0.5 increment originally
        
        x = x0+d*cos(theta); % x(1) = x0 - lperp*cosA ---> where x0 is the boundary point
        y = y0-d*sin(theta);% because row coord increases going down

        [~,idx]= min( [norm([x(1),y(1)]-cent), norm([x(end),y(end)]-cent)]  );
        lend = round(plperp*length(x)); % end point of line segment to consider
        if idx == 1
            xord = x(1:lend); 
            yord= y(1:lend);
            linner(pt,:)=[xord(1),yord(1)]; % first row of inner is coords of innermost point of line segment (d=-lperp)
            louter(pt,:)=[xord(end),yord(end)]; % first row of outer is coords of outermost point of line segment (d=+lperp)

        else % idx == 2 
            xa = fliplr(x);
            ya = fliplr(y);
            xord = xa(1:lend)';
            yord = ya(1:lend)';
            linner(pt,:)=[xord(1),yord(1)]; % first row of inner is coords of innermost point of line segment (d=-lperp)
            louter(pt,:)=[xord(end),yord(end)]; % first row of outer is coords of outermost point of line segment (d=+lperp)

        end
            
        rold = r;
        zold = z;

        tmptmp  = improfile(Ifl,xord,yord,'bilinear');
        %figure(33); imshow(imadjust(I)); hold on;  plot(x,y,'r','LineWidth',5); hold off;
        %drawnow
        tmp     = sort(tmptmp,'descend','MissingPlacement','last');
     
        if do_top3px
            % mean of just top 3 pixels (less sensitive to spatial thickening)
            fluor(pt) = mean(tmp(1:3));
        else
            % mean
            fluor(pt) = mean(tmp(1:ptop*end));  % top percentage to use    
            % median
            %fluor(pt) = double(median(tmp(1:ptop*end)));  % top percentage to use    
        end
    end
    
    
    

    % membrane gfp normalized to background and internal signal
    % (gfp_memb - mean(background)) / (mean(gfp_internal) - mean(background))
    %gfp = (gfp-mean(I(oBW)))/(mean(I(iBW))-mean(I(oBW)));
    
    % sort and take mean of bottom fraction of cytosolic (ibw) and external (obw)
    tmpibw = sort(Ifl(iBW),'MissingPlacement','last');% ascending list of gfp values in mask iBW
    tmpobw = sort(Ifl(oBW),'MissingPlacement','last');% ascending list of gfp values in mask oBW
     mibw = mean(tmpibw(1:round(0.2*end)));
     mobw = mean(tmpobw(1:round(0.2*end)));
    %mibw = double(median(tmpibw(1:round(0.2*end))));
    %mobw = double(median(tmpobw(1:round(0.2*end))));

%     if k > 1
%        fluor = (fluor-mobw)/(mibw-mobw);
%     elseif k == 1
%         fluor = (fluor-mobw)/maxfl;
%     end

 %% Option for lumen signal division
 
    if do_lumennorm
        % OPTION 1: use LUMEN division here:
        fluor = (fluor-mobw)/(mibw-mobw);
    else
        % OPTION 2: NO LUMEN DIVISION
        fluor = (fluor-mobw);
    end



%%
    getoBW{count} = oBW;
    getiBW{count} = iBW;
    getFluor{count} = fluor;
    getinner{count} = linner;
    getouter{count} = louter;

        %% Option for membrane marker normalization
    if k > 1 && do_membnorm
        getFluor{count} = fluor./getFluor{1}; % divide by membrane channel
    end

  

    
    end % if ismember k fluorch
end    %end of nfluors 

end

function [allkymos, alignall,varargout] = makekymos(dataRaw,kf)
% allkymos : cell array, 1 by nbdyquants 
%   each entry is a matrix of the kymograph to be plotted
%   rows = nbins or nframes, cols = 6 (vel, curv, dcumul, 3 fluor channels

% alignall : cell array, nframes by (nbdyquants + 1)
%   each entry is a vector of aligned boundary quantities ... 
%   (except entries in column 1: boundary points have 2 columns)
%   rows = nboundarypts or nbins or nframes, 
%   cols = 1 usually, except boundary points which have 2

% alignall_long : same as alignall, but using full nboundarypoints always
% this will be varargout, added in case where binning occurs

% kymofixedpts = dataRaw.kymofixedpts;
% kymofixedtheta = dataRaw.kymofixedtheta;

targname = dataRaw.targname;
bdyall       = dataRaw.allBoundaries; % cell array nframes by 1 boundary points
bdycurv      = dataRaw.allCurvatures; % cell array nframes by 1
bdyvel       = dataRaw.allpdist2Velocities; % cell array nframes by 1
bdypd2dcumul    = dataRaw.allpdist2Dcumuls;
bdyperim     = dataRaw.allPerims;
fluorall       = dataRaw.allIntensities; % cell array nframes by nfluor
cent         = dataRaw.allCentroids; % matrix: nframes by 2
leng         = dataRaw.allLengths; % vector: nframes by 1
im0          = dataRaw.im0;
imL          = dataRaw.imL;
nframes      = dataRaw.nframes; % note nframes does not always equal imL+1-im0 
nfluor       = dataRaw.nfluor; % number of fluorescent channels 
membch       = dataRaw.membch;
fluorch      = dataRaw.fluorch;
listFrames   = dataRaw.listFrames; % may not be consecutive frames during testing
thetacenter     = dataRaw.thetacenter; % angle for center of kymograph
fnames     = dataRaw.fnames;
savefigs     = dataRaw.savefigs;
outputfolder = dataRaw.outputfolder;
tstamp1    = dataRaw.tstamp1;

do_plotVelseries_motiontrack = dataRaw.do_plotVelseries_motiontrack;
nbins = dataRaw.nbins;
pvary = dataRaw.pvary;
allIrgb = dataRaw.allIrgb;
veltrack1 = dataRaw.veltrack1;  % fixedkymo: 0 for pdist2, 1 for motiontracking 
veltrack2 = dataRaw.veltrack2;  % variable kymo: 0 for pdist2, 1 for motiontracking
dt = dataRaw.dt;
deltaf = dataRaw.deltaf;
do_movavg = dataRaw.do_movavg;
s_movavg = dataRaw.s_movavg;
minfont = dataRaw.minfont;

% if no smoothing average, set s to 0
if ~do_movavg
    s_movavg = 0;
end

% new variables
alignbdy  = cell(nframes,1);
alignvel  = cell(nframes,1);
aligncurv = cell(nframes,1);
alignfluor  = cell(nframes,nfluor);
aligndcumul = cell(nframes,1);


nbdyquants = 3+nfluor; % vel, curv, dcumul kymos + fluor kymos  
%alignall_long = cell(nframes,nbdyquants+1); % col: bdypositions, vel, curv, dcumul, fluor kymos

if kf > 0 %kymofixedpts || kymofixedtheta 


    n = max(leng); % max # of boundary points in circular list over all frames 
    alignangles = cell(nframes,1);
    alignDtheta = cell(nframes,1);
    alignbinidx = cell(nframes,1); % bdy point bin indices for every frame
    if kf ==2 %kymofixedtheta
        edges = - pi : 2*pi/nbins : pi ; 
    else
        edges = 1 : (n-1)/nbins : n;
    end
    varargout{1} = cell(nframes,3); %for alignangles, alignDtheta, alignbinidx
    varargout{2} = edges';
else
    nvary = round(pvary*bdyperim);
end




 

%--------------------------
for j=1:nframes 
    
    movieFrame = listFrames(j);
    % snake interpolation if having fixed # of points or theta bins in kymograph
    if kf > 0 %kymofixedpts || kymofixedtheta    
        fluor = cell(1,nfluor); % pre-allocation is required   
        % interpolate bdy quantities
        % pass fluor as comma separated list, not as a cell array
        [tmp,vel,curv,pd2dcumul,fluor{:}] = interpboundary(bdyall{j},n,bdyvel{j},...
            bdycurv{j},bdypd2dcumul{j},fluorall{j,:});
        % interpBdy and interpCurv are circular; others just very close
        curv(end) = curv(1);
        vel(end) = vel(1);
        pd2dcumul(end) = pd2dcumul(1);
        % make list of fluor circular and keep as cell array
        matfluor = cell2mat(fluor);
        matfluor(end,:) = matfluor(1,:);
        fluor = mat2cell(matfluor,length(fluor{1}),ones(1,nfluor));
        % current centroid
        x0  = cent(j,1); % centroids are (col, row), (x, y)
        y0  = cent(j,2);
    else
       fluor = cell(1,nfluor); % pre-allocation is required   
        % interpolate bdy quantities
        [tmp,vel,curv,pd2dcumul,fluor{:}] = interpboundary(bdyall{j},nvary(j),bdyvel{j},...
            bdycurv{j},bdypd2dcumul{j},fluorall{j,:});
        % interpBdy and interpCurv are circular; others just very close
        curv(end) = curv(1);
        vel(end) = vel(1);
        pd2dcumul(end) = pd2dcumul(1);
        % make list of fluor circular and keep as cell array
        matfluor = cell2mat(fluor);
        matfluor(end,:) = matfluor(1,:);
        fluor = mat2cell(matfluor,length(fluor{1}),ones(1,nfluor));
        % current centroid
        x0  = cent(j,1); % centroids are (col, row), (x, y)
        y0  = cent(j,2);

%         tmp = bdyall{j}; % (row col) (y, x)
%         x0  = cent(j,1); % centroids are (col, row), (x, y)
%         y0  = cent(j,2);
%         vel = bdyvel{j}; 
%         curv = bdycurv{j};
%         pd2dcumul = bdypd2dcumul{j};
%         fluor = fluorall(j,:); %cell array
    end

    %--------------------------------------------------------------------
    % We now align the boundaries. If this is the first frame, we realign
    % it so that the angle 'theta0' corresponds roughly to the point
    % that is at the center of the kymograph (assuming symmetric changes). 
    % If this is not the first frame, we align it so that the points are closest 
    % to the previous frame, relative to the centroid.
    %--------------------------------------------------------------------
    %  thetapi = angle (rad) of initial point in aligned list, will be pi away from thetacenter
    thetastar = wrapToPi(deg2rad(thetacenter)); % reference angle, center of kymograph 
    thetapi= wrapToPi(thetastar+pi); % top of kymograph, 
    % make output from top to bottom y-axis: +pi to 0 to -pi (boundary points always go clockwise)
    newfluor = cell(size(fluor)); % pre-allocate is required here
    if j==1
        [newtmp,newvel,newcurv,newpd2dcumul,newfluor{:}] = aligninitboundary(allIrgb{j,1},tmp,[x0 y0],thetapi,vel,curv,pd2dcumul,fluor{:});

        if kf > 0 %kymofixedpts || kymofixedtheta
            newdcumul = zeros(size(newvel));
        else % kymo varying, initial frame always 0
            newdcumul = newpd2dcumul;
        end

    else % after first frame, align/track based on previous frame
        % newbdy: (y, x) relative to new centroid
        % oldbdy: (y, x) relative to old centroid

        [newtmp,newvel,newcurv,newpd2dcumul,newfluor{:}] = alignboundary_ba(tmp-[y0 x0],oldtmp-[y0old x0old],vel,curv,pd2dcumul,fluor{:});
%         if veltrack
%             pveltrack=v2struct(movieFrame,do_plotVelseries,im0,savefigs,outputfolder,tstamp1,targname,minfont);
%             [newvel] = velocitytrack(tmp,[y0 x0],oldtmp,[y0old x0old],dt,allIrgb{j,membch},pveltrack);        
%         end

        if kf > 0 && veltrack1 % (kf > 0 means kymofixedpts || kymofixedtheta)
            filetag = 'kymofixed';
            do_plotVelseries_motiontrack = 0; % no need to plot velocity for both kymofixed and kymovarying
            pveltrack=v2struct(filetag,movieFrame,do_plotVelseries_motiontrack,im0,savefigs,outputfolder,tstamp1,targname,minfont);
            [newvel] = velocity_motiontrack(tmp,[y0 x0],oldtmp,[y0old x0old],dt,allIrgb{j,membch},pveltrack);        
            newdcumul = olddcumul + dt*newvel;

        elseif kf > 0 && ~veltrack1 % (kf > 0 means kymofixedpts || kymofixedtheta)
            newdcumul = olddcumul + dt*newvel;

        elseif kf == 0 && veltrack2  % kymo varying use pmotion tracking for vel
            filetag = 'kymovary';
            pveltrack=v2struct(filetag,movieFrame,do_plotVelseries_motiontrack,im0,savefigs,outputfolder,tstamp1,targname,minfont);
            [newvel] = velocity_motiontrack(tmp,[y0 x0],oldtmp,[y0old x0old],dt,allIrgb{j,membch},pveltrack);        
            % can only use pdist2 for dcumul
            newdcumul = newpd2dcumul;

        elseif kf == 0 && ~veltrack2 % kymovarying and no motiontrack for vel
            % use pdist2 for both velocity and dcumul
            % newvel already computed above
            newdcumul = newpd2dcumul;
        end
        newtmp = newtmp+[y0 x0];
    end
    
    x0old = x0;
    y0old = y0;
    
    if kf > 0 % kymofixedpts || kymofixedtheta
    olddcumul = newdcumul;
    end
    %--------------------------------------------------------------------
    % We are ready to save this frame's boundary. - angle computations, moving average
    %--------------------------------------------------------------------
    
    % these are all circular lists

    if kf > 0 %kymofixedpts || kymofixedtheta 
        % compute angles for every boundary point in frame j
        % compute displacement relative to thetastar (center of kymograph)

        if kf == 2 %kymofixedtheta
            %figure(10); imshow(allIrgb{1,1},[]); hold on; plot(tmp(1:50:200,2),tmp(1:50:200,1),'ro');plot(tmp(1,2),tmp(1,1),'g*');hold off; drawnow
            %figure(11); imshow(allIrgb{1,1},[]); hold on; plot(newtmp(1:50:200,2),newtmp(1:50:200,1),'ro');plot(newtmp(1,2),newtmp(1,1),'g*');hold off; drawnow
            % deltay = y0 - bdyrow; delta x = bdycol - x0;
            angles = atan2(y0-newtmp(:,1),newtmp(:,2)-x0); % tmp [row, col] 
            D = wrapToPi(angles - thetastar); % angles relative to thetastar
            alignangles{j} = angles; % raw angles from -pi to pi
            alignDtheta{j} = D; % angles relative to thetastar
        else
            D = 1:n; %indices of points will be binned  based on edges defined above for kymofixedpts
        end

        % if no moving average, s = 0;
        snewfluor = cell(size(newfluor));
        [snewvel, snewcurv,snewdcumul,snewfluor{:}] = mavgsmoothboundary(s_movavg,newvel,newcurv,newdcumul,newfluor{:});
        
        
        %alignbdy{j}  = newtmp;
        binfluor = cell(1,nfluor);
        %[binidx, binvel, bincurv, bindcumul,binfluor{:}] = binbdyquant(kymofixed,D,edges,snewvel,snewcurv,snewdcumul,snewfluor{:});
        [binidx, bintmprow, bintmpcol, binvel, bincurv, bindcumul,binfluor{:}] = ...
            binbdyquant(kf,D,edges,newtmp(:,1),newtmp(:,2),snewvel,snewcurv,snewdcumul,snewfluor{:});
        alignbdy{j}  = [bintmprow,bintmpcol];
        % ============ for kymofixedtheta or fixedpts
        % we want to align quantities to be centered around a constant thetastar - done using D
        % for kymofixedpts
        % just want to bin quantities in fewer points
        alignbinidx{j} = binidx;
        alignvel{j} = binvel;
        aligncurv{j} = bincurv;
        aligndcumul{j} = bindcumul;
        alignfluor(j,:)  = binfluor; % cell array, row: frame j, col: fluor channels
        
        testfornans = any(isnan(bincurv));
        if testfornans
            kf
            j
            disp('some NaNs are present');
            %return
        end

        oldtmp       = newtmp; 
        
    else % # of points in kymograph columns vary as perimeter changes


        % if no moving average, s = 0;
        snewfluor = cell(size(newfluor));
        [snewvel, snewcurv,snewdcumul,snewfluor{:}] = mavgsmoothboundary(s_movavg,newvel,newcurv,newdcumul,newfluor{:});

        alignbdy{j}  = newtmp;  
        alignvel{j}  = snewvel;
        aligncurv{j} = snewcurv;
        aligndcumul{j} = snewdcumul;
        alignfluor(j,:)  = snewfluor; % cell array, row: frame j, col: fluor channels
        oldtmp       = newtmp; 

        %alignall_long(j,[1:3,5:nbdyquants]) = [{newtmp},{snewvel},{snewcurv},snewfluor];

    end


end % end of nframes



%compute cumulative displacement
if kf > 0 %kymofixedtheta || kymofixedpts
%     v1 = cell2mat(alignvel');
%     d1 = dt*v1;
%     dcumul = zeros(size(v1));
%     dcumul(:,2:end) = cumsum(d1(:,2:end),2);
%     aligndcumul = mat2cell(dcumul,[nbins],ones(1,nframes))'; % transpose to match column format
     varargout{1} = [alignangles,alignDtheta,alignbinidx];

end

%alignall_long(:,4) = aligndcumul; % will be empty for when usual kymograph where npoints varies 



allkymos = cell(1,nbdyquants); % every element is a matrix (max #boundary points by nframes)

alignall = cell(nframes,nbdyquants+1); % aligned  vel, curv, dcumul + fluor 
% (every element is different size matrix as # boundary points change)
alignall(:,1:4) = [alignbdy,alignvel, aligncurv, aligndcumul];
alignall(:,5:end) = alignfluor;

% ignore repeat for circ lists and circshift being used later for smoothing
ky = cell(1,nbdyquants); % temporary variable for frame j
for i = 1:nbdyquants

    % use deltaf to include blank columns in kymograph for skipped frames
    if kf > 0 %kymofixedtheta || kymofixedpts
        allkymos{i} = nan(nbins,deltaf);
        lkym = nbins*ones(nframes,1); % kymograph rows = nbins , fixed
        ell = zeros(nframes,1); % no padding needed
    else
        % using max(leng) will plot the repeat at end of circular lists
        % keeping it circular makes no difference bc so hard to see repeat 
        allkymos{i} = nan(max(nvary),deltaf); 
        lkym      = nvary; % kymograph rows varies frame to frame
        % ell: half of difference between current bdy length and max bdy length
        ell       = floor((max(nvary)-nvary)/2); % adjust every frame
    end

    for j = 1:nframes 
         % frame j, i+1 gives bdyquantity of interest: keep circ fluor list; or all bins if doing binning
        ky{i} = alignall{j,i+1};
        % use ell to center data
        jj = listFrames(j)-im0+1; % column index accounting for skipped frames
        allkymos{i}((ell(j)+1):(ell(j)+lkym(j)),jj)=ky{i}+0;   
    end
end % end of nbdyquants



end % function

function varargout = binbdyquant(kf,Dtheta,edges,varargin)
nbins = length(edges)-1;
Y = discretize(Dtheta,edges); % assign bin indices for every Dtheta, fixed for a given frame
nvar = length(varargin);
varargout = cell(1,nvar+1);

if kf == 2
    %bins go from -pi to pi (based on edges); 
    % shift to be indexed from +pi to -pi
    Y1 = nbins-Y+1; % assign bin indices in reverse order 
    varargout{1} = Y1;
else
    Y1 = Y;
    varargout{1} = Y1;
end


for i = 1:nvar
    binX = nan(nbins,1);
    X = varargin{i};

    for j = 1:nbins
        % find all data points assigned to bin index j
        binidx = (Y1==j); % same for every varargin, but changes placement of 1s as bins change
        % assign mean of data points in binidx to proper row in binned data
        binX(j) = mean(X(binidx));

    end

    varargout{i+1} = binX;
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

function [z, r, residual] = fitcircle(x, varargin)
%FITCIRCLE    least squares circle fit
%   
%   [Z, R] = FITCIRCLE(X) fits a circle to the N points in X minimising
%   geometric error (sum of squared distances from the points to the fitted
%   circle) using nonlinear least squares (Gauss Newton)
%       Input
%           X : 2xN array of N 2D points, with N >= 3
%       Output
%           Z : center of the fitted circle
%           R : radius of the fitted circle
%
%   [Z, R] = FITCIRCLE(X, 'linear') fits a circle using linear least
%   squares minimising the algebraic error (residual from fitting system
%   of the form ax'x + b'x + c = 0)
%
%   [Z, R] = FITCIRCLE(X, Property, Value, ...) allows parameters to be
%   passed to the internal Gauss Newton method. Property names can be
%   supplied as any unambiguous contraction of the property name and are 
%   case insensitive, e.g. FITCIRCLE(X, 't', 1e-4) is equivalent to
%   FITCIRCLE(X, 'tol', 1e-4). Valid properties are:
%
%       Property:                 Value:
%       --------------------------------
%       maxits                    positive integer, default 100
%           Sets the maximum number of iterations of the Gauss Newton
%           method
%
%       tol                       positive constant, default 1e-5
%           Gauss Newton converges when the relative change in the solution
%           is less than tol
%
%   [X, R, RES] = fitcircle(...) returns the 2 norm of the residual from 
%   the least squares fit
%
%   Example:
%       x = [1 2 5 7 9 3; 7 6 8 7 5 7];
%       % Get linear least squares fit
%       [zl, rl] = fitcircle(x, 'linear')
%       % Get true best fit
%       [z, r] = fitcircle(x)
%
%   Reference: "Least-squares fitting of circles and ellipses", W. Gander,
%   G. Golub, R. Strebel - BIT Numerical Mathematics, 1994, Springer

% This implementation copyright Richard Brown, 2007, but is freely
% available to copy, use, or modify as long as this line is maintained

error(nargchk(1, 5, nargin, 'struct'))

% Default parameters for Gauss Newton minimisation
params.maxits = 100;
params.tol    = 1e-5;

% Check x and get user supplied parameters
[x, fNonlinear, params] = parseinputs(x, params, varargin{:});

% Convenience variables
m  = size(x, 2);
x1 = x(1, :)';
x2 = x(2, :)';


% 1) Compute best fit w.r.t. algebraic error using linear least squares
% 
% Circle is represented as a matrix quadratic form
%     ax'x + b'x + c = 0
% Linear least squares estimate found by minimising Bu = 0 s.t. norm(u) = 1
%     where u = [a; b; c]

% Form the coefficient matrix
B = [x1.^2 + x2.^2, x1, x2, ones(m, 1)];

% Least squares estimate is right singular vector corresp. to smallest
% singular value of B
[U, S, V] = svd(B);
u = V(:, 4);

% For clarity, set the quadratic form variables
a = u(1);
b = u(2:3);
c = u(4);

% Convert to centre/radius
z = -b / (2*a);
r = sqrt((norm(b)/(2*a))^2 - c/a);

% 2) Nonlinear refinement to miminise geometric error, and compute residual
if fNonlinear
    [z, r, residual] = fitcircle_geometric(x, z, r);
else
    residual = norm(B * u);
end

% END MAIN FUNCTION BODY

% NESTED FUNCTIONS
    function [z, r, residual] = fitcircle_geometric(x, z0, r0)
        % Use a simple Gauss Newton method to minimize the geometric error
        fConverged = false;
        
        % Set initial u
        u     = [z0; r0];
        
        % Delta is the norm of current step, scaled by the norm of u
        delta = inf;
        nIts  = 0;
        
        for nIts = 1:params.maxits
            % Find the function and Jacobian 
            [f, J] = sys(u);
            
            % Solve for the step and update u
            h = -J \ f;
            u = u + h;
            
            % Check for convergence
            delta = norm(h, inf) / norm(u, inf);
            if delta < params.tol
                fConverged = true;
                break
            end
        end
        
        if ~fConverged
            warning('fitcircle:FailureToConverge', ...
                'Gauss Newton iteration failed to converge');
        end
        z = u(1:2);
        r = u(3);
        f = sys(u);
        residual = norm(f);
        
        
        function [f, J] = sys(u)
            %SYS   Nonlinear system to be minimised - the objective
            %function is the distance to each point from the fitted circle
            %contained in u

            % Objective function
            f = (sqrt(sum((repmat(u(1:2), 1, m) - x).^2)) - u(3))';
            
            % Jacobian
            denom = sqrt( (u(1) - x1).^2 + (u(2) - x2).^2 );
            J = [(u(1) - x1) ./ denom, (u(2) - x2) ./ denom, repmat(-1, m, 1)];
        end % sys
        
    end % fitcircle_geometric

% END NESTED FUNCTIONS

end % fitcircle
function [x, fNonlinear, params] = parseinputs(x, params, varargin)
% Make sure x is 2xN where N > 3
if size(x, 2) == 2
    x = x';
end

if size(x, 1) ~= 2
    error('fitcircle:InvalidDimension', ...
        'Input matrix must be two dimensional')
end

if size(x, 2) < 3
    error('fitcircle:InsufficientPoints', ...
        'At least 3 points required to compute fit')
end

% determine whether we are measuring geometric error (nonlinear), or
% algebraic error (linear)
fNonlinear = true;
switch length(varargin)
    % No arguments means a nonlinear least squares with defaul parameters
    case 0
        return
       
    % One argument can only be 'linear', specifying linear least squares
    case 1
        if strncmpi(varargin{1}, 'linear', length(varargin{1}))
            fNonlinear = false;
            return
        else
            error('fitcircle:UnknownOption', 'Unknown Option')
        end
        
    % Otherwise we're left with user supplied parameters for Gauss Newton
    otherwise
        if rem(length(varargin), 2) ~= 0
            error('fitcircle:propertyValueNotPair', ...
                'Additional arguments must take the form of Property/Value pairs');
        end

        % Cell array of valid property names
        properties = {'maxits', 'tol'};

        while length(varargin) ~= 0
            property = varargin{1};
            value    = varargin{2};
            
            % If the property has been supplied in a shortened form, lengthen it
            iProperty = find(strncmpi(property, properties, length(property)));
            if isempty(iProperty)
                error('fitcircle:UnkownProperty', 'Unknown Property');
            elseif length(iProperty) > 1
                error('fitcircle:AmbiguousProperty', ...
                    'Supplied shortened property name is ambiguous');
            end
            
            % Expand property to its full name
            property = properties{iProperty};
            
            switch property
                case 'maxits'
                    if value <= 0
                        error('fitcircle:InvalidMaxits', ...
                            'maxits must be an integer greater than 0')
                    end
                    params.maxits = value;
                case 'tol'
                    if value <= 0
                        error('fitcircle:InvalidTol', ...
                            'tol must be a positive real number')
                    end
                    params.tol = value;
            end % switch property
            varargin(1:2) = [];
        end % while

end % switch length(varargin)

end %parseinputs
function choice = timeoutdlg(delay)
    prompt='Manual over ride [y/n]?';
    f1 = findall(0, 'Type', 'figure');
    delay
    t = timer('TimerFcn', {@closeit f1}, 'StartDelay', delay);
    start(t);
    % Call the dialog
    retvals = inputdlg(prompt,'',1,{'n'});
    if numel(retvals) == 0
          choice = 'n';
    else
          choice = retvals{1};
    end
    % Delete the timer
    if strcmp(t.Running, 'on')
           stop(t);
    end
    delete(t);
    
    function closeit(src, event, f1)
        f2 = findall(0, 'Type', 'figure');
        fnew = setdiff(f2, f1);
        if ishandle(fnew);
            close(fnew);
        end
    end
end

%% 

function [shape_curvature] = curvature1_v5(Bk,Pcurv1)

%%***********************************************************************%
%*                         Curvature measure                            *%
%*                  Measure shape properties of loops.                  *%
%*                                                                      *%
%* Original author: Dr. Meghan Driscoll, Preetham Manjunatha                                *%
%* Modified by:    Bedri Abubaker-Sharif                                  *%
%************************************************************************%
%
% Usage: [shape, Icurv] = curvature(image, boundaryPoint, curvatureThresh, ...
%                                     bp_tangent, interp_resolution, loopclose)
% Inputs: % Inputs:
% binaryImage: input binary image
% Bk:   for closed curves, this should be a circular list of spaced boundary points 
% Pcurv1: packed variables for settings;

% Outputs: 
%           shape                
%           .curvature          - the boundary curvature at each boundary point (uses snakeNum) 
%                                 Curvatures above or below a cutoff are given the magnitude of the cutoff
%           .meanNegCurvature   - the mean negative curvature
%           .numIndents         - the number of boundary regions over which the curvature is negative
%           .tangentAngle       - the angle of the tangent to the boundary at each boundary point
%           .tortuosity         - the boundary tortuousity (a measure of how bendy the boundary is)
%           .uncutCurvature     - the uncut boundary curvature at each boundary point (uses snakeNum)

%           Icurv               - Output image (padded image to make 3 channel. This fixes the plot)
%--------------------------------------------------------------------------


% unpack structure variables
v2struct(Pcurv1);

% % Perimeter of the binary component
% stats = regionprops(binaryimage,'perimeter');
% perimeter = cat(1,stats(1).Perimeter);

shape_XY = fliplr(Bk); % convert to (x,y) (col, row)
% shape_XY is circular

if loopclose == 1 % loop is closed
    % if there are 5 unique bp: 1 2 3 4 5
    % shape_XY: 1 2 3 4 5 1
    % bp_positions: 5 1 2 3 4 5 1

    % Nb is # of unique boundary points; M is length of shape_XY
    Nb = size(shape_XY,1) - 1;    
    bp_positions = [shape_XY(Nb,:); shape_XY(1:Nb+1,:)]; % Mth term is already a repeat of 1st boundary point
    
else % loop is open
    % if there are 5 unique bp: 1 2 3 4 5
    % repeat neighboring elements for end points
    % shape_XY: 1 2 3 4 5 
    % bp_positions: 2 1 2 3 4 5 4

    % Nb is # of unique boundary points; M is length of shape_XY
    Nb = size(shape_XY,1) ;

    %bp_positions = [shape_XY(2,:); shape_XY(1:M,:); shape_XY(M-1,:)];
    bp_positions = [shape_XY(2,:); shape_XY(1:Nb,:); shape_XY(Nb-1,:)];
end

% initialize variables    
shape_curvature         = NaN(Nb,1);
shape_uncutCurvature    = NaN(Nb,1);
shape_meanNegCurvature  = NaN;
shape_numIndents = NaN;
shape_tortuosity = NaN;
shape_tangentAngle = NaN(Nb,1);

% calculate the curvature (by finding the radius of the osculating circle using three neaby boundary points)
for j = 1:Nb

    % assign the three points that the circle will be fit to such that the slopes are not infinite 
    p1 = bp_positions(j,:);
    p2 = bp_positions(j+1,:); pstar = p2;
    p3 = bp_positions(j+2,:); 
    s12 = (p1(2)-p2(2))/(p1(1)-p2(1));
    s23 = (p2(2)-p3(2))/(p2(1)-p3(1));
    
    % if exactly one of the two slopes is infinite, then relabel points
    if isinf(abs(s12)) && ~isinf(abs(s23))
        % switch point 2,3 for slope computations
        point0 = p2; p2 = p3; p3 = point0;
        s12 = (p1(2)-p2(2))/(p1(1)-p2(1));
        % this slope does not change
        s23 = (p2(2)-p3(2))/(p2(1)-p3(1));    
    end

    if isinf(abs(s23)) && ~isinf(abs(s12))
        % switch point 1,2 for slope computations
        point0 = p1; p1 = p2; p2 = point0;
        % this slope does not change
        s12 = (p1(2)-p2(2))/(p1(1)-p2(1));
        % this slope changes
        s23 = (p2(2)-p3(2))/(p2(1)-p3(1));    
    end

    % if the boundary is flat (s12, s23 both 0, both Inf, both -Inf)
    if s12==s23  
        shape_curvature(j) = 0;

    % if the boundary is curved
    else

        % calculate the curvature
        % mp12: midpoint between line segment 12 
        % mp23: midpoint between "     "      23
        % mp13: midpoint between "     "      13
        mp12 = (p1+p2)/2;
        mp23 = (p2+p3)/2;
        mp13 = (p1+p3)/2;

        % x coordinate of center of circumscribed circle
        xc = -s23*(mp12(1))/(-s23+s12) + s12*mp23(1)/(-s23+s12) + s12*s23*(mp23(2)-mp12(2))/(-s23+s12);
        
        if isequal(s12,0)
            yc = (-1/s23)*(xc-mp23(1))+mp23(2);
        else
            yc = (-1/s12)*(xc-mp12(1))+mp12(2);
        end
        shape_uncutCurvature(j) = 1/sqrt((p2(1)-xc)^2+(p2(2)-yc)^2);

         shape_curvature(j) = shape_uncutCurvature(j);
         % cutoff the curvature (for visualization)
%         if shape_curvature(j) > curvatureThresh
%             shape_curvature(j) = curvatureThresh;
%         end

        % determine if the curvature is positive or negative
%         [In, On] = inpolygon(mp13(1),mp13(2),shape_XY(:,1),shape_XY(:,2)); 
% 
%         if ~In              
%             shape_curvature(j) = -1*shape_curvature(j);
%             shape_uncutCurvature(j) = -1*shape_uncutCurvature(j);
%         end

        % determine if the curvature is positive or negative
        % compute tiny line segment from perimeter towards center of fitted circle
        % find if point lies outside the boundary (curvature is negative)  
        x0 = pstar(1); y0 = pstar(2);
        xp = x0 + 1e-3*(xc-x0); %order of substraction shows direction 
        yp = y0 + 1e-3*(yc-y0);
        In = inpolygon(xp,yp,shape_XY(:,1),shape_XY(:,2)); % 1 if xp,yp is inside or on polygon
        if ~In              
            shape_curvature(j) = -1*shape_curvature(j);
            shape_uncutCurvature(j) = -1*shape_uncutCurvature(j);;
        end 

    end 

end

% find the mean negative curvature (really this should use a constant dist snake)
listCurve = shape_uncutCurvature(1:Nb);
listNegCurve = abs(listCurve(listCurve < 0));
if ~isempty(listNegCurve) 
    shape_meanNegCurvature = sum(listNegCurve)/(Nb);
else
    shape_meanNegCurvature = 0;
end

% find the number of negative boundary curvature regions
curveMask = (listCurve < 0);
curveMaskLabeled = bwlabel(curveMask);
numIndents = max(curveMaskLabeled);
if curveMask(1) && curveMask(end)
    numIndents  = numIndents-1;
end
shape_numIndents = numIndents;

% find the tortuosity (should correct units)
shape_tortuosity = sum(gradient(shape_uncutCurvature(1:Nb)).^2)/perimeter;

% calculate the direction of the tangent to the boundary 
if loopclose == 1
    %Nb = M-1
    %bp_positions_tangent = [shape_XY(M-1,:); shape_XY(1:M,:)];
    bp_positions_tangent = [shape_XY(Nb,:); shape_XY(1:Nb+1,:)];
else
    %Nb = M
    %bp_positions_tangent = [shape_XY(2,:); shape_XY(1:M,:); shape_XY(M-1,:)];
    bp_positions_tangent = [shape_XY(2,:); shape_XY(1:Nb,:); shape_XY(Nb-1,:)];
end
for j=1:Nb
    pt1 = bp_positions_tangent(j,:);
    pt2 = bp_positions_tangent(j+2,:); 
    shape_tangentAngle(1,j) = mod(atan2(pt1(2)-pt2(2), pt1(1)-pt2(1)), pi);
end

% shape.XY = shape_XY;
% shape.curvature         = shape_curvature;
% shape.uncutCurvature    = shape_uncutCurvature;
% shape.meanNegCurvature  = shape_meanNegCurvature;
% shape.numIndents = shape_numIndents;
% shape.tortuosity = shape_tortuosity;
% shape.tangentAngle = shape_tangentAngle;

end

%--------------------------------------------------------------------------
function [xi,yi] = snakeinterp_ba1(x,y,RES)
% SNAKEINTERP1  Interpolate the snake to have equal distance RES
%     [xi,yi] = snakeinterp(x,y,RES)
%
%     RES: resolution desired

%     update on snakeinterp after finding a bug
%      Chenyang Xu and Jerry L. Prince, 3/8/96, 6/17/97
%      Copyright (c) 1996-97 by Chenyang Xu and Jerry L. Prince
%      image Analysis and Communications Lab, Johns Hopkins University

% convert to column vector
% make sure x, y are circular lists (last entry is repeat of first)

x = x(:); y = y(:);

N = length(x);  

dx = diff(x);
dy = diff(y);
d = sqrt(dx.*dx+dy.*dy);  % compute the distance from previous node for point 2:N+1

d = [0;d];   % point 1 to point 1 is 0 

% now compute the arc length of all the points to point 1
% we use matrix multiply to achieve summing 
M = length(d);
csumd = cumsum(d);

% now ready to reparametrize the closed curve in terms of arc length
perim = csumd(M);

if (perim/RES<3)
   error('RES too big compare to the length of original curve');
end

di = [0:RES:perim]';

xi = interp1(csumd,x,di);
yi = interp1(csumd,y,di);

N = length(xi);

if (perim - di(length(di)) <RES/2)  % deal with end boundary condition
   xi = xi(1:N-1);
   yi = yi(1:N-1);
end

% make the list circular
xi = [xi;xi(1)];
yi = [yi;yi(1)];
end


function C = plotcurvature1_v6(Irgb,boundaries,Bkcurv,Bkpoints,iflag) 

% raw or SG smoothed boundary points; circular lists
Borigs = boundaries; % (row , col )
Xorig = Borigs(:,2); % col indices of original boundary points
Yorig = Borigs(:,1); % row indices of original boundary points

%% Plot curvature
% Bkpoints, Bkcurve are circular lists



if iflag == 1 
    Xspaced = Bkpoints(1:end-1,2); % x coords for spaced boundary points used for curvature computation
    Yspaced = Bkpoints(1:end-1,1); % y coords for spaced boundary points used for curvature computation
    Cspaced = Bkcurv(1:end-1); % spaced curvatures 
    
    % %%% OPTION 1 scattered interpolation to original boundary points
    F = scatteredInterpolant(Xspaced,Yspaced,Cspaced);
    C = F(Xorig,Yorig); % circular list of curvatures interpolated to original boundary points

else
    Xspaced = Bkpoints(1:end-1,2); % x coords for spaced boundary points used for curvature computation
    Yspaced = Bkpoints(1:end-1,1); % y coords for spaced boundary points used for curvature computation
    Cspaced = Bkcurv(1:end-1); % spaced curvatures 
    
    %%% OPTION 2
    % Bkpoints(double) contains fewer boundary points than Borigs(indices; double if smoothed/interpolated)
    % Want to interpolate a curve with the same number of boundary points as Borigs
    % problem: poor X,Y b/c we're only utilizing sparse Bkpoints to reconstruct boundary; 
    % we're not including any info on Borigs, just the number of points in Borigs
    % so just ignore boundary output (yxinfo)
    % note: kinfo is circular list; same # of elements as Xorig, Yorig
    % for circular data: n = # of desired unique points + 1
    [~,kinfo] = interpboundary(Bkpoints,size(Borigs,1),Bkcurv); 
    C = kinfo; 

end

% size of data markers in scatter3
dmark = 45;



%% overlay interpolated curvature on original and interpolated/sampled boundary points
% figure(20); clf;
% s1 = subplot(1,2,1); axis equal; hold on
% title('k-circum interpolated along raw boundary points')
% imshow(Irgb)
% scatter3(Xorig,Yorig,zeros(size(Xorig)),[],C,'filled','SizeData',dmark)
% 
% cmap = parula; colormap(cmap);
% colorbar; caxis([-0.1 0.1]); 
% hold off
% 
% s2 = subplot(1,2,2); axis equal; hold on
% title('k-circum along smooth,spaced boundary points')
% imshow(Irgb)
% scatter3(Xspaced,Yspaced,zeros(size(Xspaced)),[],Cspaced,'filled','SizeData',dmark)
% 
% cmap = parula; colormap(cmap);
% colorbar; caxis([-0.1 0.1]);
% hold off


%% 

% figure(21); clf;
% subplot(2,1,1)
% plot(C); title('interpolated curvature (k)')
% xlabel('k-circum interpolated along raw boundary points')
% 
% subplot(2,1,2)
% plot(Cspaced); title('raw curvature (k)')
% xlabel('k-circum along smooth,spaced boundary points')

% figure; 
% [x1,y1,z1] = peaks;
% surf(x1,y1,z1)
end


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

function Bkpoints = defBkpoints(I,Borigs,defBk,Bk_dist,Bk_num)
%% boundary point adjustments for curvature estimation 
% Bk_num = # of UNIQUE boundary points (non-circular)
% define spaced boundary points for curvature estimation 
% Borig is (y,x) circular list or approximately circular (if using sgolay smoothed)
% Bkpoints: boundary points (row,col) for curvature computation
% assume closed loop: Bkpoints is a circular list
if strcmp(defBk,'distance')       
    [xi1,yi1] = snakeinterp_ba1(Borigs(:,2),Borigs(:,1),Bk_dist);
    Bkpoints = [yi1,xi1]; 
elseif strcmp(defBk,'number')    
    % will be circular 
    % n = # unique + 1;
    Bkpoints = interpboundary(Borigs,Bk_num+1); 
    % adjust X,Y coords if last point is too close to first point
    % only relevant for Bk points where we want equally spaced boundary
    pfirst = Bkpoints(1,:); psecond = Bkpoints(2,:); plast = Bkpoints(end,:);
    gaptest = norm(pfirst-plast);
    gap1 = norm(pfirst-psecond);
    if gaptest < 1e-8 % Bkpoints is perfectly circular, just finalize
        Bkpoints(end,:) = Bkpoints(1,:);
    elseif gaptest < 0.5*gap1
        Bkpoints = Bkpoints(1:end-1,:);
        Bkpoints = [Bkpoints; Bkpoints(1,:)];
    end
else
    Bkpoints = Borigs;
end

% f10 = figure(10); clf(f10);
% imshow(I);
% axis equal, hold on
% plot(Bkpoints(1,2),Bkpoints(1,1),'r*','LineWidth',10)
% plot(Bkpoints(:,2),Bkpoints(:,1),'b*','LineWidth',3)
% title('Smoothed and spaced boundary')
% hold off
% drawnow

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
            sbdy = sbdy/(2*s+1);
            varargout{k} = sbdy;
        end    
        % figure()
        %plot(bdy(1:200,1),bdy(1:200,2),'x',sbdy(1:200,1),sbdy(1:200,2),'.r')
    else
        varargout = varargin;
    end
end

function K = curve_fitcircle(bdy,delta)
% INPUT
% bdy: list of (y,x) or (row,col) positions of boundary points (circular list)
% delta is # of points to the right/left to consider fitting

% OUTPUT
% C1: local curvatures based on fitting a circle that minimizes geometric error
% same length as bdy (circular list)

% make equal spaced boundary segments for consistency to delta around the perimeter
n = size(bdy,1); % keep same length as original boundary
bdyI = interpboundary(bdy,n); % circular

B = fliplr(bdyI); % now it's [x, y] coords
Pn = B(1:end-1,:); % only uniqe boundary points considered for circshift to work properly

curve = zeros(length(Pn),1);

for i=1:length(Pn)
    s=circshift(Pn,delta-i+1);%not delta-i 
    % Fit to a circle the 2*delta+1 points centered at perimeter index "i"
    [zz,rr]=fitcircle(s(1:2*delta+1,:)'); % not -1
    xc = zz(1); yc = zz(2);     % circle center coords (x, y)
    x0 = Pn(i,1); y0 = Pn(i,2); % boundary index coords (x, y)

    if isnan(rr) % rr is NaN; 
        % so curvature is 1/Inf = zero
        curve(i)=0;
    else
        curve(i) = 1/rr; % unsigned curvature
        % determine if the curvature is positive or negative
        % compute tiny line segment from perimeter towards center of fitted circle
        % find if point lies outside the boundary (curvature is negative)  

        xp = x0 + 1e-3*(xc-x0); %order of substraction shows direction 
        yp = y0 + 1e-3*(yc-y0);
        In = inpolygon(xp,yp,Pn(:,1),Pn(:,2)); % 1 if xp,yp is inside or on polygon
        if ~In              
            curve(i) = -curve(i);
        end  
    end
end 


% make circular
K = [curve; curve(1,:)];


%ncurve=curve/(abs(mean(curve)));
end


%---------------
function C2 = plotcurvature2_v1(Irgb,boundaries,c2a,Bkpoints,iflag) 
% varargin when using spaced points: Bkpoints, iflag
% otherwise only 3 inputs


% raw or SG smoothed boundary points; circular lists
Borigs = boundaries; % (row , col )
Xorig = Borigs(:,2); % col indices of original boundary points
Yorig = Borigs(:,1); % row indices of original boundary points

%% Plot curvature
% Bkpoints, Bkcurve are circular lists


if nargin > 3
    if iflag == 1 
        Xspaced = Bkpoints(1:end-1,2); % x coords for spaced boundary points used for curvature computation
        Yspaced = Bkpoints(1:end-1,1); % y coords for spaced boundary points used for curvature computation
        C2spaced = c2a(1:end-1); % spaced curvatures 
        
        % %%% OPTION 1 scattered interpolation to original boundary points
        F = scatteredInterpolant(Xspaced,Yspaced,C2spaced);
        C2 = F(Xorig,Yorig); % circular list of curvatures interpolated to original boundary points
    
    else
        Xspaced = Bkpoints(1:end-1,2); % x coords for spaced boundary points used for curvature computation
        Yspaced = Bkpoints(1:end-1,1); % y coords for spaced boundary points used for curvature computation
        C2spaced = c2a(1:end-1); % spaced curvatures 
        
        %%% OPTION 2
        % Bkpoints(double) contains fewer boundary points than Borigs(indices; double if smoothed/interpolated)
        % Want to interpolate a curve with the same number of boundary points as Borigs
        % problem: poor X,Y b/c we're only utilizing sparse Bkpoints to reconstruct boundary; 
        % we're not including any info on Borigs, just the number of points in Borigs
        % so just ignore boundary output (yxinfo)
        % note: kinfo is circular list; same # of elements as Xorig, Yorig
        % for circular data: n = # of desired unique points + 1
        [~,kinfo] = interpboundary(Bkpoints,size(Borigs,1),c2a); 
        C2 = kinfo; 
    
    end
else
    C2 = c2a; % stay circular
end
% size of data markers in scatter3
dmark = 45;

%% overlay interpolated curvature on original and interpolated/sampled boundary points
% figure(40); 
% s1 = tiledlayout('flow'); s1.Padding='compact';s1.TileSpacing='compact';
% nexttile(1,[1 3]); axis equal; hold on
% if nargin>3
%     title('k-fitcircle interpolated along spaced boundary points')
% else
%     title('k-fitcircle along boundary points')
% end
% imshow(Irgb)
% f40ax1 = scatter3(Xorig,Yorig,zeros(size(Xorig)),[],C2,'filled','SizeData',dmark);
% cmap = parula; colormap(cmap);
% colorbar; f40ax1.Parent.CLim = [-0.22 0.22]; %caxis([-0.03 0.03]); 
% hold off
% 
% if nargin>3
% nexttile(4,[1 3]); axis equal; hold on
% title('k-fitcircle along spaced boundary points')
% imshow(Irgb)
% f40ax2 = scatter3(Xspaced,Yspaced,zeros(size(Xspaced)),[],C2spaced,'filled','SizeData',dmark);
% cmap = parula; colormap(cmap);
% colorbar; f40ax2.Parent.CLim = [-0.03, 0.03]; %caxis([-0.03 0.03]);
% hold off
% end

%% 

% figure(41); 
% s2 = tiledlayout('flow'); s2.Padding='compact';s2.TileSpacing='compact';
% nexttile(1,[1 3]);
% plot(C2); 
% if nargin > 3
%     title('interpolated curvature (k)')
%     xlabel('k-fitcircle interpolated along raw boundary points')
% else
%     title('curvature (k-fitcircle)')
%     xlabel('k-fitcircle along raw boundary points')
% end
% 
% if nargin>3
%     nexttile(4, [1 3]);
%     plot(C2spaced); title('raw curvature (k)')
%     xlabel('k-fitcircle along smooth,spaced boundary points')
% end

% figure; 
% [x1,y1,z1] = peaks;
% surf(x1,y1,z1)
end

function [A , c] = MinVolEllipse(P, tolerance)
% [A , c] = MinVolEllipse(P, tolerance)
% Finds the minimum volume enclsing ellipsoid (MVEE) of a set of data
% points stored in matrix P. The following optimization problem is solved: 
%
% minimize       log(det(A))
% subject to     (P_i - c)' * A * (P_i - c) <= 1
%                
% in variables A and c, where P_i is the i-th column of the matrix P. 
% The solver is based on Khachiyan Algorithm, and the final solution 
% is different from the optimal value by the pre-spesified amount of 'tolerance'.
%
% inputs:
%---------
% P : (d x N) dimnesional matrix containing N points in R^d.
% tolerance : error in the solution with respect to the optimal value.
%
% outputs:
%---------
% A : (d x d) matrix of the ellipse equation in the 'center form': 
% (x-c)' * A * (x-c) = 1 
% c : 'd' dimensional vector as the center of the ellipse. 
% 
% example:
% --------
%      P = rand(5,100);
%      [A, c] = MinVolEllipse(P, .01)
%
%      To reduce the computation time, work with the boundary points only:
%      
%      K = convhulln(P');  
%      K = unique(K(:));  
%      Q = P(:,K);
%      [A, c] = MinVolEllipse(Q, .01)
%
%
% Nima Moshtagh (nima@seas.upenn.edu)
% University of Pennsylvania
%
% December 2005
% UPDATE: Jan 2009



%%%%%%%%%%%%%%%%%%%%% Solving the Dual problem%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% ---------------------------------
% data points 
% -----------------------------------
[d N] = size(P);

Q = zeros(d+1,N);
Q(1:d,:) = P(1:d,1:N);
Q(d+1,:) = ones(1,N);


% initializations
% -----------------------------------
count = 1;
err = 1;
u = (1/N) * ones(N,1);          % 1st iteration


% Khachiyan Algorithm
% -----------------------------------
while err > tolerance,
    X = Q * diag(u) * Q';       % X = \sum_i ( u_i * q_i * q_i')  is a (d+1)x(d+1) matrix
    M = diag(Q' * inv(X) * Q);  % M the diagonal vector of an NxN matrix
    [maximum j] = max(M);
    step_size = (maximum - d -1)/((d+1)*(maximum-1));
    new_u = (1 - step_size)*u ;
    new_u(j) = new_u(j) + step_size;
    count = count + 1;
    err = norm(new_u - u);
    u = new_u;
end



%%%%%%%%%%%%%%%%%%% Computing the Ellipse parameters%%%%%%%%%%%%%%%%%%%%%%
% Finds the ellipse equation in the 'center form': 
% (x-c)' * A * (x-c) = 1
% It computes a dxd matrix 'A' and a d dimensional vector 'c' as the center
% of the ellipse. 

U = diag(u);

% the A matrix for the ellipse
% --------------------------------------------
A = (1/d) * inv(P * U * P' - (P * u)*(P*u)' );


% center of the ellipse 
% --------------------------------------------
c = P * u;
end

function Ellipse_plot(A, C, varargin)
%
%  Ellipse_Plot(A,C,N) plots a 2D ellipse or a 3D ellipsoid 
%  represented in the "center" form:  
%               
%                   (x-C)' A (x-C) <= 1
%
%  A and C could be the outputs of the function: "MinVolEllipse.m",
%  which computes the minimum volume enclosing ellipsoid containing a 
%  set of points in space. 
% 
%  Inputs: 
%  A: a 2x2 or 3x3 matrix.
%  C: a 2D or a 3D vector which represents the center of the ellipsoid.
%  N: the number of grid points for plotting the ellipse; Default: N = 20. 
%
%  Example:
%  
%       P = rand(3,100);
%       t = 0.001;
%       [A , C] = MinVolEllipse(P, t)
%       figure
%       plot3(P(1,:),P(2,:),P(3,:),'*')
%       hold on
%       Ellipse_plot(A,C)
%  
%
%  Nima Moshtagh
%  nima@seas.upenn.edu
%  University of Pennsylvania
%  Feb 1, 2007
%  Updated: Feb 3, 2007

%%%%%%%%%%%  start  %%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 20; % Default value for grid

% See if the user wants a different value for N.
%----------------------------------------------
if nargin > 2
 	N = varargin{1};
end


% check the dimension of the inputs: 2D or 3D
%--------------------------------------------
if length(C) == 3,
    Type = '3D';
elseif length(C) == 2,
    Type = '2D';
else
    display('Cannot plot an ellipse with more than 3 dimensions!!');
    return
end

% "singular value decomposition" to extract the orientation and the
% axes of the ellipsoid
[U D V] = svd(A);

if strcmp(Type, '2D'),
    % get the major and minor axes
    %------------------------------------
    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));

    theta = [0:1/N:2*pi+1/N];

    % Parametric equation of the ellipse
    %----------------------------------------
    state(1,:) = a*cos(theta); 
    state(2,:) = b*sin(theta);

    % Coordinate transform 
    %----------------------------------------
    X = V * state;
    X(1,:) = X(1,:) + C(1);
    X(2,:) = X(2,:) + C(2);
    
elseif strcmp(Type,'3D'),
    % generate the ellipsoid at (0,0,0)
    %----------------------------------
    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));
    c = 1/sqrt(D(3,3));
    [X,Y,Z] = ellipsoid(0,0,0,a,b,c,N);
    
    %  rotate and center the ellipsoid to the actual center point
    %------------------------------------------------------------
    XX = zeros(N+1,N+1);
    YY = zeros(N+1,N+1);
    ZZ = zeros(N+1,N+1);
    for k = 1:length(X),
        for j = 1:length(X),
            point = [X(k,j) Y(k,j) Z(k,j)]';
            P = V * point;
            XX(k,j) = P(1)+C(1);
            YY(k,j) = P(2)+C(2);
            ZZ(k,j) = P(3)+C(3);
        end
    end
end


% Plot the ellipse
%----------------------------------------
if strcmp(Type,'2D'),
    plot(X(1,:),X(2,:),'c','LineWidth',2);
    hold on;
    plot(C(1),C(2),'m*','MarkerSize',10,LineWidth=3);
    axis equal
    %grid
else
    mesh(XX,YY,ZZ);
    axis equal
    %hidden off
end

end


function [aspratio,A0,A0c,a0majax] = aspectratio(bdyraw)
    % prm0 = v2struct(minfont,movieFrame,do_plotEllipseries)
    %v2struct(prm0);

    %byraw = rows - bdyraw(:,1) + 1;
    byraw = bdyraw(:,1);
    bxraw = bdyraw(:,2);
    xpts = [bxraw,byraw];
    %xmean = xpts-mean(xpts);
%     if do_plotEllipseries
%         f2 = figure(2); imshow(Irgb); axis equal; hold on;
%         plot(xpts(:,1),xpts(:,2),'r','LineWidth',3); 
%         FontSize = minfont; %min(minfont,floor(min(cols,rows)/10));
%         text(6,6,num2str(movieFrame),'Color','white','FontSize',FontSize)
% 
%         set(gca,'YDir','reverse');
%         
%         [A0,A0c] = MinVolEllipse(xpts',0.01);
%         a0 = sqrt(1/min(A0(1),A0(4)));
%         b0 = sqrt(1/max(A0(1),A0(4)));
%         eccent = sqrt(1-(b0/a0)^2);
%         aspratio = a0/b0;
%         %stats.Eccentricity
%     
%         Ellipse_plot(A0,A0c);
%         varargout{1} = f2;
%     else
        [A0,A0c] = MinVolEllipse(xpts',0.01);
        a0 = sqrt(1/min(A0(1),A0(4)));
        b0 = sqrt(1/max(A0(1),A0(4)));
        eccent = sqrt(1-(b0/a0)^2);
        aspratio = a0/b0;

        a0majax = 2*a0;

% 
%     end

end

function v = recordvid(currentframe,vidparams)
% save captured frame as avi

% recall vidparams = v2struct(fname,vidrate,basystem)
v2struct(vidparams) 


if strcmp('linux',basystem)
    v = VideoWriter(fname,'Motion JPEG AVI');
else
    v = VideoWriter(fname,'MPEG-4');
end

v.FrameRate = vidrate; %playback framerate
open(v);
  

axis equal    
writeVideo(v,currentframe.cdata);
%pause(0.1)
close(v);


end