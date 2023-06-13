format compact
clearvars, clc
close all
warning('off','all')
%set(0,'DefaultFigureWindowStyle','normal') 


%% SETTINGS 
%--------------------------------------------------------------------------   
aligncase       = 0; % TOGGLE THIS TO DO ALIGNMENT 
submaskcase     = 0; % TOGGLE THIS TO DO SUB-MASKING
savefigs        = 1;
fixAllFrames    = 0; % when false: specify in frames2fix 
guvname         = 'Global6_DPPC';
listFrames      = [1:89];   % choose frames to include in analysis, in ascending order

% inputfolder     = strcat(guvname,'_input/L04-30_raw/');
% outputfolder    = strcat(guvname,'_input/',guvname,'_acs_v0/'); % name of output folder

inputfolder     = strcat('Supp_',guvname,'_input/');
outputfolder    = strcat('Supp_',guvname,'_input/processed/'); % name of output folder


% ALWAYS place membrane marker in fname1 or reference frame
% fname1          = strcat(guvname,'_membprox'); % membprox
% fname2          = strcat(guvname,'_actin'); 
% fname3          = strcat(guvname,'_ActAprox'); % ActAprox
% fnames          = {fname1;fname2;fname3};

fname1          = strcat(guvname,'_mCh'); % membrane marker file name
fname2          = strcat(guvname,'_YFP'); % actin marker file name
fname3          = strcat(guvname,'_CFP');
fnames          = {fname1;fname2;fname3};

% if strcmp(guvname(1),'L')
%     fname1          = strcat(guvname,'_CmCh'); % membrane marker file name
%     fname2          = strcat(guvname,'_CYFP'); % actin marker file name
%     fname3          = strcat(guvname,'_CCFP');
%     % fname4          = strcat(guvname,'_BF');
%     % fname5          = strcat(guvname,'_647_dye');
%     fnames          = {fname1;fname2;fname3};
%     needlech            = 4 ; % BF channel for needle
% elseif strcmp(guvname(1),'G')
%     fname1          = strcat(guvname,'_CmCh'); % membrane marker file name
%     fname2          = strcat(guvname,'_CYFP'); % actin marker file name
%     fname3          = strcat(guvname,'_CCFP');
%     fnames          = {fname1;fname2;fname3};
% end

maskchannel     = 1; % channel number of mask
maskfname       = strcat(guvname,'_mask');
initmag         = 400; 


% Pick one of these - default is 'circle'
    fixFrameCommand = 'circle';
    %fixFrameCommand = 'polygon';
    %fixFrameCommand = 'rectangle';

% ---------------------------------------------------------
nchannels       = length(fnames);

fullFileNames   = cell(1,nchannels);
allr0           = cell(1,nchannels);
allr1           = cell(1,nchannels);
for i = 1:nchannels
    fullFileNames{i}  = strcat(inputfolder,fnames{i},'_aligned.tif');
    if ~isfile(fullFileNames{i})
        fullFileNames{i} = strcat(inputfolder,fnames{i},'.tif');
    end
end


tstamp1 = string(datetime('now','Format','dd_MMM_yyyy_HH_mm_ss'));
% gather general info from channel 1
finfo        = imfinfo(fullFileNames{1});
nallframes   = length(finfo);
initcols         = finfo.Width; % full cols
initrows         = finfo.Height;% full rows
bits         = finfo(1).BitDepth; 

% temporal cropping
im0          = listFrames(1);         
imL          = min(nallframes,listFrames(end));       
nframes      = length(listFrames);


% frames for segmentation
if fixAllFrames
    frames2fix = [1:30];   
    frames2fix_second = [1:30];
else
    frames2fix = [1:5]; % specify which frames to fix here
    %frames2fix = [1]; % specify which frames to fix here
    frames2fix_second = [1:5];
end

delay = 10; % delay for manual checkings


if exist(outputfolder,'dir') == 0
    mkdir(outputfolder);
end

%--------------------------------------------------------------------------

varpass = v2struct;

%% align channels
% case 1: one channel is off from others: 
%      keep dimensions of recttempcrop same, just shift what you access from original image coordinates
%      constant new rectpos over all frames in one channel only
% case 2: all channels shift over time:
%      keep dimensions of recttempcrop same, just shift what you access
%      variable new rectpos over all frames, but same over all channels
% case 3: want to center a moving cell
%      keep dimensions of recttempcrop same, just shift wh
% at you access
%      variable new rectpos over all frames, but same over all channels


switch aligncase
    case 0
        % no alignment needed
        frows = initrows; fcols = initcols;

    case 1

        % choose target frame for channel (only 1 frame will be used for alignment)
        % compute mask for channels
        % find centroid xci,yci (ith channel) and compare disp (x,y) relative to reference frame centroid
        % adjust rectpos for misaligned channels based on centroid displacement relative to good channels
        % make rectpos into matrix with length nframes
        refchan = 1; 
        alignchan = 2; % single channel to align
        targlist = 34; % frame number to use for alignment
        [allCent,~] = centercoords(alignchan,refchan,targlist,varpass);
        
        % fix any other channels
        alignrestchans = [5];
%         allDisp{alignrestchans} = allDisp{alignchan};
        allCent{alignrestchans} = allCent{alignchan};
        
%         D = cell2mat(allDisp);
%         Dsmall = (D <= 0.01*min(initrows,initcols));
%         D(Dsmall) = 0;
%     
%         allDispadj = mat2cell(D,ones(1,nchannels),[2]);

        allImages_align = alignframes(allCent,alignchan,alignrestchans,...
            refchan,varpass);
        frows = size(allImages_align{1,maskchannel},1);
        fcols = size(allImages_align{1,maskchannel},2);
end

%% draw recttemp crop:
% initial temporal and rect crop, initial mask
% rectangular crop initial based on membrane marker first frame
if aligncase == 0
    [I1gs0,I1orig0] = readframe(fullFileNames{maskchannel},1); %uint16 output (not rgb)    
else
    %make sure frows and fcols have been updated to use indices of aligned images
    maskch_align = cat(3,allImages_align{:,maskchannel});
    [I1gs0,I1orig0] = readframe(maskch_align,1); %uint16 output (not rgb)    
end

I1rectp = rectcropThisFrame(I1gs0);
% xmin ymin width height
I1x0 = max(1,round(I1rectp(1)));
I1x1 = min(fcols,round(I1x0+I1rectp(3)));
I1y0 = max(1,round(I1rectp(2)));
I1y1 = min(frows,round(I1y0+I1rectp(4)));

rectrows = I1y1 - I1y0+1;
rectcols = I1x1 - I1x0+1;

r0 = repmat([I1x0, I1y0],nframes,1);
r1 = repmat([I1x1, I1y1],nframes,1);

for j = 1:nchannels
    allr0{j} = r0;
    allr1{j} = r1;
end


%% final mask and images
% offer chance to define a sub-mask starting from allBW_init
% allBW_init and Xinit have same dimensions based on recttempcrop

if submaskcase
    % recall
    % listFrames = list of temporal cropping
    % nframes = length of listFrames
    
    % generate initial rectangular crop from previous drawrectangle --> allImages_init
    % select initial mask automatically or manually inside rectangular crop
    % save corresponding binarization
    % temporal cropping: frame numbers in original tiff
    varpass.frames2fix = frames2fix;
    saveseg = 0;
    [allImages_init, allBW_init] = segmentframes(listFrames,rectrows,...
        rectcols,allr0,allr1,saveseg,varpass);

    % can use allBW_init to force later sub-mask (secondary segmentation) to only consider desired region, 
    % it's a way to force a boundary not to extend past desired region,
    % even if background signal or artifacts come nearby
    saveseg = 1;
    varpass.frames2fix = frames2fix_second;
    [allImages, allBW] = segmentframes(listFrames,rectrows,rectcols,allr0,allr1,...
        saveseg,varpass,allBW_init);
else
    saveseg = 1;
    [allImages, allBW] = segmentframes(listFrames,rectrows,rectcols,allr0,allr1,...
        saveseg,varpass);
end


%% -----------
%% Auxiliary functions
function [allCent,allDisp] = centercoords(alignchan,refchan,list,varpass,varargin)
% input: 
% alignchan,refchan: select channels to find centroid (at least 2: 1 reference, 1 test)
% list: frames to find centroids
% varpass: general setting variables from main workspace
% varargin: optional

% output:
% allDisp: displace centroid (x,y) for every channel (will be zeros for unselected channels)
% refCent: reference centroid position (x,y) - rounded to grid integers

v2struct(varpass)
nlist = length(list);

%---------------------------------------------

% initialize variables for centroid calculation
allDisp = cell(1,nchannels); % displacements
allCent = cell(1,nchannels); % centroid 
for i = 1:nchannels
    allDisp{i} = NaN(nlist,1); % channel i, only TARGET frames
    allCent{i} = NaN(nlist,2); % centroid only for TARGET frames
end

for i = 1:nlist

    for k = 1:nchannels   
    
    if k == refchan || ismember(k,alignchan)    
        %---------------------------------------------
        movieFrame=list(i); % TIF frame indices for temporal cropping 
        frame = i; % indices for saving data
    
        [Xgs,Xorig] = readframe(fullFileNames{k},movieFrame); %uint16 output (not rgb)    
        % define initial images for computations
        % segment based on I (gaussian smoothed membrane maker) 
        % or segment based on Iorig (original membrane marker)
        Xinit = Xgs;
        
        % make boundary marker into rgb image for plotting
        % use imadjust to increase contrast
        Xrgb = repmat(imadjust(Xinit),1,1,3);
        %----------------------------------
                  
        % autosegment without oBW bc no previous frame
        if k == 1 %refchan
            BW        = segment(Xinit);   
            BWref = BW;
        else
            BW = segment(Xinit,BWref);
        end
        f4 = figure(4); clf(f4);
        f4.Position(1:2) = [300,600];
        imshow(Xinit,[],'InitialMagnification',initmag);
        axis equal, hold on
        visboundaries(BW,'Color','y');
        choice = timeoutdlg(delay);
        if ((choice~='n')&&(choice~='N'))&&(choice~=0) %((choice=='y')||(choice=='Y'))
            choice2 = secondtimeoutdlg(delay);
            BW = fixThisFrame(Xinit,BW,choice2);
        end
        title('initializing mask based on only first frame')
        hold off
        drawnow
    
        [~,L]     = bwboundaries(BW,'noholes');
        stats     = regionprops(L,'Centroid');
        allCent{k}(i,:)  = stats(1).Centroid;    
    end
    
    end % end of nchannels

end % end of nlist

for k = 1:nchannels
    % fill ref centroids for channels not in alignchan
    % will be adjusted later
    if ~ismember(k,alignchan) && k ~= refchan
        allCent{k} = allCent{refchan};
    end
end

for k = 1:nchannels
    allDisp{k} = allCent{k}-allCent{refchan};
end

end

function [allImages0] = alignframes(allCent,alignchan,alignrestchans,...
    refchan,varpass)

v2struct(varpass)

xcref = min(max(round(allCent{refchan}(1)),1),initcols); % b/w  1 and initcols
ycref = min(max(round(allCent{refchan}(2)),1),initrows); % b/w 1 and initrows
xcalign = min(max(round(allCent{alignchan}(1)),1),initcols); % b/w  1 and initcols
ycalign = min(max(round(allCent{alignchan}(2)),1),initrows); % b/w 1 and initrows

dx1 = min(xcref-1,xcalign-1); % cols to the left of center
dx2 = min(initcols-xcref,initcols-xcalign); % cols to the right
dy1 = min(ycref-1,ycalign-1); % rows above center point
dy2 = min(initrows-ycref,initrows-ycalign); % rows below

outrows = dy2+dy1+1;
outcols = dx2+dx1+1;
allImages0 = cell(nframes,nchannels);
% if bits==8
%     for i = 1:nchannels
%         ubits = 'uint8';
%         allImages0{i}  = zeros(outrows,outcols,ubits);   % channel i, all frames
%     end
% elseif bits==16
%     for i = 1:nchannels
%         ubits = 'uint16';
%         allImages0{i}  = zeros(outrows,outcols,nframes,ubits);   % channel i, all frames    
%     end
% elseif bits==32
%     for i = 1:nchannels
%         ubits = 'uint32';
%         allImages0{i}  = zeros(outrows,outcols,nframes,ubits);   % channel i, all frames    
%     end    
% end

for k = 1:nchannels
    for i = 1:nframes
        movieFrame = listFrames(i);

        if ismember(k,alignchan) || ismember(k,alignrestchans)
            [~,Im0] = readframe(fullFileNames{k},movieFrame); %uint16 output (not rgb)    
            Im1 = Im0(ycalign-dy1:ycalign+dy2,xcalign-dx1:xcalign+dx2);                    
        else
            [~,Im0] = readframe(fullFileNames{k},movieFrame); %uint16 output (not rgb)    
            Im1 = Im0(ycref-dy1:ycref+dy2,xcref-dx1:xcref+dx2);   
        end
        allImages0{i,k} = Im1;
        % save mask and image as tiff series if did segmentation in matlab
        if savefigs
            if movieFrame == im0 
                imwrite(Im1,strcat(outputfolder,fnames{k},'_aligned.tif'));
            else
                imwrite(Im1,strcat(outputfolder,fnames{k},'_aligned.tif'),'WriteMode','append');
            end
        end 
    end
end
end

function [allImages0,allBW0,allfused] = segmentframes(list,outrows,outcols,p0,p1,saveseg,...
    varpass,varargin)

% input: 
% list: temporally crop frames listed
% varpass: general setting variables from main workspace
% varargin: optional, row/col coordinates of upper left and lower right of rectangle for spatial cropping

% if varargin is empty, function will not do spatial cropping
% so it will ask only once for initial fixed mask based on first frame of mask channel (membrane marker preferred)
% if varargin is not empty, function will do spatial cropping

% output:
% allImages0: every channel, all frames in temporal crop, 
% allBW0: mask based on membrane channel, 

v2struct(varpass)
nlist = length(list);
submask = ~isempty(varargin);
x0 = cell(nchannels,1); y0 = cell(nchannels,1);
x1 = cell(nchannels,1); y1 = cell(nchannels,1);
if submask
    BW0 = varargin{1};
    for j = 1:nchannels
        x0{j} = p0{j}(:,1); y0{j} = p0{j}(:,2);
        x1{j} = p1{j}(:,1); y1{j} = p1{j}(:,2);
    end
%     outrows = min(y1{1}-y0{1}+1);%all identical so min is just picking one
%     outcols = min(x1{1}-x0{1}+1);%all identical so min is just picking one
else
    % use list of upperleft and bottom right coordinates
    for j = 1:nchannels
        x0{j} = p0{j}(:,1); y0{j} = p0{j}(:,2);
        x1{j} = p1{j}(:,1); y1{j} = p1{j}(:,2);
    end
%     outrows = min(y1{1}-y0{1}+1);%all identical so min is just picking one
%     outcols = min(x1{1}-x0{1}+1);%all identical so min is just picking one
end


%---------------------------------------------

% initialize variables for temporal cropping and initial segmentation
% use frows, fcols (not newrows, newcols) if no spatial cropping yet
allImages0 = cell(nlist,nchannels);
Igs0 = cell(nlist,nchannels);
Iorig0 = cell(nlist,nchannels);
Iinit0 = cell(nlist,nchannels);
Irgb0 = cell(nlist,nchannels);

allBW0 = zeros(outrows,outcols,nlist,'logical'); % Segmentation all temporal crop frames - only based on 1 channel
allfused = zeros(outrows,outcols,3,nlist);
% if bits==8
%     for i = 1:nchannels
%         ubits = 'uint8';
%         allImages0{i}  = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Igs0{i} = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Iorig0{i} = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Iinit0{i} = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Irgb0{i}  = zeros(outrows,outcols,3,nlist,ubits);   % channel i, all frames  
%     end
% elseif bits==16
%     for i = 1:nchannels
%         ubits = 'uint16';
%         allImages0{i}  = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Igs0{i} = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Iorig0{i} = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Iinit0{i} = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Irgb0{i}  = zeros(outrows,outcols,3,nlist,ubits);   % channel i, all frames        
%     end
% elseif bits==32
%     for i = 1:nchannels
%         ubits = 'uint32';
%         allImages0{i}  = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Igs0{i} = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Iorig0{i} = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Iinit0{i} = zeros(outrows,outcols,nlist,ubits);   % channel i, all frames
%         Irgb0{i}  = zeros(outrows,outcols,3,nlist,ubits);   % channel i, all frames     
%     end    
% end


for k = 1:nchannels   

    for i = 1:nlist
        %---------------------------------------------
        movieFrame=list(i); % TIF frame indices for temporal cropping 
        frame = i; % indices for saving data
        x0i = x0{k}(i); x1i = x1{k}(i);
        y0i = y0{k}(i); y1i = y1{k}(i);

        [Xgs0,Xorig0] = readframe(fullFileNames{k},movieFrame); %uint16 output (not rgb)    
        Xgs = Xgs0(y0i:y1i,x0i:x1i); % spatial cropping or full area 
        Xorig = Xorig0(y0i:y1i,x0i:x1i);
        Igs0{i,k} = Xgs;
        Iorig0{i,k}= Xorig;
        % define initial images for computations
        % segment based on I (gaussian smoothed membrane maker) 
        % or segment based on Iorig (original membrane marker)
        Xinit = Xgs;
        Iinit0{i,k} = Xinit; 
        
        % make boundary marker into rgb image for plotting
        % use imadjust to increase contrast
        Xrgb = repmat(imadjust(Xinit),1,1,3);
        Irgb0{i,k} = Xrgb;
        %----------------------------------
                  
        if ~submask && movieFrame==im0 && k == maskchannel 
            % if no spatial cropping/segmentation,
            % only do initial segment with first frame of membrane marker channel
            % autosegment without oBW bc no previous frame
            BW        = segment(Xinit);        
            f4 = figure(4); clf(f4);
            f4.Position(1:2) = [300,600];
            imshow(Xinit,[],'InitialMagnification',initmag);
            axis equal, hold on
            visboundaries(BW,'Color','y');
            FontSize = 18; %min(minfont,floor(min(cols,rows)/10));
            text(10,10,num2str(movieFrame),'Color','white','FontSize',FontSize)
            drawnow
            choice = timeoutdlg(delay);
            if ((choice~='n')&&(choice~='N'))&&(choice~=0)
                choice2 = secondtimeoutdlg(delay);
                BW = fixThisFrame(Xinit,BW,choice2);
            end
            title('initializing mask based on only first frame')
            hold off
            %drawnow

        elseif ~submask && movieFrame > im0 && k == maskchannel
            oBW = BW;
            BW = segment(Xinit,oBW);
            if ismember(movieFrame,frames2fix)
                f4 =figure(4); clf(f4);
                f4.Position(1:2) = [300,600];
                imshow(Xinit,[],'InitialMagnification',initmag);
                axis equal, hold on
                visboundaries(BW,'Color','y');
                FontSize = 18; %min(minfont,floor(min(cols,rows)/10));
                text(10,10,num2str(movieFrame),'Color','white','FontSize',FontSize)
                drawnow
                choice = timeoutdlg(delay);
                if ((choice~='n')&&(choice~='N'))&&(choice~=0)
                    choice2 = secondtimeoutdlg(delay);
                    BW = fixThisFrame(Xinit,BW,choice2);
                end
                title('finalizing mask')
                hold off
                %drawnow
            else
                % just visualize 
                figure(4); clf;
                imshow(Xinit,[],'InitialMagnification',initmag);
                axis equal, hold on
                visboundaries(BW,'Color','y');
                FontSize = 18; %min(minfont,floor(min(cols,rows)/10));
                text(10,10,num2str(movieFrame),'Color','white','FontSize',FontSize)
                title('finalizing mask')
                hold off
                drawnow
            end
 

        elseif submask && movieFrame == im0 && k == maskchannel
            % Xinit with segmentation from the initial mask
            BWinit = BW0(:,:,i);
            Xinit_seg = immultiply(BWinit,Xinit);
            BW = segment(Xinit_seg);
            if ismember(movieFrame,frames2fix)
                % autosegment without oBW bc no previous frame
                figure(4); clf;
                f4.Position(1:2) = [300,600];
                imshow(Xinit_seg,[],'InitialMagnification',initmag);
                axis equal, hold on
                visboundaries(BW,'Color','y');
                FontSize = 18; %min(minfont,floor(min(cols,rows)/10));
                text(10,10,num2str(movieFrame),'Color','white','FontSize',FontSize)
                drawnow
                choice = timeoutdlg(delay);
                if ((choice~='n')&&(choice~='N'))&&(choice~=0)
                    choice2 = secondtimeoutdlg(delay);
                    BW = fixThisFrame(Xinit_seg,BW,choice2);
                end
                title('finalizing mask')
                hold off
                %drawnow
            else
                % just visualize in first round
                figure(4); clf;
                imshow(Xinit_seg,[],'InitialMagnification',initmag);
                axis equal, hold on
                visboundaries(BW,'Color','y');
                FontSize = 18; %min(minfont,floor(min(cols,rows)/10));
                text(10,10,num2str(movieFrame),'Color','white','FontSize',FontSize)
                title('finalizing mask')
                hold off
                drawnow
            end

        elseif submask && movieFrame > im0 && k == maskchannel
            oBW = BW;
            BWinit = BW0(:,:,i);
            Xinit_seg = immultiply(BWinit,Xinit);
            BW  = segment(Xinit_seg,oBW);
                %BW  = segment(I1init);
            if ismember(movieFrame,frames2fix)
                f4 =figure(4); clf(f4);
                f4.Position(1:2) = [300,600];
                imshow(Xinit_seg,[],'InitialMagnification',initmag);
                axis equal, hold on
                visboundaries(BW,'Color','y');
                FontSize = 18; %min(minfont,floor(min(cols,rows)/10));
                text(10,10,num2str(movieFrame),'Color','white','FontSize',FontSize)
                drawnow
                choice = timeoutdlg(delay);
                if ((choice~='n')&&(choice~='N'))&&(choice~=0)
                    choice2 = secondtimeoutdlg(delay);
                    BW = fixThisFrame(Xinit_seg,BW,choice2);
                end
                title('finalizing mask')
                hold off
               
                %drawnow
            else
                % just visualize 
                figure(4); clf;
                imshow(Xinit_seg,[],'InitialMagnification',initmag);
                axis equal, hold on
                visboundaries(BW,'Color','y');
                FontSize = 18; %min(minfont,floor(min(cols,rows)/10));
                text(10,10,num2str(movieFrame),'Color','white','FontSize',FontSize)
                title('finalizing mask')
                hold off
                drawnow
            end
        end        

%     % add noise to the mask channel
%     if k == maskchannel || k == 3 % ActA is channel 3
%         Xorig_noise = imnoise(Xorig,"gaussian",0,0.1);
%         BW2 = imerode(BW,strel("disk",20));
%         BW3 = logical(BW-BW2);
%         Xorig = Xorig + immultiply(Xorig_noise,BW3);
%     end

  
    allImages0{i,k} = Xorig;
    if k == maskchannel
        allBW0(:,:,i) = BW;
        bmask = bwperim(BW);
        Xfused = labeloverlay(Xrgb,bmask,'colormap',[1 1 0]);
        allfused(:,:,:,i) = Xfused;
    end
    
%     figure(4);
%     FontSize = 18; %min(minfont,floor(min(cols,rows)/10));
%     text(10,10,num2str(movieFrame),'Color','white','FontSize',FontSize)
        
    if saveseg % create tiff series segmentation
        % save mask and image as tiff series if did segmentation in matlab
        if bits == 8; bwim = im2uint8(BW); 
        elseif bits == 16; bwim = im2uint16(BW);
        elseif bits == 32; bwim = im2uint32(BW);
        end
        if movieFrame == im0 
            imwrite(Xorig,strcat(outputfolder,fnames{k},'_final.tif'));
            if k == maskchannel
                    imwrite(bwim,strcat(outputfolder,maskfname,'_final.tif'));
                    imwrite(Xfused,strcat(outputfolder,maskfname,'_fused.tif'));
            end
        else
            imwrite(Xorig,strcat(outputfolder,fnames{k},'_final.tif'),...
                'WriteMode','append');
            if k == maskchannel
                imwrite(bwim,strcat(outputfolder,maskfname,'_final.tif'),...
                    'WriteMode','append');
                imwrite(Xfused,strcat(outputfolder,maskfname,'_fused.tif'),...
                    'WriteMode','append');
            end
        end
  
    end % end of savefigs
    
    end % end of nlist

end % end of nchannels


if savefigs
    % save mfile if saving figs 
    currentfile = strcat(mfilename,'.m');
    destinationfile = strcat(outputfolder,mfilename,'_',tstamp1,'.m');
    copyfile(currentfile,destinationfile);   
end

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

function nBW = fixThisFrame(I,BW,fixFrameCommand)
    if nargin<3
        fixFrameCommand = 'circle';
    end
    f4 = figure(4); clf(f4);
    imshow(I,[],'InitialMagnification',600);
    axis equal, hold on
    if strcmp(fixFrameCommand,'rectangle')
        text(48,48,'Use the mouse to draw a rectangle around the desired cell',...
            'Color','y','FontSize',16)
        BW = bwareaopen(BW,200);
        visboundaries(BW,'Color','y');
        drawnow
        roi = drawrectangle;
        wait(roi);
    elseif strcmp(fixFrameCommand,'polygon')
        text(48,48,'Use the mouse to draw a polygon around the desired cell',...
            'Color','y','FontSize',16)
        BW = bwareaopen(BW,200);
        visboundaries(BW,'Color','y');
        drawnow
        roi = drawpolygon;
        wait(roi);
    else %strcmp(fixFrameCommand,'circle')
        text(48,48,'Use the mouse to draw a circle around the desired cell',...
            'Color','y','FontSize',16)
        BW = bwareaopen(BW,200);
        visboundaries(BW,'Color','y');
        drawnow
        roi = drawcircle;
        wait(roi);
    end

    msk = createMask(roi);
    nBW = activecontour(immultiply(I,msk),msk,200,'Chan-Vese','SmoothFactor',1,'ContractionBias',0);
    nBW = bwareaopen(nBW,100);
    f4 = figure(4); clf(f4);
    imshow(I,[],'InitialMagnification',600);
    axis equal, hold on
    visboundaries(nBW,'Color','r');
    drawnow
    pause
    hold off
end

function rectpos = rectcropThisFrame(I)
    figure(1); clf;
    imshow(I,[],'InitialMagnification',400);
    axis equal, hold on
    text(48,48,'Use the mouse to draw a rectangle around the desired cell',...
        'Color','y','FontSize',16)
    roi = drawrectangle;
    wait(roi);
    rectpos = roi.Position;
    hold off
end

function BW = segment(I0,oBW,varargin)

    prop = 0.05; %0.05
    if isempty(varargin)
        minArea = 500; %100
        maxDistance = 30; %30;
    else
        minArea = varargin{1};
        maxDistance = varargin{2};
    end
    I = im2uint16(I0);
    T = graythresh(I);
    if nargin==2
        % dBW = distance to nearest "logical 1" pixel;
        % dBW is high for background far from the target in oBW
        dBW  = bwdist(oBW); 
        %figure(5); imshow(dBW,[]);
        % Itmp1 is ON for target and decreases as pixels get farther away
        % from target
        Itmp1 = 1./(1+prop*dBW); 
        %figure(6); imshow(Itmp1,[]);
        % make a cutoff for Itmp to be 0 beyond a certain distance
        Itmp2 = immultiply(dBW<maxDistance,Itmp1);
        %figure(7); imshow(Itmp2,[]);
        Iplus    = uint16(immultiply(single(I),Itmp2));
    else
        Iplus = I;
    end
    %BW = imbinarize(Iplus);
    %BW = imbinarize(Iplus,1.75*T); % for L03
    BW = imbinarize(Iplus,.75*T); % original
    % fill holes
    BW = imfill(BW,'holes');
    BW = bwareaopen(BW,minArea); % remove objects containing less than minArea pixels
    % 10 iterations
    BW = activecontour(Iplus,BW,10,'Chan-Vese','SmoothFactor',1,'ContractionBias',0);
    % fills gaps
    BW = imclose(BW,strel('disk',3,0));
    BW = imfill(BW,'holes');
    BW = bwareaopen(BW,minArea);
    %figure(5); imshow(BW,[])
end

% function BW = segment(I,oBW,varargin)
%     if isempty(varargin)
%         maxDistance = 5; %30;
%         minArea = 200; %100
%         prop = 0.05; %0.05, later 0.5
% 
%     else
%         minArea = varargin{1};
%         maxDistance = varargin{2};
%     end
%     T = graythresh(I);
%     if nargin==2
%         dBW  = bwdist(oBW); % distance to nearest "1/true" pixel; is 0 for mask area
%         Itmp = 1./(1+prop*dBW); % semi-mask 
%         Itmp = immultiply(dBW<maxDistance,Itmp);
%         I    = uint16(immultiply(single(I),Itmp));
%     end
%     BW = imbinarize(I,.75*T);
%     BW = imfill(BW,'holes');
%     BW = bwareaopen(BW,minArea); % remove objects containing less than minArea pixels
%     BW = activecontour(I,BW,10,'Chan-Vese','SmoothFactor',1,'ContractionBias',0);
%     BW = imclose(BW,strel('disk',3,0));
%     BW = imfill(BW,'holes');
%     BW = bwareaopen(BW,minArea);
% end

function choice = timeoutdlg(delay)

    prompt='Manual over ride [y/n]?';
    f1 = findall(0, 'Type', 'figure');
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

function choice = secondtimeoutdlg(delay)

    prompt='draw shape [0 = rectangle, 1 = circle, 2 = polygon]';
    f1 = findall(0, 'Type', 'figure');
    %delay;
    t = timer('TimerFcn', {@closeit f1}, 'StartDelay', delay);
    start(t);
    % Call the dialog
    retvals = inputdlg(prompt,'',1,{'1'});
    if numel(retvals) == 0
          choice = '1';
    else
          if strcmp(retvals{1},'0') 
              choice = 'rectangle';

          elseif strcmp(retvals{1},'2') 
              choice = 'polygon';

          else
              choice = 'circle';

          end

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


