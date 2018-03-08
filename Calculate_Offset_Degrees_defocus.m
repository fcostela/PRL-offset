cd('/Users/FranciscoCostela/Desktop/magnification')

% clear
t0=clock;
nvPath='/Users/FranciscoCostela/Desktop/GazeContigentDisplay/matlab_undersample';
nvPath2='/Users/FranciscoCostela/Desktop/magnification/Video clip eyetracking data';

etPath='/Users/FranciscoCostela/Desktop/PRL study/CVL';

movPath='/Users/FranciscoCostela/Desktop/magnification/ClipsForNorming';
croiPath = '/Users/FranciscoCostela/Desktop/magnification/crois';%croi_subject';
EditPath = '/Users/FranciscoCostela/Desktop/magnification/ClipsEdited';
hemioPath = '/Users/FranciscoCostela/Desktop/hemianopes';
defocusPath = '/Users/FranciscoCostela/Desktop/PRL study/defocus';

clear etdata;
load 'video number lookup.mat'

screenDist=100;       % distance in cm
screenSize=[60 34];   % [horz vert] screen dimensions in cm
screenRes=[2560 1440];%[1920 1080] for projector %screen dims in pixels

nv = 1;
defocus = 1;
hemios = 0;
hemiosnew = 1;
dofigures = 1;
same_videos = 0;
%First find which clips the CVL subjects saw...

% subjects taken from pelishare/projects/ CVL IA/ Data / raw / cvl
%oldcvlcode = {'1507ye','2888ea','2765ma', '2678nr', '2498ar','2449ai','2165rs','2113ad','1995dr','1730yy','1508in','0785dn','0761rg','0212dy'};
%newcvlcode = {'1208dd', '1732ri', '1783az','2338ms','2339be','2484ao','3024mn'};

cvlcode = {'1507ye','2888ea','2765ma', '2678nr', '2498ar','2449ai','2165rs','2113ad','1995dr','1730yy','1508in','0785dn','0761rg','0212dy','1208dd', '1732ri', '2338ms','2484ao','3024mn', '1783rv'}; %2339be
oldYesNo = [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0];

%remember that 0 is to the right, 90 is north
theta = [96.7436 nan 19.2763 nan nan 108.2694 5.9781 nan nan nan 175.1178 nan 99.4825 73.905 86.0094 79.6773 75.9121 18.6412 nan 87.8261  ];
eccentricity = [0.7823 nan 0.9576 nan nan 0.6288 0.7193 nan nan nan 0.9376 nan 0.974 0.9839 0.8485 0.8587 0.8279 0.878 nan 0.7525  ];

% cvl_age = [71, 74, 32, 39, 57, 80, 85, 72, 87, 29, 45, 67, 48, 63];
% figure;plot(30:10:90,hist(cvl_age,30:10:90));

nvcode = {'2005dh','2152so','2764hn', '0811ne', '2237ns','2157ss','2158yn','2278te','0972rs','2753nt','2037sn','2086ty','2153ng', '2072ds', '2176as','1613ey','2130oa','2030nr','1848ne', '2703le'}; %'2081lo'
nv_age = [71, 73, 27, 46, 58, 80, 85, 72, 83, 27, 67, 48, 63];
oldYesNoNV = ones(1,21); %[ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];


defocuscode = {'66ss','1061sm','1992au', '2499sm', '2654fl','2664nn','2723rr','2724xz','2726cn','2727wi','2728lr','2729ao','2742yn','2743ne','2751wd'}; %'2081lo'


%nvcode = {'2764hn','2086ty'};
%nvcode = {'2764hn', '0811ne', '2158yn','2278te','0972rs','2086ty'};

%load 'video number lookup.mat'
%This line below are actually hemianopes: first 7 are left, last 4 are right
%newer hemianope starts with 2302kr
hemiocode = { '1897db', '2604ns','2220er','2750hy','2769tn', '2821sy','2801tn', '1723pn', '2820sl','2851ye','2500ns', '2396ny',  '2918ty' , '2813kn', '2898ss', '2302kr', '2549nt', '2644en', '2830', '2876as' };
hemiolabel = {'L', 'L', 'L' , 'L', 'L', 'L', 'L', 'R', 'R', 'R', 'R' , 'R' , 'R', 'L', 'R',   'R', 'L', 'L', 'L', 'L'};

if hemios
    cvlcode = hemiocode;
else if nv
        if defocus
            cvlcode = defocuscode;
        else
            cvlcode = nvcode;
        end
    oldYesNo = oldYesNoNV;
    end
end



 cvl_hists = {};
 nv_hists = {};
 offset_mean = {};

% We iterate through all the subjects (there are also the same number in the NV group)
for index=1:length(cvlcode)
    
    % files = dir([etPath '/freeNorm ' cvlcode{index} '*.mat']);
    index
    disp(cvlcode{index});
    % We may want to get all the videos ~40 from the NV subjects. In that
    % case    
   
    % Subjects from old data watched with this resolution (IMAC screen)
    screenRes=[2560 1440];
    % New subjects watched with the projector configuration    
    
    if ~nv
        [videos matfile ] = extractVideosSeen(1,'2338ms');
        if ~oldYesNo(index)
            disp('It''s a new subject');
            screenRes=[1920 1080];
        end
        %load cvl_hist
        
        if hemios
            [myvideos matfile ] = extractVideosSeen(2,cvlcode{index});
            numVideos = length(myvideos);
        else
            [myvideos matfile ] = extractVideosSeen(1,cvlcode{index});
            numVideos = length(videos);
        end
        
    else if same_videos
            %load nv_hist_same
            [myvideos matfile ] = extractVideosSeen(2,nvcode{index});
            numVideos = length(videos);
        else
            %load nv_hist_all
            if defocus
                [myvideos matfile ] = extractVideosSeen(3,defocuscode{index});
            else
                [myvideos matfile ] = extractVideosSeen(0,nvcode{index});
            end
                numVideos = min(40,length(myvideos));
        end
    end
    
    hist_offset_mean = {};
    for video_index = 1:numVideos
        
        if same_videos || ~nv
            myvideo_index = find(myvideos == videos(video_index));
            if isempty(myvideo_index)
                continue;
            end
        end
        
        % This line also changes depending on plotting nv or CVL
        % subjects
        if ~nv
            if hemios
                 % load([hemioPath '/no_guide/' char(cvlcode{index}) '/' char(strrep(matfile(myvideo_index),' ', '_'))]);
                 % load ( [nvPath '/' char(matfile(myvideo_index(1)))]); 
                  load ( [nvPath '/' char(matfile(video_index(1)))]);
            else
%                 if exist([nvPath filesep char(cvlcode{index}) '/' char(matfile(myvideo_index))])>0
%                     load([nvPath filesep char(cvlcode{index}) '/' char(matfile(myvideo_index))]);
%                 else
%                     load([etPath filesep char(cvlcode{index}) '/' char(matfile(myvideo_index))]);
%                 end                
                 load( [ etPath '/' char(cvlcode{index}) '/' char(matfile(myvideo_index))]);
            end
        else
            % Uncomment this line if we want to match the same videos
            % watched by CVL (but not all the videos would be in that pool)
            if same_videos
                load ( [nvPath '/' char(matfile(myvideo_index))]);
            else
                if exist([nvPath filesep char(matfile(video_index))])>0
                    load([nvPath filesep char(matfile(video_index))]);
                else if defocus
                        load([defocusPath filesep char(cvlcode{index}) filesep 'TV_Watching Data/' char(matfile(video_index))]);
                        
                    else
                        load([nvPath2 filesep char(matfile(video_index))]);
                    end
                end
                
                %load ( [nvPath '/' char(matfile(video_index))]);
            end
        end
        
        % get movie file name
        [pathstr name ext] = fileparts(movieFileName);
        
        %these two lines to read uncompressed video version
        name = strrep(name, '_c 2','');
        name = strrep(name, '_c','');
        
        temp.x = eyetrackRecord.x;
        temp.y = eyetrackRecord.y;
        temp.t = eyetrackRecord.t;
        temp.missing = eyetrackRecord.missing;
        etdata = temp;
        
        disp(['Clip number:' num2str(video_index) ' index: ' num2str(videoNumber) 'File:' name] );
        % get ET data for the current subject that viewed this clip
        
        subind=[1:length(etdata)];
        %frametime=1000*linspace(0,30,length(movmat));
        frametime=1000*linspace(0,30,720);
        
        etSubjAll=cell(length(frametime)-1,1);
        sdims=[2560 1440]; % hor x vert
        % each subject
        % sort coords and time data
        if sum(etdata.missing) == size(etdata.x,2)
            disp('Video without gaze data');
            continue;
        end
        
        etxy=[etdata.x' etdata.y'];
        ettime=etdata.t';
        % trim coords that are outside of screen
        indx=etxy(:,1)>=0 & etxy(:,1)<=sdims(1);
        indy=etxy(:,2)>=0 & etxy(:,2)<=sdims(2);
        etxy=etxy((indx+indy)>=2,:);
        if isempty(etxy)
            disp('Video without gaze data');
            continue;
        end
        % trim time indices
        ettime=ettime((indx+indy)>=2,:);
        % get time elapsed since beginning of clip
        % for each recorded eye position
        eltime=ettime-ettime(1);
        eltime=eltime(1:end-1);
        % detect and replace saccades
        [etxy_new sacind]=detectSaccades(etxy,eltime);
        for i=1:length(frametime)-1 % each frame
            etind=find(eltime>=frametime(i) & eltime<frametime(i+1));
            etSubjAll{i}=vertcat(etSubjAll{i},[etxy_new(etind,1) etxy_new(etind,2)]);
        end
        
        % If this is a CVL subject from old data, load directly the croi we saved for his
        % movies
        if ~nv 
            disp('loading CROI...');
%             load([croiPath '/croi' char(cvlcode{index}) '-' num2str(videoNumber) '.mat'],'croi');
            load([croiPath '/' num2str(videoNumber) '.mat'],'croi');
        else %For nv subjects, we have to recalculate the croi without their gaze
            
            if ~defocus
                load([croiPath '/croi' char(cvlcode{index}) '-' num2str(videoNumber)],'croi');
            else
                load([croiPath '/' num2str(videoNumber) '.mat'],'croi');
            end
            
            % Only uncomment the code below if no crois have been
            % calculated for NV subjects where the CROI does not include
            % their own data
            
          %  load ([croiPath '/croi' char(cvlcode{index}) '-' num2str(videoNumber)]);
%             subs4vid = find(videoNumbers == videoNumber);            
%             
%             skipped = 0;
%             for i = 1:length(subs4vid)
%                 
%                 load([nvPath2 filesep eyetrackFiles{subs4vid(i)}]);
%                 % we exclude nv subjects if they watched the clip. we
%                 % do not want their data in the roi
%                 if isempty(strfind(eyetrackFiles{subs4vid(i)},char(nvcode(index))))
%                     %load([nvpath filesep eyetrackFiles{subs4vid(i)}]);
%                     temp.x = eyetrackRecord.x;
%                     temp.y = eyetrackRecord.y;
%                     temp.t = eyetrackRecord.t;
%                     temp.missing = eyetrackRecord.missing;
%                     etdata(i) = temp;
%                 else
%                     %  disp('skipping one subject');
%                     skipped = i;
%                 end
%             end
%             if skipped>0 && skipped<size(etdata,2)
%                 etdata(skipped) = [];
%             end
%             
%             
%             %% get all eye positions from all subjects for each frame
%             subind=[1:length(etdata)];
%             %frametime=1000*linspace(0,30,length(movmat));
%             frametime=1000*linspace(0,30,720);
%             
%             etAll=cell(length(frametime)-1,1);
%             sdims=[2560 1440]; % hor x vert
%             for sub=1:length(etdata) % each subject
%                 % sort coords and time data
%                 etxy=[etdata(sub).x' etdata(sub).y'];
%                 ettime=etdata(sub).t';
%                 % trim coords that are outside of screen
%                 indx=etxy(:,1)>=0 & etxy(:,1)<=sdims(1);
%                 indy=etxy(:,2)>=0 & etxy(:,2)<=sdims(2);
%                 etxy=etxy((indx+indy)>=2,:);
%                 % trim time indices
%                 ettime=ettime((indx+indy)>=2,:);
%                 % get time elapsed since beginning of clip
%                 % for each recorded eye position
%                 eltime=ettime-ettime(1);
%                 eltime=eltime(1:end-1);
%                 % detect and replace saccades
%                 [etxy_new sacind]=detectSaccades(etxy,eltime);
%                 for i=1:length(frametime)-1 % each frame
%                     etind=find(eltime>=frametime(i) & eltime<frametime(i+1));
%                     etAll{i}=vertcat(etAll{i},[etxy_new(etind,1) etxy_new(etind,2)]);
%                 end
%             end
%             
%             %% integrate the kdes over each frame
%             disp('creating ROI');
%             destdims=sdims;%(end:-1:1); % [vert horz]
%             sigma=400;
%             mag=2;
%             croi=zeros(length(etAll),2);
%             
%             for i=1:length(etAll)
%                 [croi(i,:) val]=integROI(etAll{i},destdims,sigma,mag);
%             end
%             
%             save ([croiPath '/croi' char(cvlcode{index}) '-' num2str(videoNumber)], 'croi');
        end
        
        %         disp('Time to compare with COI!');
        
        croiDeg = pix2deg(croi,screenDist,screenSize,screenRes);
        for frame=1:size(croi,1)
            cvl_croi(frame,1) = nanmean(etSubjAll{frame}(:,1));
            cvl_croi(frame,2) = nanmean(etSubjAll{frame}(:,2));
            
            cvl_croiDeg = pix2deg(cvl_croi,screenDist,screenSize,screenRes);
            
            %My definition of offset is the COI from subject - democratic
            %COI
            offset(frame,:) =cvl_croiDeg(frame,1:2) - croiDeg(frame,1:2);
            
        end
        
        %Using a 2d Histogram to compare the location of the PRL relative
        %to the democratic COI
        %         if dofigures
        %             if same_videos || ~nv
        %                 subplot(4,6,video_index);
        %             else
        %                 subplot(4,6,video_index);
        %             end
        %         end
        
        % this was to control gazes out of the screen
        %offset = offset(abs(offset(:,1))<sdims(1) & abs(offset(:,2))<sdims(2),:);
        
        %for pixels
        %         vXEdge = linspace(-1400, 1400, 20); % min, max, resolution respectively
        %         vYEdge = linspace(-1000, 1000, 20);
        % for degrees
        vXEdge = linspace(-19, 19, 20); 
        vYEdge = linspace(-19, 19, 20); 
        mHist2d = hist2d([offset(:,2) offset(:,1)], vYEdge, vXEdge);
        
        %         if dofigures
        %             nXBins = length(vXEdge);
        %             nYBins = length(vYEdge);
        %             vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
        %             vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
        %             pcolor(vXLabel, vYLabel,mHist2d); %colorbar
        %         end
        
        hist_offset_mean{video_index} = mHist2d;
        offset_mean{index}(video_index,1) = nanmean(offset(:,1));
        offset_mean{index}(video_index,2) = nanmean(offset(:,2));
        list_videos{index} = myvideos;
        myhist{index} = mHist2d;
        
        
        %     if dofigures
        %         if ~nv
        %             figname = [char(cvlcode(index)) '.png'];
        %         else if same_videos
        %                 figname = [char(nvcode(index)) '.png'];
        %             else
        %                 figname = [char(nvcode(index)) 'All.png'];
        %             end
        %         end
        %         print('-dpng',figname);
        %         close all
        %     end
        all_hists{index}{video_index} = mHist2d;
        % Very useful way to get the average from all the matrixes in the
        % struct
        
        Mcell = arrayfun(@(x) hist_offset_mean{x}, 1:length(hist_offset_mean), 'uni', 0);
        if ~nv
            cvl_hists{index} = nanmean( reshape(cell2mat(Mcell), 19, 19, []), 3 );
        else
            nv_hists{index} = nanmean( reshape(cell2mat(Mcell), 19, 19, []), 3 );
        end
        
        if nv
            if same_videos
                save('nv_hist_sameDegrees.mat', 'nv_hists', 'offset_mean');
            else
                if ~defocus
                    save('nv_hist_allDegrees.mat', 'nv_hists', 'offset_mean');
                else
                     save('defocus_hist_allDegrees.mat', 'nv_hists', 'offset_mean', 'all_hists');
                end
            end
        else
            if hemios
                save('hemios_hist.mat', 'cvl_hists', 'offset_mean' , 'list_videos');
            else
                save('cvl_histDegrees.mat', 'cvl_hists', 'offset_mean');
            end
        end
    end
end

if dofigures
    for myindex=1:length(cvlcode)
        
        myindex        
        if defocus
            subplot(3,5,myindex);
        else
            subplot(4,5,myindex);
        end
        nXBins = length(vXEdge);
        nYBins = length(vYEdge);
        vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
        vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
        if ~nv
            pcolor(vXLabel, vYLabel,cvl_hists{myindex});
            if hemios
                title([ char(hemiocode(myindex)) ' ' hemiolabel{myindex}],'FontName','Arial','FontSize',20);
            else
               % title(char(cvlcode(myindex)),'FontName','Arial','FontSize',20);
                 title(['P' num2str(myindex)],'FontName','Arial','FontSize',20);
            end
        else
             pcolor(vXLabel, vYLabel,nv_hists{myindex});
            hold on
            plot(0,0,'MarkerSize',8,'Marker','+','LineWidth',2,'Color',[1 1 1]);
            
           % pcolor(vXLabel, vYLabel,nv_hists{myindex});            
          %  title(char(nvcode(myindex)),'FontName','Arial','FontSize',20);
          title(['Defocus ' num2str(myindex)],'FontName','Arial','FontSize',20);
        end
        set(gca, 'FontName','Arial','FontSize',16);
        
    end
    if same_videos || ~nv
        print('-dpng','global.png');
    else
        print('-dpng','globalAll.png');
    end    
    
end
tend=etime(clock,t0)/60;
% 
% 
distance=[];
distance2 = [];
distance_cvl= [];
distance_left = [];
distance_right = [];
% Code to analyze
load ('/Users/FranciscoCostela/Desktop/magnification/nv_hist_allDegreesGOOD.mat');
 all_means = offset_mean ;
 NV_means = offset_mean;
 nv_xy = [];
 
for i=1:length(all_means)
   
    nv_xy(i,1) = nanmean(all_means{i}(:,1));
    nv_xy(i,2) = nanmean(all_means{i}(:,2));
    
    for j=1:length(all_means{i})
        angle(i,j) = radtodeg(atan2(all_means{i}(j,2) , all_means{i}(j,1)));
        distance(i,j) = sqrt(all_means{i}(j,1) ^2 + all_means{i}(j,2)^2 ); 
        offset_x (i,j) = offset_mean{i}(j,1);
        offset_y (i,j) = offset_mean{i}(j,2);
        distance_x(i,j) = abs(all_means{i}(j,1));
        distance_left(i,j) = -all_means{i}(j,1);
        distance_right(i,j) = all_means{i}(j,1);
   %scatter( all_means{i}(j,1), all_means{i}(j,2),5,[randi(i)/13 0.5 randi(i)/13],'o','filled');  
  hold on
    end
end

offset_x(~offset_x) = nan;
offset_x(:,41:46) = [];
offset_y(~offset_y) = nan;
offset_y(:,41:46) = [];
offset_x(20,:) = [];
offset_y(20,:) = [];
u = ksdensity(offset_x(:),offset_x(:),'function','cdf');
v = ksdensity(offset_y(:),offset_y(:),'function','cdf');
%[Rho,nu, paramci] = copulafit('t',[u v],'Method','ApproximateML');
[paramhat, paramci] = copulafit('Clayton', [u v]);

x1 = ksdensity(x,u1,'function','icdf');
y1 = ksdensity(y,v1,'function','icdf');


for i=1:20
[hx120(i) p(i)] = ttest(offset_x(i,:), 0, 0.01,'both');
[hy120(i) p(i)] = ttest(offset_y(i,:), 0, 0.01,'both');
end


hx120 | hy120
%It gives the normal distribution mu and sigma
pd = fitdist(distance(:),'Normal')
threshold = pd.sigma * 1.96 + pd.mu;
 
distance(~distance) = nan;
mean_angle = nanmean(angle,2);
mean_distance = nanmean(distance,2);

PRLYesNo = isnan(theta);
for PRL_index = 1:length(all_means)
    if ~PRLYesNo(PRL_index)
        disp(['Subject ' cvlcode{PRL_index}]);
        disp(['Angle:' num2str(mean_angle(PRL_index)) ' - ' num2str(theta(PRL_index))]);
        disp(['Ecc:' num2str(mean_distance(PRL_index)) ' - ' num2str(eccentricity(PRL_index))]);
        disp(['Diff:' num2str(180 - abs(abs(mean_angle(PRL_index) - theta(PRL_index)) - 180))]);
    end
    
end

% % for i=1:length(all_means20)
% %    
% %     for j=1:length(all_means20{i})
% %       
% %         if all_means20{i}(j,1)~= 0
% %         distance2(i,j) = sqrt(all_means20{i}(j,1) ^2 + all_means20{i}(j,2)^2 );
% %         end
% %    %scatter( all_means{i}(j,1), all_means{i}(j,2),5,[randi(i)/13 0.5 randi(i)/13],'o','filled');
% %     
% %   hold on
% %     end
% % end
% left_index = 1;
% right_index = 1;
% distance_cvl_left_left= [];
% distance_cvl_left_right= [];
% distance_cvl_right_left= [];
% distance_cvl_right_right= [];

load ('/Users/FranciscoCostela/Desktop/magnification/defocus_hist_allDegrees.mat');
 all_means = offset_mean ;
 defocus_means = offset_mean;
 % There are 5 defocus level. Time to create 5 categories and order videos
 % and offsets according to these levels. Info has been stored in
 % defocusdata.mat
 
 load ('/Users/FranciscoCostela/Desktop/PRL study/defocus/defocusdata.mat');
   
        
 distanceDefocus = [];
 distanceDefocus1 = zeros(length(all_means),4);
 distanceDefocus2 = zeros(length(all_means),4);
 distanceDefocus3 = zeros(length(all_means),4);
 distanceDefocus4 = zeros(length(all_means),4);
 distanceDefocus5 = zeros(length(all_means),4);
 defocus_xy = [];
 
 hists_levels1 = {};
 hists_levels2 = {};
 hists_levels3 = {};
 hists_levels4 = {};
 hists_levels5 = {};
 histIndex = ones(1,5);
 
 
    
for i=1:length(all_hists)
   i1 = 1;
 i2 = 1;
 i3 = 1;
 i4 = 1;
 i5 = 1;
   
    defocus_xy(i,1) = nanmean(all_means{i}(:,1));
    defocus_xy(i,2) = nanmean(all_means{i}(:,2));
    
    for j=1:min(20,length(all_hists{i}))
        
        switch(defocuslevel(i,j))
            
            case 1
                distanceDefocus1(i,i1) = sqrt(all_means{i}(j,1) ^2 + all_means{i}(j,2)^2 );
                hists_levels1{histIndex(1)} = all_hists{i}{j};
                i1 = i1+1;
                histIndex(1) = histIndex(1) +1;
            case 2
                distanceDefocus2(i,i2) = sqrt(all_means{i}(j,1) ^2 + all_means{i}(j,2)^2 );
                hists_levels2{histIndex(2)} = all_hists{i}{j};
                i2 = i2+1;
                histIndex(2) = histIndex(2) +1;
                
            case 3
                distanceDefocus3(i,i3) = sqrt(all_means{i}(j,1) ^2 + all_means{i}(j,2)^2 );
                hists_levels3{histIndex(3)} = all_hists{i}{j};
                i3 = i3+1;
                histIndex(3) = histIndex(3) +1;
                
            case 4
                distanceDefocus4(i,i4) = sqrt(all_means{i}(j,1) ^2 + all_means{i}(j,2)^2 );
                hists_levels4{histIndex(4)} = all_hists{i}{j};
                i4 = i4+1;
                histIndex(4) = histIndex(4) +1;
                
            case 5
                distanceDefocus5(i,i5) = sqrt(all_means{i}(j,1) ^2 + all_means{i}(j,2)^2 );
                hists_levels5{histIndex(5)} = all_hists{i}{j};
                i5 = i5+1;
                histIndex(5) = histIndex(5) +1;
                
        end
         offset_x (i,j) = offset_mean{i}(j,1);
        offset_y (i,j) = offset_mean{i}(j,2);
        angleDefocus(i,j) = radtodeg(atan2(all_means{i}(j,2) , all_means{i}(j,1)));
        distanceDefocus(i,j) = sqrt(all_means{i}(j,1) ^2 + all_means{i}(j,2)^2 );  
        distanceDefocus_x(i,j) = abs(all_means{i}(j,1));
        distanceDefocus_left(i,j) = -all_means{i}(j,1);
        distanceDefocus_right(i,j) = all_means{i}(j,1);
   %scatter( all_means{i}(j,1), all_means{i}(j,2),5,[randi(i)/13 0.5 randi(i)/13],'o','filled');
    
  hold on
    end
end

nonormal = [];
nonormal = cat(2, distanceDefocus2, distanceDefocus3, distanceDefocus4, distanceDefocus5);
nonormal(16:20,:) = [];
defocus_averages = nanmean(nonormal,2);

for i=1:20
[hx120(i) p(i)] = ttest(offset_x(i,:), 0, 0.01/20);
[hy120(i) p(i)] = ttest(offset_y(i,:), 0, 0.01/20);
end

hx120 | hy120

defocus_hists = {};
for levels = 1:5
    subplot(2,3,levels);
    switch levels
        case 1 
    Mcell = arrayfun(@(levels) hists_levels1{levels}, 1:59, 'uni', 0);
        case 2
         Mcell = arrayfun(@(levels) hists_levels2{levels}, 1:59, 'uni', 0);
        case 3
         Mcell = arrayfun(@(levels) hists_levels3{levels}, 1:59, 'uni', 0);
        case 4
         Mcell = arrayfun(@(levels) hists_levels4{levels}, 1:59, 'uni', 0);
        case 5
       Mcell = arrayfun(@(levels) hists_levels5{levels}, 1:59, 'uni', 0);      
    end
    defocus_hists{levels} = nanmean( reshape(cell2mat(Mcell), 19, 19, []), 3 );
    nXBins = length(vXEdge);
    nYBins = length(vYEdge);
    vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
    vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
    pcolor(vXLabel, vYLabel,defocus_hists{levels});
    hold on
            plot(0,0,'MarkerSize',8,'Marker','+','LineWidth',2,'Color',[1 1 1]);
            
            
    title(['Level ' num2str(levels)],'FontName','Arial','FontSize',20);
    
end

       
distanceDefocus(~distanceDefocus) = nan;

distanceDefocus1(~distanceDefocus1) = nan;
distanceDefocus2(~distanceDefocus2) = nan;
distanceDefocus3(~distanceDefocus3) = nan;
distanceDefocus4(~distanceDefocus4) = nan;
distanceDefocus5(~distanceDefocus5) = nan;

mean_distance_defocus1 = nanmean(distanceDefocus1,2);
mean_mean_distance_defocus1 = nanmean(mean_distance_defocus1);
std_mean_distance_defocus1 = nanstd(mean_distance_defocus1)/sqrt(length(defocuscode));

mean_distance_defocus2 = nanmean(distanceDefocus2,2);
mean_mean_distance_defocus2 = nanmean(mean_distance_defocus2);
std_mean_distance_defocus2 = nanstd(mean_distance_defocus2)/sqrt(length(defocuscode));

mean_distance_defocus3 = nanmean(distanceDefocus3,2);
mean_mean_distance_defocus3 =  nanmean(mean_distance_defocus3);
std_mean_distance_defocus3 =  nanstd(mean_distance_defocus3)/sqrt(length(defocuscode));

mean_distance_defocus4 =  nanmean(distanceDefocus4,2);
mean_mean_distance_defocus4 =  nanmean(mean_distance_defocus4);
std_mean_distance_defocus4 =  nanstd(mean_distance_defocus4)/sqrt(length(defocuscode));

mean_distance_defocus5 =  nanmean(distanceDefocus5,2);
mean_mean_distance_defocus5 =  nanmean(mean_distance_defocus5);
std_mean_distance_defocus5 =  nanstd(mean_distance_defocus5)/sqrt(length(defocuscode));

load ('/Users/FranciscoCostela/Desktop/magnification/cvl_histDegrees.mat');

cvl_means = offset_mean;
for i=1:length(offset_mean)
    
    cvl_xy(i,1) = nanmean(offset_mean{i}(:,1));
    cvl_xy(i,2) = nanmean(offset_mean{i}(:,2));
    
    if strcmp(hemiolabel{i},'L')
        colortype = [1 0 0] ;
    else
        colortype = [0 0 1];
    end
    
    if i>13
        colortype = [0 1 0] ;
    end
   
    for j=1:length(offset_mean{i})
        
        distance_cvl(i,j) = sqrt(offset_mean{i}(j,1)^2 + offset_mean{i}(j,2)^2 );
        distance_cvl_x(i,j) = abs(offset_mean{i}(j,1));   
        offset_cvl_x (i,j) = offset_mean{i}(j,1);
        offset_cvl_y (i,j) = offset_mean{i}(j,2);
%         if strcmp(hemiolabel{i},'L')
%             distance_cvl_left_left(left_index,j) = -offset_mean{i}(j,1);
%             distance_cvl_right_left(left_index,j) = offset_mean{i}(j,1);
%         else
%             distance_cvl_right_right(right_index,j) = offset_mean{i}(j,1);
%             distance_cvl_left_right(right_index,j) = -offset_mean{i}(j,1);
%             
%         end
%         %scatter( all_means{i}(j,1), all_means{i}(j,2),5,[randi(i)/15 0.5 randi(i)/15],'o','filled');
%         scatter( offset_mean{i}(j,1), offset_mean{i}(j,2), 20,[colortype(1) colortype(2) colortype(3)],'o');
%         
%          hold on
%     end
%     if strcmp(hemiolabel{i},'L')
%         left_index = left_index +1;
%     else
%         right_index = right_index +1;
     end
    meandistance(i,1) = nanmean(offset_mean{i}(:,1));
    meandistance(i,2) = nanmean(offset_mean{i}(:,2));
    scatter( meandistance(i,1), meandistance(i,2), 100,[colortype(1) colortype(2) colortype(3)],'o','filled');
        
end
 xlabel('Horizontal gaze offset (pixel)','FontName', 'Arial', 'FontSize', 18);
 ylabel('Vertical gaze offset (pixel)','FontName', 'Arial', 'FontSize', 18);

line([-400 400], [0 0], 'Color', 'k');
line([0 0], [-200 200], 'Color', 'k');
set(gca,'xlim',[-400, 400]);
set(gca,'ylim',[-200, 200], 'FontName', 'Arial', 'FontSize', 18);
xlabel('Horizontal gaze offset (pixel)');
ylabel('Vertical gaze offset (pixel)','FontName', 'Arial', 'FontSize', 18);
% 
% 
% % 
% mean_distance = mean(distance_x,2);
 for s=1:length(nvcode)
    std_distance(s) = nanstd(distance(s,:)); 
 end
% distance_left(~distance_left) = nan;
mean_mean_distance = nanmean(mean_distance);
% mean_distance_left = nanmean(distance_left,2);
% mean_mean_distance_left = nanmean(mean_distance_left);
% std_distance_left = nanstd(distance_left)/sqrt(length(nvcode));
% mean_std_distance_left = mean(std_distance_left);
% mean_distance_right = mean(distance_right,2);
% std_distance_right = std(distance_right);
% mean_std_distance_right = mean(std_distance_right);
std_mean_distance = std(mean_distance)/sqrt(length(nvcode));
% %distance2(~distance2) = nan;
% % mean_distance2 = nanmean(distance2,2);
% % for s=1:length(nvcode)
% %    std_distance2(s) = nanstd(distance2(s,:)); 
% % end
% % mean_mean_distance2 = nanmean(mean_distance2);
% % std_mean_distance2 = nanstd(mean_distance2)/sqrt(length(nvcode));
% 
distance_cvl(~distance_cvl) = nan;
mean_distance_cvl = mean(distance_cvl,2);
 for s=1:length(cvlcode)
    std_distance_cvl(s) = std(distance_cvl(s,:)); 
 end
 mean_mean_distance_cvl = mean(mean_distance_cvl);
 std_mean_distance_cvl = std(mean_distance_cvl)/sqrt(length(cvlcode));
% 
% mean_distance_cvl_left_left = mean(distance_cvl_left_left,2);
% mean_mean_distance_cvl_left_left = mean(mean_distance_cvl_left_left);
% std_distance_cvl_left_left = std(distance_cvl_left_left)/sqrt(7);
% mean_std_distance_cvl_left_left = mean(std_distance_cvl_left_left);
% 
% mean_distance_cvl_left_right = mean(distance_cvl_left_right,2);
% mean_mean_distance_cvl_left_right = mean(mean_distance_cvl_left_right);
% std_distance_cvl_left_right = std(distance_cvl_left_right)/sqrt(4);
% mean_std_distance_cvl_left_right = mean(std_distance_cvl_left_right);
% 
% 
% mean_distance_cvl_right_left = mean(distance_cvl_right_left,2);
% std_distance_cvl_right_left = std(distance_cvl_right_left);
% mean_distance_cvl_right_right = mean(distance_cvl_right_right,2);
% std_distance_cvl_right_right = std(distance_cvl_right_right);
% 

mean_distance_defocus = mean(distanceDefocus,2);
mean_mean_distance_defocus = mean(mean_distance_defocus);
std_mean_distance_defocus = std(mean_distance_defocus)/sqrt(length(defocuscode));


%% this code below useful to plot the scatter of average points comparing
% NV, defocus, and CVL subjects
figure
for i=1:length(nv_xy)
    scatter(nv_xy(i,1), nv_xy(i,2), 70,'b','o', 'filled');
    hold on
    scatter(cvl_xy(i,1), cvl_xy(i,2), 70,'r','o', 'filled');
    if i<16
        scatter(defocus_xy(i,1), defocus_xy(i,2), 70,'g','o', 'filled');    
    end
    
end
xlabel('Horizontal gaze offset distance (degrees)',  'FontSize', 18, 'FontName', 'Arial');
ylabel('Vertical gaze offset distance (degrees)', 'FontSize', 18, 'FontName', 'Arial');
set(gca,'xlim',[-10.5 3]);
line ([-10.5 3] , [0 0], 'LineStyle', '--', 'Color', 'k');
line ([0 0] , [-8 10], 'LineStyle', '--', 'Color', 'k');
box off


%%

% %
% % % figure
% % % scatter([1:length(all_means)],mean_distance,15,'b','o','filled');hold on;
% % % scatter(14,mean_mean_distance,20,'b','o', 'filled');
% % % scatter([1:length(all_means20)],mean_distance2,15,'g','o','filled');hold on;
% % % scatter(14,mean_mean_distance2,20,'g','o', 'filled');
% % % scatter([1:length(offset_mean)-1],mean_distance_cvl(1:end-1),15,'r','o','filled');hold on;
% % % scatter(14,mean_mean_distance_cvl,20,'r','o', 'filled');
% % 
errorbar(1:length(nvcode), mean_distance, std_distance, 'bx');
hold on;
%errorbar(1.3:length(nvcode)+.3, mean_distance2, std_distance2, 'gx');
errorbar(1.2:length(nvcode)+.2, mean_distance_cvl(1:end), std_distance_cvl, 'rx');
box off
xlabel('Subject index',  'FontSize', 18, 'FontName', 'Arial');
ylabel('Average gaze offset distance (degrees)', 'FontSize', 18, 'FontName', 'Arial');
set(gca, 'FontSize', 16', 'FontName', 'Arial');
legend({'NV subjects', 'CVL subject w/ similar age'});
legend boxoff
% % 
% % 
% % 
% % 
figure;
bar(1,mean_mean_distance, 'b');
hold on;
%bar(2, mean_mean_distance_defocus, 'g');
%bar(2, mean_mean_distance_defocus1, 'g');
%bar(3, mean_mean_distance_defocus2, 'g');
%bar(4, mean_mean_distance_defocus3, 'g');
bar(2, mean_mean_distance_defocus4, 'g');
%bar(6, mean_mean_distance_defocus5, 'g');

bar(3, mean_mean_distance_cvl, 'r');
%errorbar(2, mean_mean_distance2, std_mean_distance2, 'kx');
errorbar(1, mean_mean_distance, std_mean_distance, 'kx');
%errorbar(2, mean_mean_distance_defocus, std_mean_distance_defocus, 'kx');
%errorbar(2, mean_mean_distance_defocus1, std_mean_distance_defocus1, 'kx');
%errorbar(3, mean_mean_distance_defocus2, std_mean_distance_defocus2, 'kx');
%errorbar(4, mean_mean_distance_defocus3, std_mean_distance_defocus3, 'kx');
errorbar(2, mean_mean_distance_defocus4, std_mean_distance_defocus4, 'kx');
%errorbar(6, mean_mean_distance_defocus5, std_mean_distance_defocus5, 'kx');


errorbar(3, mean_mean_distance_cvl, std_mean_distance_cvl, 'kx');
box off
xlabel('Subject index', 'FontSize', 18, 'FontName', 'Arial');
ylabel('Average gaze offset distance (degrees)', 'FontSize', 18, 'FontName', 'Arial');
set(gca, 'FontSize', 16', 'FontName', 'Arial', 'xTick',[]);
legend({'NV subjects', 'Defocus subjects', 'CVL subjects'});
%legend({'NV subjects' , 'Hemianopes'});
set(gca, 'XTickLabel', {'NV', 'Defocus', 'CVL'});
set(gca, 'XTick', [1 2 3]);
legend boxoff
% % % 
% % % figure;
% % % for i=1:length(all_means)
% % %     
% % %    
% % %   mymean= mean(all_means{i})  
% % %   
% % %    % for j=1:length(all_means{1})
% % %   % scatter( all_means{i}(j,1), all_means{i}(j,2),5,[randi(i)/13 0.5 randi(i)/13],'o','filled');
% %    scatter( mymean(1,1), mymean(1,2),7,[randi(i)/13 0.5 randi(i)/13],'o','filled');
% %     
% %   hold on
% %    % end
% % end
% 
% figure;
% bar(1,mean_mean_distance_left, 'b');
% hold on;
% %bar(2, mean_mean_distance2, 'g');
% bar(2, mean_mean_distance_cvl_left_left, 'r');
% bar(3,mean_mean_distance_cvl_left_right, 'g');
% %errorbar(2, mean_mean_distance2, std_mean_distance2, 'kx');
% errorbar(1, mean_mean_distance_left, mean_std_distance_left, 'kx');
% errorbar(2, mean_mean_distance_cvl_left_left, mean_std_distance_cvl_left_left, 'kx');
% % bar(4, mean_distance_right_cvl, 'r');
% errorbar(3, mean_mean_distance_cvl_left_right, mean_std_distance_cvl_left_right, 'kx');
% % errorbar(4, mean_distance_cvl_right, std_distance_cvl_right, 'kx');
% box off
% xlabel('Group category', 'FontSize', 18, 'FontName', 'Arial');
% ylabel('Left shift in horizontal gaze offset (pixels)', 'FontSize', 18, 'FontName', 'Arial');
% set(gca, 'FontSize', 16', 'FontName', 'Arial', 'xTick',[]);
% %legend({'NV subject - Analysis 1', 'NV subject - Analysis 2', 'CVL subject w/ similar age'});
% legend({'NV subjects - Analysis 1', 'Left hemianopes', 'Right hemianopes'});
% 
% legend boxoff
