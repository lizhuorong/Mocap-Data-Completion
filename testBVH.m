clear all,

addpath('./matlab');
addpath('./NDLUTIL0p161');
% name = 'Melingkar.bvh';
% name = 'Maju Undur.bvh';
name = 'fastsong7ed_Take_001.bvh';

frame_step = 1;
write_video = 1;

[skel, channels, frameLength] = bvhReadFile(name);
PlaySkelData(skel, channels, frameLength, strcat(name(1:end-3),'mp4'));

% alljoint_xyz_frames = [];
% for i=1:size(channels,1)
%     tmp = bvh2xyz(skel, channels(i,:), 0);
%     tmp = [tmp(:,1),tmp(:,3),tmp(:,2)]';
%     alljoint_xyz_frames = [alljoint_xyz_frames,tmp(:)];
% end

alljoint_xyz_frames = dlmread('D:\codes\BU\Mycodes\motion\results_rendering\EXP1_ms_entries\Fast_song_7\updated\T0\rec_10.txt')';% 3 values/matrix
% alljoint_xyz_frames = dlmread('D:\codes\BU\Mycodes\motion\results_rendering\EXP2_missing_joints_Chaimue\Maju\RWrist\Maju Undur.txt')';% check txt
% alljoint_xyz_frames = dlmread('D:\codes\BU\Mycodes\motion\results_rendering\EXP5_2D3D\3D\F_ms_frame\combined\35.txt')';% rendering



%%%%%%%% Note that all the xyz of joints/markers are stored in the matrix
%%%%%%%% of "alljpoint_xyz_frames". Thus you should process it by your
%%%%%%%% interpolation experiments and then replace it with your results here. 
%%%%%%%% After that, you continue to render your updated motion data. The
%%%%%%%% program outputs two videos that are associated with the original
%%%%%%%% bvh file and your updated motion data respectively. Finally, a
%%%%%%%% figure will show the original motion and your results
%%%%%%%% simultaneously for comparison.

%% Initialise figure

xlim = get(gca, 'xlim');
minY1 = xlim(1);
maxY1 = xlim(2);
ylim = get(gca, 'ylim');
minY3 = ylim(1);
maxY3 = ylim(2);
zlim = get(gca, 'zlim');
minY2 = zlim(1);
maxY2 = zlim(2);
cla,
xlim = [minY1 maxY1];
ylim = [minY3 maxY3];
zlim = [minY2 maxY2];
set(gca, 'xlim', xlim, 'ylim', ylim, 'zlim', zlim);

Njoints = length(skel.tree);

for ff = 1:1
      for nn = 1:Njoints
            id = skel.tree(nn).id;
            pose = alljoint_xyz_frames(id*3+1:id*3+3,ff);
            hp(nn) = plot3(pose(1),pose(2),pose(3),'b.','MarkerSize',18);
            parent = skel.tree(nn).parent;
            if parent > 0
                parent_pose = alljoint_xyz_frames((parent-1)*3+1:parent*3,ff);
                hl(nn) = plot3([parent_pose(1) pose(1)],[parent_pose(2) pose(2)],[parent_pose(3) pose(3)],'r','LineWidth',2);
            end
      end
      drawnow
end

%% Run animation

if write_video, vidObj = VideoWriter(strcat(name(1:end-4),'updated.mp4'),'MPEG-4'); open(vidObj); end 

for ff = 2:frame_step:size(channels,1)
%       title(sprintf('%1.2f seconds',frameLength*ff))
      title(sprintf('frame %d',ff))
      for nn = 1:Njoints
            id = skel.tree(nn).id;
            pose = alljoint_xyz_frames(id*3+1:id*3+3,ff);
            hp(nn).XData = pose(1);
            hp(nn).YData = pose(2);
            hp(nn).ZData = pose(3);
            parent = skel.tree(nn).parent;
            if parent > 0
                    parent_pose = alljoint_xyz_frames((parent-1)*3+1:parent*3,ff);
                    hl(nn).XData = [parent_pose(1) pose(1)];
                    hl(nn).YData = [parent_pose(2) pose(2)];
                    hl(nn).ZData = [parent_pose(3) pose(3)];
            end
      end
      drawnow
      if write_video,
          Frame = getframe(gcf);
          writeVideo(vidObj,Frame);
          [X,map] = frame2im(Frame);
          imwrite(X,strcat(name(1:end-4),'update-frame-',num2str(ff),'.png'));
      end
end
if write_video, close(vidObj); end
close all,
VideoInCustomGUI(strcat(name(1:end-3),'mp4'),strcat(name(1:end-4),'updated.mp4'));


function PlaySkelData(skelStruct, channels, frameLength, videoname)

% Play skel MoCap data.
% FORMAT 
% ARG skel : the skeleton for the motion.
% ARG channels : the channels for the motion.
% ARG frameLength : the framelength for the motion.
% ARG videoname : output video.
if nargin < 3
  frameLength = 1/120;
  write_video=0;
elseif nargin == 4
    write_video = 1;
end
clf

handle = skelVisualise(channels(1, :), skelStruct);

% Get the limits of the motion.
xlim = get(gca, 'xlim');
minY1 = xlim(1);
maxY1 = xlim(2);
ylim = get(gca, 'ylim');
minY3 = ylim(1);
maxY3 = ylim(2);
zlim = get(gca, 'zlim');
minY2 = zlim(1);
maxY2 = zlim(2);
for i = 1:size(channels, 1)
  Y = skel2xyz(skelStruct, channels(i, :));
  minY1 = min([Y(:, 1); minY1]);
  minY2 = min([Y(:, 2); minY2]);
  minY3 = min([Y(:, 3); minY3]);
  maxY1 = max([Y(:, 1); maxY1]);
  maxY2 = max([Y(:, 2); maxY2]);
  maxY3 = max([Y(:, 3); maxY3]);
end
xlim = [minY1 maxY1];
ylim = [minY3 maxY3];
zlim = [minY2 maxY2];
set(gca, 'xlim', xlim, ...
         'ylim', ylim, ...
         'zlim', zlim);

% Play the motion
if write_video, vidObj = VideoWriter(videoname,'MPEG-4'); open(vidObj); end 
for j = 1:size(channels, 1)
  pause(frameLength)
  skelModify(handle, channels(j, :), skelStruct);
  if write_video,
      Frame = getframe;
      writeVideo(vidObj,Frame);
      [X,map] = frame2im(Frame);
      imwrite(X,strcat(videoname(1:end-4),'-frame-',num2str(j),'.png'));
  end
end
if write_video, close(vidObj); end
end


function VideoInCustomGUI(vname1,vname2)

videoSrc1 = vision.VideoFileReader(vname1, 'ImageColorSpace', 'Intensity');
videoSrc2 = vision.VideoFileReader(vname2, 'ImageColorSpace', 'Intensity');
[hFig, hAxes] = createFigureAndAxes();

insertButtons(hFig, hAxes, videoSrc1, videoSrc2);

% Initialize the display with the first frame of the video
[frame1,frame2] = getAndProcessFrame(videoSrc1,videoSrc2);
% Display input video frame on axis
showFrameOnAxis(hAxes.axis1, frame1);
showFrameOnAxis(hAxes.axis2, frame2);
end


    function [hFig, hAxes] = createFigureAndAxes()
        % Close figure opened by last run
        figTag = 'CVST_VideoOnAxis_9804532';
        close(findobj('tag',figTag));

        % Create new figure
        hFig = figure('numbertitle', 'off', 'name', 'Video In Custom GUI', ...
               'menubar','none', 'toolbar','none', 'resize', 'on', 'tag',figTag, ...
               'renderer','painters', 'position',[680 678 480 240],...
               'HandleVisibility','callback');
        % hide the handle to prevent unintended modifications of our custom UI

        % Create axes and titles
        hAxes.axis1 = createPanelAxisTitle(hFig,[0.1 0.2 0.36 0.6],'Original Motion'); % [X Y W H]
        hAxes.axis2 = createPanelAxisTitle(hFig,[0.5 0.2 0.36 0.6],'Interploation results');
    end

    
    function hAxis = createPanelAxisTitle(hFig, pos, axisTitle)
        % Create panel
        hPanel = uipanel('parent',hFig,'Position',pos,'Units','Normalized');

        % Create axis
        hAxis = axes('position',[0 0 1 1],'Parent',hPanel);
        hAxis.XTick = [];
        hAxis.YTick = [];
        hAxis.XColor = [1 1 1];
        hAxis.YColor = [1 1 1];
        % Set video title using uicontrol. uicontrol is used so that text
        % can be positioned in the context of the figure, not the axis.
        titlePos = [pos(1)+0.02 pos(2)+pos(3)+0.3 0.3 0.07];
        uicontrol('style','text', 'String', axisTitle, 'Units','Normalized',...
            'Parent',hFig,'Position', titlePos, 'BackgroundColor',hFig.Color);
    end

    
    function insertButtons(hFig,hAxes,videoSrc1,videoSrc2)
        % Play button with text Start/Pause/Continue
        uicontrol(hFig,'unit','pixel','style','pushbutton','string','Start',...
                'position',[10 10 75 25], 'tag','PBButton123','callback',...
                {@playCallback,videoSrc1,videoSrc2,hAxes});
        % Exit button with text Exit
        uicontrol(hFig,'unit','pixel','style','pushbutton','string','Exit',...
                'position',[100 10 50 25],'callback', ...
                {@exitCallback,videoSrc1,videoSrc2,hFig});
    end

    
    function playCallback(hObject,~,videoSrc1,videoSrc2,hAxes)
       try
            % Check the status of play button
            isTextStart = strcmp(hObject.String,'Start');
            isTextCont  = strcmp(hObject.String,'Continue');
            if isTextStart
               % Two cases: (1) starting first time, or (2) restarting
               % Start from first frame
               if isDone(videoSrc1)
                  reset(videoSrc1);
                  reset(videoSrc2);
               end
            end
            if (isTextStart || isTextCont)
                hObject.String = 'Pause';
            else
                hObject.String = 'Continue';
            end
            % Rotate input video frame and display original and rotated frames on figure
            
            while strcmp(hObject.String, 'Pause') && ~isDone(videoSrc1)
                % Get input video frame and rotated frame
                [frame,rotatedImg] = getAndProcessFrame(videoSrc1,videoSrc2);
                % Display input video frame on axis
                showFrameOnAxis(hAxes.axis1, frame);
                % Display rotated video frame on axis
                showFrameOnAxis(hAxes.axis2, rotatedImg);
            end

            % When video reaches the end of file, display "Start" on the play button.
            if isDone(videoSrc1)
               hObject.String = 'Start';
            end
       catch ME
           % Re-throw error message if it is not related to invalid handle
           if ~strcmp(ME.identifier, 'MATLAB:class:InvalidHandle')
               rethrow(ME);
           end
       end
    end

    
    function [frame1,frame2] = getAndProcessFrame(videoSrc1,videoSrc2)
        % Read input video frame
        frame1 = step(videoSrc1);
        frame2 = step(videoSrc2);
    end

    
    function exitCallback(~,~,videoSrc1,videoSrc2,hFig)
        % Close the video file
        release(videoSrc1);
        release(videoSrc2);
        % Close the figure window
        close(hFig);
    end    

% function [skeleton,time] = loadbvh(fname,varargin)
% %% LOADBVH  Load a .bvh (Biovision) file.
% %
% % Loads BVH file specified by FNAME (with or without .bvh extension)
% % and parses the file, calculating joint kinematics and storing the
% % output in SKELETON.
% %
% % Optional argument 'delim' allows for setting the delimiter between
% % fields. E.g., for a tab-separated BVH file:
% %
% %    skeleton = loadbvh('louise.bvh','delim','\t')
% %
% % By default 'delim' is set to the space character.
% %
% % SKELETON is a structure with N elements, where N is the number of joints.
% % Each element consists of the following fields:
% %
% %     * `name` -- human-readable description of the joint
% %     * `nestdepth` -- how many joints away from the origin
% %     * `parent` -- index to parent joint
% %     * `offset` -- translation from previous joint to current joint
% %     * `Nchannels` -- number of channels describing pose
% %                      (usually 3 for rotations only or 6 for 6DOF pose)
% %     * `Nframes` -- number of time samples
% %     * `order` -- the Euler angle order; e.g. [3 1 2] = ZXY
% %     * `Dxyz` -- XYZ displacements (directly from `trans` matrix)
% %     * `rxyz` -- XYZ rotations (to calculate `trans` matrix)
% %     * `trans` -- transformation matrix 
% %
% % Some details on the BVH file structure are given in "Motion Capture File
% % Formats Explained": http://www.dcs.shef.ac.uk/intranet/research/resmes/CS0111.pdf
% % But most of it is fairly self-evident.
% 
% %% Options
% 
% p = inputParser;
% p.addParameter('delim',' ');
% parse(p,varargin{:});
% opt = p.Results;
% 
% %% Load and parse header data
% %
% % The file is opened for reading, primarily to extract the header data (see
% % next section). However, I don't know how to ask Matlab to read only up
% % until the line "MOTION", so we're being a bit inefficient here and
% % loading the entire file into memory. Oh well.
% 
% % add a file extension if necessary:
% if ~strncmpi(fliplr(fname),'hvb.',4)
%   fname = [fname,'.bvh'];
% end
% 
% fid = fopen(fname);
% C = textscan(fid,'%s');
% fclose(fid);
% C = C{1};
% 
% 
% %% Parse data
% %
% % This is a cheap tokeniser, not particularly clever.
% % Iterate word-by-word, counting braces and extracting data.
% 
% % Initialise:
% skeleton = [];
% ii = 1;
% nn = 0;
% brace_count = 1;
% 
% while ~strcmp( C{ii} , 'MOTION' )
%   
%   ii = ii+1;
%   token = C{ii};
%   
%   if strcmp( token , '{' )
%     
%     brace_count = brace_count + 1;
%     
%   elseif strcmp( token , '}' )
%     
%     brace_count = brace_count - 1;
%     
%   elseif strcmp( token , 'OFFSET' )
%     
%     skeleton(nn).offset = [str2double(C(ii+1)) ; str2double(C(ii+2)) ; str2double(C(ii+3))];
%     ii = ii+3;
%     
%   elseif strcmp( token , 'CHANNELS' )
%     
%     skeleton(nn).Nchannels = str2double(C(ii+1));
%     
%     % The 'order' field is an index corresponding to the order of 'X' 'Y' 'Z'.
%     % Subtract 87 because char numbers "X" == 88, "Y" == 89, "Z" == 90.
%     if skeleton(nn).Nchannels == 3
%       skeleton(nn).order = [C{ii+2}(1),C{ii+3}(1),C{ii+4}(1)]-87;
%     elseif skeleton(nn).Nchannels == 6
%       skeleton(nn).order = [C{ii+5}(1),C{ii+6}(1),C{ii+7}(1)]-87;
%     else
%       error('Not sure how to handle not (3 or 6) number of channels.')
%     end
%     
%     if ~all(sort(skeleton(nn).order)==[1 2 3])
%       error('Cannot read channels order correctly. Should be some permutation of [''X'' ''Y'' ''Z''].')
%     end
% 
%     ii = ii + skeleton(nn).Nchannels + 1;
% 
%   elseif strcmp( token , 'JOINT' ) || strcmp( token , 'ROOT' )
%     % Regular joint
%     
%     nn = nn+1;
%     
%     skeleton(nn).name = C{ii+1};
%     skeleton(nn).nestdepth = brace_count;
% 
%     if brace_count == 1
%       % root node
%       skeleton(nn).parent = 0;
%     elseif skeleton(nn-1).nestdepth + 1 == brace_count;
%       % if I am a child, the previous node is my parent:
%       skeleton(nn).parent = nn-1;
%     else
%       % if not, what is the node corresponding to this brace count?
%       prev_parent = skeleton(nn-1).parent;
%       while skeleton(prev_parent).nestdepth+1 ~= brace_count
%         prev_parent = skeleton(prev_parent).parent;
%       end
%       skeleton(nn).parent = prev_parent;
%     end
%     
%     ii = ii+1;
%             
%   elseif strcmp( [C{ii},' ',C{ii+1}] , 'End Site' )
%     % End effector; unnamed terminating joint
%     %
%     % N.B. The "two word" token here is why we don't use a switch statement
%     % for this code.
%     
%     nn = nn+1;
%     
%     skeleton(nn).name = ' ';
%     skeleton(nn).offset = [str2double(C(ii+4)) ; str2double(C(ii+5)) ; str2double(C(ii+6))];
%     skeleton(nn).parent = nn-1; % always the direct child
%     skeleton(nn).nestdepth = brace_count;
%     skeleton(nn).Nchannels = 0;
%         
%   end
%   
% end
% 
% %% Initial processing and error checking
% 
% Nnodes = numel(skeleton);
% Nchannels = sum([skeleton.Nchannels]);
% Nchainends = sum([skeleton.Nchannels]==0);
% 
% % Calculate number of header lines:
% %  - 5 lines per joint
% %  - 4 lines per chain end
% %  - 4 additional lines (first one and last three)
% Nheaderlines = (Nnodes-Nchainends)*5 + Nchainends*4 + 4;
% 
% rawdata = importdata(fname,opt.delim,Nheaderlines);
% 
% if ~isstruct(rawdata)
%    error('Could not parse BVH file %s. Check the delimiter.',fname)
% end
% 
% index = strncmp(rawdata.textdata,'Frames:',7);
% Nframes = sscanf(rawdata.textdata{index},'Frames: %f');
% 
% index = strncmp(rawdata.textdata,'Frame Time:',10);
% frame_time = sscanf(rawdata.textdata{index},'Frame Time: %f');
% 
% time = frame_time*(0:Nframes-1);
% 
% if size(rawdata.data,2) ~= Nchannels
%   error('Error reading BVH file: channels count does not match.')
% end
% 
% if size(rawdata.data,1) ~= Nframes
%   warning('LOADBVH:frames_wrong','Error reading BVH file: frames count does not match; continuing anyway.')
%   Nframes = size(rawdata.data,1);
% end
% 
% %% Load motion data into skeleton structure
% %
% % We have three possibilities for each node we come across:
% % (a) a root node that has displacements already defined,
% %     for which the transformation matrix can be directly calculated;
% % (b) a joint node, for which the transformation matrix must be calculated
% %     from the previous points in the chain; and
% % (c) an end effector, which only has displacement to calculate from the
% %     previous node's transformation matrix and the offset of the end
% %     joint.
% %
% % These are indicated in the skeleton structure, respectively, by having
% % six, three, and zero "channels" of data.
% % In this section of the code, the channels are read in where appropriate
% % and the relevant arrays are pre-initialised for the subsequent calcs.
% 
% channel_count = 0;
% 
% for nn = 1:Nnodes
%     
%   if skeleton(nn).Nchannels == 6 % root node
%     
%     % assume translational data is always ordered XYZ
%     skeleton(nn).Dxyz = repmat(skeleton(nn).offset,[1 Nframes]) + rawdata.data(:,channel_count+[1 2 3])';
%     skeleton(nn).rxyz(skeleton(nn).order,:) = rawdata.data(:,channel_count+[4 5 6])';
%         
%     % Kinematics of the root element:
%     skeleton(nn).trans = nan(4,4,Nframes);
%     for ff = 1:Nframes
%       skeleton(nn).trans(:,:,ff) = transformation_matrix(skeleton(nn).Dxyz(:,ff) , skeleton(nn).rxyz(:,ff) , skeleton(nn).order);
%     end
%     
%   elseif skeleton(nn).Nchannels == 3 % joint node
%         
%     skeleton(nn).rxyz(skeleton(nn).order,:) = rawdata.data(:,channel_count+[1 2 3])';
%     skeleton(nn).Dxyz  = nan(3,Nframes);
%     skeleton(nn).trans = nan(4,4,Nframes);
%     
%   elseif skeleton(nn).Nchannels == 0 % end node
%     skeleton(nn).Dxyz  = nan(3,Nframes);
%   end
%   
%   channel_count = channel_count + skeleton(nn).Nchannels;
%   skeleton(nn).Nframes = Nframes;
%   
% end
% 
% 
% %% Calculate kinematics
% %
% % No calculations are required for the root nodes.
% 
% % For each joint, calculate the transformation matrix and for convenience
% % extract each position in a separate vector.
% for nn = find([skeleton.parent] ~= 0 & [skeleton.Nchannels] ~= 0)
%   
%   parent = skeleton(nn).parent;
%   
%   for ff = 1:Nframes
%     transM = transformation_matrix( skeleton(nn).offset , skeleton(nn).rxyz(:,ff) , skeleton(nn).order );
%     skeleton(nn).trans(:,:,ff) = skeleton(parent).trans(:,:,ff) * transM;
%     skeleton(nn).Dxyz(:,ff) = skeleton(nn).trans([1 2 3],4,ff);
%   end
% 
% end
% 
% % For an end effector we don't have rotation data;
% % just need to calculate the final position.
% for nn = find([skeleton.Nchannels] == 0)
%   
%   parent = skeleton(nn).parent;
%   
%   for ff = 1:Nframes
%     transM = skeleton(parent).trans(:,:,ff) * [eye(3), skeleton(nn).offset; 0 0 0 1];
%     skeleton(nn).Dxyz(:,ff) = transM([1 2 3],4);
%   end
% 
% end
% 
% end
% 
% 
% function transM = transformation_matrix(displ,rxyz,order)
% % Constructs the transformation for given displacement, DISPL, and
% % rotations RXYZ. The vector RYXZ is of length three corresponding to
% % rotations around the X, Y, Z axes.
% %
% % The third input, ORDER, is a vector indicating which order to apply
% % the planar rotations. E.g., [3 1 2] refers applying rotations RXYZ
% % around Z first, then X, then Y.
% %
% % Years ago we benchmarked that multiplying the separate rotation matrices
% % was more efficient than pre-calculating the final rotation matrix
% % symbolically, so we don't "optimise" by having a hard-coded rotation
% % matrix for, say, 'ZXY' which seems more common in BVH files.
% % Should revisit this assumption one day.
% %
% % Precalculating the cosines and sines saves around 38% in execution time.
% 
% c = cosd(rxyz);
% s = sind(rxyz);
% 
% RxRyRz(:,:,1) = [1 0 0; 0 c(1) -s(1); 0 s(1) c(1)];
% RxRyRz(:,:,2) = [c(2) 0 s(2); 0 1 0; -s(2) 0 c(2)];
% RxRyRz(:,:,3) = [c(3) -s(3) 0; s(3) c(3) 0; 0 0 1];
% 
% rotM = RxRyRz(:,:,order(1))*RxRyRz(:,:,order(2))*RxRyRz(:,:,order(3));
% 
% transM = [rotM, displ; 0 0 0 1];
% 
% end
% 
% 
% % [sk,Ti] = loadbvh('Bokmue.bvh');
% % for j=2:length(Ti)
% %     
% %     cla
% %     for i=1:length(sk)
% %         % -------------------------------------------------------------------------
% %         % Display the bones
% %         if sk(i).parent
% %             s=sk(i).t_xyz(:,j); sp = sk(sk(i).parent).t_xyz(:,j);
% %             plot3([s(1) sp(1)],[s(2) sp(2)],[s(3) sp(3)],'r','LineWidth',5);
% %         end
% %         plot3(sk(i).t_xyz(1,j),sk(i).t_xyz(2,j),sk(i).t_xyz(3,j),'g.','MarkerSize',10);
% %         if sk(i).Nchannels
% %             plotJoint(sk(i).t_xyz(:,j)',sk(i).T(1:3,1:3,j)*sk(i).R0,15)
% %         end
% %     end
% %     axis off;
% %     axis equal
% %     shading interp;
% %     light
% %     lighting phong
% %     
% %     drawnow
% %     
% % end
% % Ti
% % 
% % 
% % function plotJoint(pt,R,sz)
% % R = permute(R,[3 1 2]);
% % plot3([pt(:,1) pt(:,1)+R(:,1,1)*sz],[pt(:,2) pt(:,2)+R(:,2,1)*sz],[pt(:,3) pt(:,3)+R(:,3,1)*sz],'m','LineWidth',2);
% % plot3([pt(:,1) pt(:,1)+R(:,1,2)*sz],[pt(:,2) pt(:,2)+R(:,2,2)*sz],[pt(:,3) pt(:,3)+R(:,3,2)*sz],'g','LineWidth',2);
% % plot3([pt(:,1) pt(:,1)+R(:,1,3)*sz],[pt(:,2) pt(:,2)+R(:,2,3)*sz],[pt(:,3) pt(:,3)+R(:,3,3)*sz],'b','LineWidth',2);
% % end
% 
% 
