% Plot a phylogenetic tree. Like Matlab's 'plot', but enable different branch colors.
function handles = phytree_plot_load(tr,varargin)
%PLOT renders a phylogenetic tree.
%
%   PLOT(TREE) renders a phylogenetic tree object into a MATLAB figure as a
%   phylogram. The significant distances between branches and nodes are in
%   horizontal direction, vertical coordinates are accommodated only for
%   display purposes. Handles to graph elements are stored in the
%   'UserData' figure field, such that graphic properties can be easily
%   modified.
%
%   PLOT(TREE,ACTIVEBRANCHES) hides the non'active branches and all their
%   descendants. ACTIVEBRANCHES is a logical array of size
%   [numBranches x 1] indicating the active branches.
%
%   PLOT(...,'TYPE',type) selects the method to render the phylogenetic
%   tree. Options are: 'square' (default), 'angular', and 'radial'.
%
%   PLOT(...,'ORIENTATION',orient) will orient the phylogenetic tree within
%   the figure window. Options are: 'top', 'bottom', 'left' (default), and,
%   'right'. Orientation parameter is valid only for phylograms or
%   cladograms.
%
%   PLOT(...,'BRANCHLABELS',value) hides/unhides branch labels. Options are
%   true or false. Branch labels are placed next to the branch node.
%   Defaults to false (true) when TYPE is (is not) 'radial'.
%
%   PLOT(...,'LEAFLABELS',value) hides/unhides leaf labels. Options are
%   true or false. Leaf labels are placed next to the leaf nodes. Defaults
%   to false (true) when TYPE is (is not) 'radial'.
%
%   PLOT(...,'TERMINALLABELS',value) hides/unhides terminal labels. Options
%   are true (default) or false. Terminal labels are placed over the axis
%   tick labels, ignored when 'radial' type is used.
%
%   H = PLOT(...) returns a structure with handles to the graph elements.
%
%   Example:
%
%       tr = phytreeread('pf00002.tree')
%       plot(tr,'type','radial')
%
%       % Graph element properties can be modified as follows:
%
%       h=get(gcf,'UserData')
%       set(h.branchNodeLabels,'FontSize',6,'Color',[.5 .5 .5])
%
%   See also PHYTREE, PHYTREE/VIEW, PHYTREEREAD, PHYTREETOOL, SEQLINKAGE.

% Copyright 2003-2006 The MathWorks, Inc.
% $Revision: 1.1.6.10 $ $Author: batserve $ $Date: 2006/06/16 20:06:45 $
%
% Changed by Or Zuk 06/2008, to enable drawing colored branches

if numel(tr)~=1
    error('Bioinfo:phytree:plot:NoMultielementArrays',...
        'Phylogenetic tree must be an 1-by-1 object.');
end

% set defaults
dispBranchLabels = NaN;
dispLeafLabels = NaN;
dispTerminalLabels = false;
renderType = 'square';
orientation = 'left';
rotation = 0;
new_fig_flag = 1; % default is having a new figure 

tr = struct(tr); %tr is the tree
tr.numBranches = size(tr.tree,1); %does not include leafs

%Sets active or inactive leafs
if nargin>1 && islogical(varargin{1})
    activeBranches = varargin{1};
    argStart = 2;
else
    activeBranches = true(tr.numBranches,1);
    argStart = 1;
end

%Get input arguements
if nargin - argStart > 0
    if rem(nargin - argStart,2) == 1
        error('Bioinfo:phytree:plot:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'type','orientation','rotation',...
        'branchlabels','leaflabels','terminallabels', 'weights', 'new_fig'}; % New! added weights! 
    for j = argStart:2:nargin-argStart
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error('Bioinfo:phytree:plot:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:phytree:plot:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1 % type
                    oktypes={'square','angular','radial'};
                    l = strmatch(lower(pval),oktypes); %#ok
                    if isempty(l)
                        error('Bioinfo:phytree:plot:UnknownTypeName',...
                            'Unknown option for %s.',upper(okargs{k}));
                    else
                        if l==4
                            l=1;
                        end
                        renderType = oktypes{l};
                    end
                case 2 % orientation
                    oktypes={'left','right','top','bottom'};
                    l = strmatch(lower(pval),oktypes); %#ok
                    if isempty(l)
                        error('Bioinfo:phytree:plot:UnknownOrientation',...
                            'Unknown option for %s.',upper(okargs{k}));
                    else
                        orientation = oktypes{l};
                    end
                case 3 % rotation
                    if isreal(pval(1))
                        rotation = double(pval(1));
                    else
                        error('Bioinfo:phytree:plot:NotValidType',...
                            'ROTATION must be numeric and real');
                    end
                case 4 % branch labels
                    dispBranchLabels = opttf(pval);
                case 5 % leaf labels
                    dispLeafLabels = opttf(pval);
                case 6 % terminal labels
                    dispTerminalLabels = opttf(pval);
                case 7 % weights 
                    weights = pval; % numeric 
                case 8 % new figure false
                    new_fig_flag = 0;
            end
        end
    end
end

% set dependent defaults
if isnan(dispBranchLabels)
    if isequal(renderType,'radial')
        dispBranchLabels = true;
    else
        dispBranchLabels = false;
    end
end
if isnan(dispLeafLabels)
    if isequal(renderType,'radial')
        dispLeafLabels = true;
    else
        dispLeafLabels = false;
    end
end

%Computes branch locations in axes
tr = doBasicCalculations(tr,activeBranches,renderType);

nodeIndex   = 1:tr.numLabels;
leafIndex   = 1:tr.numLeaves;
branchIndex = tr.numLeaves+1:tr.numLabels;


% check empty names
for ind = nodeIndex
    if isempty(tr.names{ind})
        if ind > tr.numLeaves
            tr.names{ind} = ['Branch ' num2str(ind-tr.numLeaves)];
        else
            tr.names{ind} = ['Leaf ' num2str(ind)];
        end
    end
end

% rendering graphic objects
if(new_fig_flag)
    fig = figure('Renderer','ZBuffer');
else 
    fig = gcf;
end
h.axes = axes; hold on;  % This is problematic: ruins subplot 
sepUnit = max(tr.x)*[-1/20 21/20];

% setting the axes
switch renderType
    case {'square','angular'}
        switch orientation
            case 'left'
            	%Just consider this case
                set(h.axes,'YTick',1:numel(tr.terminalNodes),'Ydir','reverse',...
                    'YtickLabel','','YAxisLocation','Right')
                if dispTerminalLabels
                    set(h.axes,'Position',[.05 .10 .7 .85])
                else
                    set(h.axes,'Position',[.05 .10 .9 .85])
                end
                xlim(sepUnit);
                ylim([0 numel(tr.terminalNodes)+1]);
            case 'right'
                set(h.axes,'YTick',1:numel(tr.terminalNodes),'Xdir','reverse','Ydir','reverse',...
                    'YtickLabel','','YAxisLocation','Left')
                if dispTerminalLabels
                    set(h.axes,'Position',[.25 .10 .7 .85])
                else
                    set(h.axes,'Position',[.05 .10 .9 .85])
                end
                xlim(sepUnit);
                ylim([0 numel(tr.terminalNodes)+1]);
            case 'top'
                set(h.axes,'XTick',1:numel(tr.terminalNodes),...
                    'XtickLabel','','XAxisLocation','Top')
                if dispTerminalLabels
                    set(h.axes,'Position',[.10 .05 .85 .7])
                else
                    set(h.axes,'Position',[.10 .05 .85 .9])
                end
                ylim(sepUnit);
                xlim([0 numel(tr.terminalNodes)+1]);
            case 'bottom'
                set(h.axes,'XTick',1:numel(tr.terminalNodes),'Ydir','reverse',...
                    'XtickLabel','','XAxisLocation','Bottom')
                if dispTerminalLabels
                    set(h.axes,'Position',[.10 .25 .85 .7])
                else
                    set(h.axes,'Position',[.10 .05 .85 .9])
                end
                ylim(sepUnit);
                xlim([0 numel(tr.terminalNodes)+1]);
        end
    case 'radial'
        set(h.axes,'XTick',[],'YTick',[])
        set(h.axes,'Position',[.05 .05 .9 .9])
        dispTerminalLabels = false;
        axis equal
end

%Set up coloring scheme for lineage segments
lineages = tr.numBranches + tr.numLeaves;
lineageSegs = tr.orderedLineageSegs;
allSegStates =[];
for lin = 1:lineages
    allSegStates = [allSegStates, lineageSegs.segStates{lin}];
end
maxState = max(allSegStates)+1;
colorMatrix = colormap(cool(maxState)); %jet also works nice

% Set up for drawing lines
switch renderType
    case 'square'
        X = tr.x([nodeIndex;repmat([tr.par(1:tr.numLabels-1) tr.numLabels],2,1)]);
        Y = tr.y([repmat(nodeIndex,2,1);[tr.par(1:tr.numLabels-1) tr.numLabels]]);
        
        %Compute branch segment coordinates
        lineageSegs = getLineageSegmentCoordinates(lineageSegs,X,Y);
        
        hold on;
    case 'angular'
        X = tr.x([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]);
        Y = tr.y([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]);
    case 'radial'
        R = tr.x;
        A = tr.y / numel(tr.terminalNodes)*2*pi+rotation*pi/180;
        tr.x = R .* sin(A);
        tr.y = R .* cos(A);
        X = tr.x([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]);
        Y = tr.y([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]);
end

% plotting lines
switch orientation % Note: do not use bottom/top orientation with radial type!!!
    case {'left','right'}
        %h.BranchLines = plot(X,Y,'-k'); //Matlab's way
        lineCounter = 1;
        for lin = 1:(lineages-1)
            % plot horizontal line segments 
            numSegs = size(lineageSegs.segLengths{lin},2);
            for seg = 1: numSegs
                %Color all black
                %h.BranchLines(lineCounter) = plot(lineageSegs.x{lin}(1:2,seg), Y(1:2,lin), 'color', 'k', 'LineWidth', 2.0); lineCounter = lineCounter + 1;
                %Color according to color scheme
                h.BranchLines(lineCounter) = plot(lineageSegs.x{lin}(1:2,seg), Y(1:2,lin), 'color', colorMatrix(lineageSegs.segStates{lin}(seg)+1,:), 'LineWidth', 2.0); lineCounter = lineCounter + 1;
            end
            % plot vertical line segments
            h.BranchLines(lineCounter) = plot(X(2:3,lin), Y(2:3,lin), 'color', colorMatrix(lineageSegs.segStates{lin}(numSegs)+1,:), 'LineWidth', 2.0); lineCounter = lineCounter + 1;
            %h.BranchLines(lineCounter) = plot(X(2:3,lin), Y(2:3,lin), 'color', 'k', 'LineWidth', 2.0); lineCounter = lineCounter + 1;
        end
        
    case {'top','bottom'}
        if(strmatch('radial', renderType, 'exact'))
            h.BranchLines = plot(X,Y,'-k');
        else
            h.BranchLines = plot(Y,X,'-k');
            tmp_X = Y; Y = X; X = tmp_X;
        end
end
delete(h.BranchLines(~tr.activeNodes));
h.BranchLines = h.BranchLines(tr.activeNodes);

% Set line widths based on weights
if(exist('weights', 'var')) % New: check if we have also weights information.
    n = size(X, 2); linewidth = 0.5; 
    color_mat = colormap; ncolors = size(color_mat,1); % 64 colors
    min_w = min(weights(weights>0)), max_w = max(weights(:)), range_w = max_w - min_w; 
    k_max = min(range_w / min_w, 10); % set maximum line width
    
    for i=1:n % loop over branches and color them
        line([X(1,i), X(2,i)], [Y(1,i), Y(2,i)], ...
            'Color', color_mat(   max(1,ceil( ncolors * (weights(i) - min_w) / range_w  )),: ), ...
            'LineWidth', linewidth + (weights(i) - min_w) * (k_max-1)*linewidth / range_w  ); % draw the horizonal lines
        switch renderType
            case 'square'
                line([X(2,i), X(3,i)], [Y(2,i), Y(3,i)], ...
                    'Color', color_mat(   max(1,ceil( ncolors * (weights(i) - min_w) / range_w  )),: ), ...
                    'LineWidth', linewidth + (weights(i) - min_w) * (k_max-1)*linewidth / range_w   ); % draw the vertical lines
        end
    end
end

% drawing nodes
switch renderType
    case {'square','angular'}
        switch orientation
            case {'left','right'}
                %h.BranchDots = plot(tr.x(branchIndex(tr.activeNodes(branchIndex))),...
                    %tr.y(branchIndex(tr.activeNodes(branchIndex))),'o',...
                    %'MarkerSize',nodeMarkerSize,'MarkerEdgeColor','k',...
                    %'MarkerFaceColor','b');
                %h.LeafDots = plot(tr.x(leafIndex(tr.activeNodes(leafIndex))),...
                    %tr.y(leafIndex(tr.activeNodes(leafIndex))),'square',...
                    %'MarkerSize',4,'MarkerEdgeColor','k',...
                    %'MarkerFaceColor','w');
                 leafCounter = 1;
                 for i = 1:tr.numLeaves
                    %leafName = tr.names(i);
                    %leafNameBlocks = regexp(leafName, '-', 'split');
                    %leafState = str2num(leafNameBlocks{1,1}{1,1}); %should be two if only one '_' in name
                    
                    %Paint black
                    %h.leafDots(leafCounter) = plot(tr.x(i), tr.y(i), 'square', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
                    
                    %Paint by state
                    leafState = lineageSegs.segStates{i}(1)+1;
                    h.leafDots(leafCounter) = plot(tr.x(i), tr.y(i), 'square', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorMatrix(leafState,:));
                    
                    %if (leafState == 1) 
                     %   h.leafDots(leafCounter) = plot(tr.x(i), tr.y(i), 'square', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorMatrix(20,:));
                    %else
                     %   h.leafDots(leafCounter) = plot(tr.x(i), tr.y(i), 'square', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorMatrix(1,:));
                    %end
                    leafCounter = leafCounter + 1;
                 end
            case {'top','bottom'}
                h.BranchDots = plot(tr.y(branchIndex(tr.activeNodes(branchIndex))),...
                    tr.x(branchIndex(tr.activeNodes(branchIndex))),'o',...
                    'MarkerSize',5,'MarkerEdgeColor','k',...
                    'MarkerFaceColor','b');
                h.LeafDots = plot(tr.y(leafIndex(tr.activeNodes(leafIndex))),...
                    tr.x(leafIndex(tr.activeNodes(leafIndex))),'square',...
                    'MarkerSize',4,'MarkerEdgeColor','k',...
                    'MarkerFaceColor','w');
        end
    case 'radial'
        h.BranchDots = plot(tr.x(branchIndex(tr.activeNodes(branchIndex))),...
            tr.y(branchIndex(tr.activeNodes(branchIndex))),'o',...
            'MarkerSize',5,'MarkerEdgeColor','k',...
            'MarkerFaceColor','b');
        h.LeafDots = plot(tr.x(leafIndex(tr.activeNodes(leafIndex))),...
            tr.y(leafIndex(tr.activeNodes(leafIndex))),'square',...
            'MarkerSize',4,'MarkerEdgeColor','k',...
            'MarkerFaceColor','w');
end

if(exist('weights', 'var'))
    colorbar('WestOutside', 'YTickLabel', [min_w:(range_w/10):max_w]);
end

% resize figure if needed
switch renderType
    case {'square','angular'}
        switch orientation
            case {'left','right'}
                correctFigureSize(fig, 15 * numel(tr.terminalNodes),0);
                fontRatio = max(get(fig,'Position').*[0 0 0 1])/numel(tr.terminalNodes);
            case {'top','bottom'}
                correctFigureSize(fig, 0, 15 * numel(tr.terminalNodes));
                fontRatio = max(get(fig,'Position').*[0 0 1 0])/numel(tr.terminalNodes);
        end
    case 'radial'
        temp = 10/pi*numel(tr.terminalNodes);
        correctFigureSize(fig,temp,temp);
        fontRatio = max(get(fig,'Position').*[0 0 1 0])/numel(tr.terminalNodes);
end

set(h.axes,'Fontsize',min(9,ceil(fontRatio/1.5)));

% set branch node labels
X = tr.x(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels)));
Y = tr.y(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels)));
switch renderType
    case {'square','angular'}
        switch orientation
            case {'left'}
                h.branchNodeLabels = text(X+sepUnit(1)/2,Y,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
                set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
                set(h.branchNodeLabels,'vertical','bottom')
                set(h.branchNodeLabels,'horizontal','right')
                set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
            case {'right'}
                h.branchNodeLabels = text(X+sepUnit(1)/2,Y,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
                set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
                set(h.branchNodeLabels,'vertical','bottom')
                set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
            case {'top'}
                h.branchNodeLabels = text(Y,X-sepUnit(1)/2,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
                set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
                set(h.branchNodeLabels,'vertical','bottom','Rotation',30)
                set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
            case {'bottom'}
                h.branchNodeLabels = text(Y,X+sepUnit(1)/2,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
                set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
                set(h.branchNodeLabels,'vertical','bottom','Rotation',30)
                set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
        end
    case 'radial'
        h.branchNodeLabels = text(X,Y,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
        set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
        set(h.branchNodeLabels,'vertical','bottom')
        set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio*1.2)));
        for ind = 1:numel(h.branchNodeLabels)
            if X(ind)<0
                set(h.branchNodeLabels(ind),'horizontal','right')
                set(h.branchNodeLabels(ind),'Position',get(h.branchNodeLabels(ind),'Position')+[sepUnit(1)/2 0 0])
            else
                set(h.branchNodeLabels(ind),'horizontal','left')
                set(h.branchNodeLabels(ind),'Position',get(h.branchNodeLabels(ind),'Position')-[sepUnit(1)/2 0 0])
            end
        end
end

% set leaf nodes labels
X = tr.x(leafIndex(tr.activeNodes(1:tr.numLeaves)));
Y = tr.y(leafIndex(tr.activeNodes(1:tr.numLeaves)));
switch renderType
    case {'square','angular'}
        switch orientation
            case {'left'}
                h.leafNodeLabels = text(X-sepUnit(1)/2,Y,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
                set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
                set(h.leafNodeLabels,'horizontal','left')
                set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
            case {'right'}
                h.leafNodeLabels = text(X-sepUnit(1)/2,Y,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
                set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
                set(h.leafNodeLabels,'horizontal','right')
                set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
            case {'top'}
                h.leafNodeLabels = text(Y,X-sepUnit(1)/2,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
                set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
                set(h.leafNodeLabels,'horizontal','left','Rotation',60)
                set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
            case {'bottom'}
                h.leafNodeLabels = text(Y,X-sepUnit(1),tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
                set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
                set(h.leafNodeLabels,'horizontal','right','Rotation',60)
                set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
        end
    case 'radial'
        h.leafNodeLabels = text(X,Y,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
        set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
        set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio*1.2)));
        % textHeight = mean(cell2mat(get(h.leafNodeLabels,'Extent')))*[0 0 0 1]';
        for ind = 1:numel(h.leafNodeLabels)
            if X(ind)<0
                set(h.leafNodeLabels(ind),'horizontal','right')
                set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')+[sepUnit(1) 0 0])
            else
                set(h.leafNodeLabels(ind),'horizontal','left')
                set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')-[sepUnit(1) 0 0])
            end
            %             a=atan(Y(ind)/X(ind))*180/pi;
            %             if a > 0  a = max(0,a-60)/2; else
            %                       a = min(0,a+60)/2; end
            %             set(h.leafNodeLabels(ind),'Rotation',a)
        end
        [sortedY,hsY]=sort(Y);
        idx=hsY(X(hsY)>0 & sortedY>0);
        if numel(idx)
            extentY = get(h.leafNodeLabels(idx(1)),'Extent')*[0;0;0;1];
            positionY = get(h.leafNodeLabels(idx(1)),'Position')*[0;1;0];
            for i = 2:numel(idx)
                position = get(h.leafNodeLabels(idx(i)),'Position');
                positionY = max(positionY+extentY,position(2));
                position(2) = positionY;
                set(h.leafNodeLabels(idx(i)),'Position',position)
            end
        end
        idx=hsY(X(hsY)<0 & sortedY>0);
        if numel(idx)
            extentY = get(h.leafNodeLabels(idx(1)),'Extent')*[0;0;0;1];
            positionY = get(h.leafNodeLabels(idx(1)),'Position')*[0;1;0];
            for i = 2:numel(idx)
                position = get(h.leafNodeLabels(idx(i)),'Position');
                positionY = max(positionY+extentY,position(2));
                position(2) = positionY;
                set(h.leafNodeLabels(idx(i)),'Position',position)
            end
        end
        idx=flipud(hsY(X(hsY)>0 & sortedY<0));
        if numel(idx)
            extentY = get(h.leafNodeLabels(idx(1)),'Extent')*[0;0;0;1];
            positionY = get(h.leafNodeLabels(idx(1)),'Position')*[0;1;0];
            for i = 2:numel(idx)
                position = get(h.leafNodeLabels(idx(i)),'Position');
                positionY = min(positionY-extentY,position(2));
                position(2) = positionY;
                set(h.leafNodeLabels(idx(i)),'Position',position)
            end
        end
        idx=flipud(hsY(X(hsY)<0 & sortedY<0));
        if numel(idx)
            extentY = get(h.leafNodeLabels(idx(1)),'Extent')*[0;0;0;1];
            positionY = get(h.leafNodeLabels(idx(1)),'Position')*[0;1;0];
            for i = 2:numel(idx)
                position = get(h.leafNodeLabels(idx(i)),'Position');
                positionY = min(positionY-extentY,position(2));
                position(2) = positionY;
                set(h.leafNodeLabels(idx(i)),'Position',position)
            end
        end

end

% correct axis limits given the extent of labels
if dispBranchLabels
    E = cell2mat(get(h.branchNodeLabels,'Extent'));
    if strcmp(get(gca,'XDir'),'reverse')
        E(:,1) = E(:,1) - E(:,3);
    end
    if strcmp(get(gca,'YDir'),'reverse')
        E(:,2) = E(:,2) - E(:,4);
    end
    E=[E;[xlim*[1;0] ylim*[1;0] diff(xlim) diff(ylim)]];
    mins = min(E(:,[1 2]));
    maxs = max([sum(E(:,[1 3]),2) sum(E(:,[2 4]),2)]);
    axis([mins(1) maxs(1) mins(2) maxs(2)])
end

if dispLeafLabels
    E = cell2mat(get(h.leafNodeLabels,'Extent'));
    if strcmp(get(gca,'XDir'),'reverse')
        E(:,1) = E(:,1) - E(:,3);
    end
    if strcmp(get(gca,'YDir'),'reverse')
        E(:,2) = E(:,2) - E(:,4);
    end
    E=[E;[xlim*[1;0] ylim*[1;0] diff(xlim) diff(ylim)]];
    mins = min(E(:,[1 2]));
    maxs = max([sum(E(:,[1 3]),2) sum(E(:,[2 4]),2)]);
    axis([mins(1) maxs(1) mins(2) maxs(2)])
end

% set terminal nodes labels
switch renderType
    case {'square','angular'}
        X = tr.x(tr.terminalNodes) * 0;
        Y = tr.y(tr.terminalNodes);
        switch orientation
            case {'left'}
                X = X + max(xlim) - sepUnit(1)/2;
                h.terminalNodeLabels = text(X,Y,tr.names(tr.terminalNodes));
            case {'right'}
                X = X + max(xlim) - sepUnit(1)/2;
                h.terminalNodeLabels = text(X,Y,tr.names(tr.terminalNodes));
                set(h.terminalNodeLabels,'Horizontal','right')
            case {'top'}
                X = X + max(ylim) - sepUnit(1)/2;
                h.terminalNodeLabels = text(Y,X,tr.names(tr.terminalNodes));
                set(h.terminalNodeLabels,'Rotation',90)
            case {'bottom'}
                X = X + max(ylim) - sepUnit(1)/2;
                h.terminalNodeLabels = text(Y,X,tr.names(tr.terminalNodes));
                set(h.terminalNodeLabels,'Rotation',270)
        end
    case 'radial'
        h.terminalNodeLabels = text(0,0,' ');
end

if dispTerminalLabels
    set(h.terminalNodeLabels,'Fontsize',min(9,ceil(fontRatio/1.5)));
end

if ~dispBranchLabels
    set(h.branchNodeLabels,'visible','off');
end
if ~dispLeafLabels
    set(h.leafNodeLabels,'visible','off');
end
if ~dispTerminalLabels
    set(h.terminalNodeLabels,'visible','off');
end

%Reset x axis ticks
%trueDates = zeros(1, tr.numLeaves);
%for j = 1:tr.numLeaves
 %   leafName = tr.names(j);
 %   leafNameBlocks = regexp(leafName, ' ', 'split');
 %   leafDate = str2num(leafNameBlocks{1,1}{1,2}); %should be two if only one '_' in name
 %   trueDates(j) = leafDate;
%end
%trueRootTime = max(trueDates) - max(tr.x);
%trueRootTimeYears = trueRootTime / 365.25;
%xTickTrueTimes = 365.25*(ceil(trueRootTimeYears):1:20);
%xTickAxisLocs = xTickTrueTimes - trueRootTime;
%xTickAxisStart = ceil(trueRootTimeYears) - 16; % for root to 0 and 0 to 4
%xTickAxisLabels = (xTickAxisStart:1:4);

%For entropy trees where all time is negative
%xTickAxisStart = ceil(trueRootTimeYears) - 20; % for root to 0 and 0 to 4
%xTickAxisLabels = (xTickAxisStart:1:0);

%set(h.axes, 'XTick', xTickAxisLocs);
%set(h.axes, 'XTickLabel', xTickAxisLabels);

% Some additional plotting features
%handleNow = gca;
%keyboard
set(h.axes, 'Color', [0.6, 0.6, 0.6]) %set plot background to black
%set(gcf, 'color', 'k') %set figure background to black
hcbar = colorbar('EastOutside'); %add colorbar to plot
%set(hcbar, 'YTick', [0:5:20]); %set YTicks on colorbar
%set(hcbar, 'YTickLabel', [0:0.25:1.0]); %set labels
ylabel(hcbar, 'Mutation load', 'FontSize',12)
xlabel(h.axes, 'Time (years)', 'FontSize', 12); %set xlabel
set(h.axes, 'YTick', []); %turn off YTicks
set(h.axes,'Fontsize',12); %set axes font size
set(h.axes, 'XColor', [0.0, 0.0, 0.0]); %set to grey so visible with both white or black background

%legend('Patch 1','Patch 2'); legend('boxoff');

box off
hold off

% store handles
set(fig,'UserData',h)
if nargout
    handles = h;
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tr = doBasicCalculations(tr,activeBranches,renderType)

% helper function to compute and find some features of the tree
tr.numLeaves = tr.numBranches + 1;
tr.numLabels = tr.numBranches + tr.numLeaves;

% remove uderscores from names
for ind = 1:tr.numLabels
    tr.names{ind}(tr.names{ind}=='_')=' ';
end

% obtain parents for every node
tr.par(tr.tree(:)) = tr.numLeaves + [1:tr.numBranches 1:tr.numBranches];

% find active nodes
tr.activeNodes = true(tr.numLabels,1);
for ind =tr.numBranches:-1:1
    tr.activeNodes(tr.tree(ind,:)) = tr.activeNodes(tr.numLeaves+ind) & activeBranches(ind);
end

% propagate last leaf
tr.lastleaf = 1:tr.numLabels;
for ind = tr.numBranches:-1:1
    if ~tr.activeNodes(tr.tree(ind,1))
        tr.lastleaf(tr.tree(ind,:))=tr.lastleaf(ind+tr.numLeaves);
    end
end

tr.activeBranches = tr.activeNodes(tr.numLeaves+1:tr.numLabels)&activeBranches;
tr.activeLeaves = tr.activeNodes(1:tr.numLeaves);

% find x coordinates of branches
tr.x = tr.dist;
for ind = tr.numBranches:-1:1
    tr.x(tr.tree(ind,:)) = tr.x(tr.tree(ind,:)) + tr.x(ind+tr.numLeaves);
end

% find y coordinates of branches
tr.terminalNodes = tr.lastleaf([true,diff(tr.lastleaf(1:tr.numLeaves))~=0]);
tr.y=zeros(tr.numLabels,1);
tr.y(tr.terminalNodes)=1:length(tr.terminalNodes);
switch renderType
    case 'square'
        for ind = 1:tr.numBranches
            if tr.activeBranches(ind)
                tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
            end
        end
    case {'angular','radial'}
        for ind = 1:tr.numBranches
            if tr.activeBranches(ind)
                if tr.x(tr.tree(ind,1))/tr.x(tr.tree(ind,2))>3
                    tr.y(ind+tr.numLeaves) = tr.y(tr.tree(ind,1));
                elseif tr.x(tr.tree(ind,2))/tr.x(tr.tree(ind,1))>3
                    tr.y(ind+tr.numLeaves) = tr.y(tr.tree(ind,2));
                else
                    tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
                end
            end
        end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function correctFigureSize(fig,recommendedHeight,recommendedWidth)
% helper function to increase initial figure size depending on the screen &
% tree sizes
screenSize = diff(reshape(get(0,'ScreenSize'),2,2),[],2)-[0;100];
% 100 gives extra space for the figure header and win toolbar
position = get(fig,'Position');
if recommendedHeight > position(4)
    if recommendedHeight < sum(position([2 4]))
        position(2) = sum(position([2 4])) - recommendedHeight;
        position(4) = recommendedHeight;
    elseif recommendedHeight < screenSize(2)
        position(2) = 30;
        position(4) = recommendedHeight;
    else
        position(2) = 30;
        position(4) = screenSize(2);
    end
end
if recommendedWidth > position(3)
    if recommendedWidth < sum(position([1 3]))
        position(1) = sum(position([1 3])) - recommendedWidth;
        position(3) = recommendedWidth;
    elseif recommendedWidth < screenSize(1)
        position(1) = 0;
        position(3) = recommendedHeight;
    else
        position(1) = 0;
        position(3) = screenSize(1);
    end
end
set(fig,'Position',position)

