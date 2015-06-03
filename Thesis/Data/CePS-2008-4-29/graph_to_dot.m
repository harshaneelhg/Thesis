function graph_to_dot(adj, varargin)

% graph_to_dot(adj, VARARGIN)  Creates a GraphViz (AT&T) format file representing 
%                     a graph given by an adjacency matrix.
%  Optional arguments should be passed as name/value pairs [default]
%
%   'filename'  -  if omitted, writes to 'tmp.dot'
%  'arc_label'  -  arc_label{i,j} is a string attached to the i-j arc [""]
% 'node_label'  -  node_label{i} is a string attached to the node i ["i"]
%  'width'      -  width in inches [10]
%  'height'     -  height in inches [10]
%  'leftright'  -  1 means layout left-to-right, 0 means top-to-bottom [0]
%  'directed'   -  1 means use directed arcs, 0 means undirected [1]
%
% For details on dotty, See http://www.research.att.com/sw/tools/graphviz
%
% by Dr. Leon Peshkin, Jan 2004      inspired by Kevin Murphy's  BNT
%    pesha @ ai.mit.edu /~pesha
                   
node_label = [];   arc_label = [];   % set default args
width = 10;        height = 10;
leftright = 0;     directed = 1;     filename = 'tmp.dot';
[m,n]=size(adj);
isQuery = zeros(1,m);   
Rank = [];
for i = 1:2:nargin-1              % get optional args
    switch varargin{i}
        case 'filename', filename = varargin{i+1};
        case 'node_label', node_label = varargin{i+1};
        case 'arc_label', arc_label = varargin{i+1};
        case 'width', width = varargin{i+1};
        case 'height', height = varargin{i+1};
        case 'leftright', leftright = varargin{i+1};
        case 'directed', directed = varargin{i+1};
        case 'isQuery', isQuery = varargin{i+1};
        case 'setRank', Rank = varargin{i+1}; 
    end
end
fid = fopen(filename, 'w');
if directed
    fprintf(fid, 'digraph G {\n');
    arctxt = '->'; 
    if isempty(arc_label)
        labeltxt = '';
    else
        labeltxt = '[label="%s"]';
    end
else
    fprintf(fid, 'graph G {\n');
    arctxt = '--'; 
    if isempty(arc_label)
        labeltxt = '[dir=none]';
    else
        labeltxt = '[label="%s",dir=none]';
    end
end
fprintf(fid, 'center = 1;\n');
fprintf(fid, 'size=\"%d,%d\";\n', width, height);
if leftright
    fprintf(fid, 'rankdir=LR;\n');
end
pos = find(isQuery>0);%with two query
if length(pos)==2
    Rank = isQuery;
    Rank(pos(1)) = 1;
    Rank(pos(2)) = 2;
end
Nnds = length(adj);
if isempty(Rank)
else%set the rank of the node
    for node = 1:Nnds
        if Rank(node)==0%not set
        elseif Rank(node)==1%Set source
            fprintf(fid,'{rank = source;\n "%d";};\n',node);
        elseif Rank(node)==2%Set Sink
            fprintf(fid,'{rank = sink;\n "%d";};\n ',node);
        end
    end
end
%fprintf(fid, 

for node = 1:Nnds               % process NODEs 
    if isempty(node_label)
        fprintf(fid, '%d;\n', node);
    elseif isQuery(node)==0
        fprintf(fid, '%d [ label = "%s"];\n', node, node_label{node});
    else
        fprintf(fid, '%d [ label = "%s", shape = "%s", style ="%s", color ="%s" ];\n', node, node_label{node},'diamond','filled','skyblue');
        %if isempty(Rank)
        %    fprintf(fid, '%d [ label = "%s", shape = "%s", style ="%s", color ="%s" ];\n', node, node_label{node},'diamond','filled','skyblue');
        %elseif Rank(node)==0%not set
        %    fprintf(fid, '%d [ label = "%s", shape = "%s", style ="%s", color ="%s" ];\n', node, node_label{node},'diamond','filled','skyblue');
        %elseif Rank(node)==1%Set source
        %    fprintf(fid, '%d [ rank="source", label = "%s", shape = "%s", style ="%s", color ="%s"];\n', node, node_label{node},'diamond','filled','skyblue');
        %elseif Rank(node)==2%Set sink
        %    fprintf(fid, '%d [ rank="sink", label = "%s", shape = "%s", style ="%s", color ="%s" ];\n', node, node_label{node},'diamond','filled','skyblue');
        %end
    end
end
edgeformat = strcat(['%d ',arctxt,' %d ',labeltxt,';\n']);
for node1 = 1:Nnds              % process ARCs
    if directed
        arcs = find(adj(node1,:));         % children(adj, node);
    else
        arcs = find(adj(node1,node1+1:Nnds)) + node1; % remove duplicate arcs
    end
    for node2 = arcs
        if isempty(arc_label)     % thanks to Nicholas Wayne Henderson nwh@owlnet.rice.edu
            fprintf(fid, edgeformat, node1, node2);  
        else
            fprintf(fid, edgeformat, node1, node2, arc_label{node1,node2});
        end
    end
end
fprintf(fid, '}'); 
fclose(fid);