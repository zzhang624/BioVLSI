function ACLfig(type,opt)
switch lower(type)
    case 'fullslide'; ifig = 1;
    case 'halfslide'; ifig = 2;
    case '4slide';    ifig = 3;
    case 'jpaper';    ifig = 4;
    otherwise;  
        msgbox(sprintf('Unknown figure type ''%s'', using ''fullslide''',type),'ACLfig');
        ifig = 1;
end

% set varous figure properties
% *** add more columns for more types ***
% type =    full  half  4sli  jpap
% ifig =     1     2     3     4
p.figW =   [8     4.55  3.5   3.25]; % figure width
p.figH =   [6     3.35  2.25  2.70]; % figure height
p.plotLW = [1.4   1.4   1.4   1.0 ]; % plotted line linewidth
p.axisLW = [1.4   1.25  1.25  0.9 ]; % axis linewidth
p.tickF =  [ 14    12    11     9 ]; % tick label font size
p.axisF =  [ 16    13    11     9 ]; % axis title font size

% check for opt
if nargin > 1 && ~isempty(fieldnames(opt))
   fNames = fieldnames(opt);
   for ii = 1:length(fNames)
       p.(fNames{ii})(ifig) = opt.(fNames{ii});
   end
end

% remove most of whitespace produced by default
set(gca,'LooseInset',get(gca,'TightInset')+(0.02));

% set figure size (and screen placement, approximately centered)
oldUnits = get(0,'Units'); % get curent units
set(0,'Units','inches'); % set units to inches
screenSize = get(0,'ScreenSize'); % get screen size
set(0,'Units',oldUnits); % revert units to previous setting
left = screenSize(3)/2 - p.figW(ifig)/2; % left edge of fig
bottom = screenSize(4)/2 - p.figH(ifig)/2; % bottom edge of fig

set(gcf,'Color','w',...
        'Units','inches',...
        'Position',     [left bottom p.figW(ifig) p.figH(ifig)],...
        'PaperSize',                [p.figW(ifig) p.figH(ifig)],...
        'PaperPosition',[0    0      p.figW(ifig) p.figH(ifig)])

% set font and font sizes
set(gca,'FontName','Arial')
set(gca,'FontSize',p.tickF(ifig))               % tick labels, legend labels...
set(get(gca,'xlabel'),'FontSize',p.axisF(ifig)) % xlabel
set(get(gca,'ylabel'),'FontSize',p.axisF(ifig)) % ylabel
set(get(gca,'title'),'FontSize',p.axisF(ifig))  % title

% set linewidths
%set(get(gca,'Children'),'LineWidth',p.plotLW(ifig)) % plot linewidth
set(gca,'LineWidth',p.axisLW(ifig))                 % axis linewidth
