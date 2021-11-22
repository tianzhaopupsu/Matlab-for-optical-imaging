function h = pl3D(filename, varargin)
close all
% pl3D('/Users/TianZhao/Desktop/D21/tracking21-02-21_225004')
% color_line3 plots a 3-D "line" with c-data as color
%
%       h = color_line(x, y, z, c)
%       by default: 'LineStyle','-' and 'Marker','none'
%
%          or
%       h = color_line(x, y, z, c, mark) 
%          or
%       h = color_line(x, y, z, c, 'Property','value'...) 
%             with valid 'Property','value' pairs for a surface object
%
%  in:  x      x-data
%       y      y-data
%       z      z-data
%       c      4th dimension for colouring
%       mark   for scatter plots with no connecting line
%
% out:  h   handle of the surface object
fid = fopen(filename);
chan1 = fread(fid,[11, inf],'double');
 chan1=chan1';
trackDataFilt=zeros(length(chan1),11);
trackDataFilt(:,1)=chan1(:,1);
trackDataFilt(:,2)=chan1(:,2);
trackDataFilt(:,3)=chan1(:,3);
trackDataFilt(:,4)=chan1(:,4);
trackDataFilt(:,5)=chan1(:,5);
trackDataFilt(:,6)=chan1(:,6);
trackDataFilt(:,7)=chan1(:,7);
trackDataFilt(:,8)=chan1(:,8);
trackDataFilt(:,9)=chan1(:,9);
trackDataFilt(:,10)=chan1(:,10);
trackDataFilt(:,11)=chan1(:,11);
a=length(trackDataFilt);

trackDataFilt=trackDataFilt(0.05*a:0.95*a,1:11);
intensity=trackDataFilt(:,7)+trackDataFilt(:,8)+trackDataFilt(:,9)+trackDataFilt(:,10)+trackDataFilt(:,11);

b=length(trackDataFilt);
c=1:length(trackDataFilt);
x=trackDataFilt(:,1);
y=trackDataFilt(:,2);
z=trackDataFilt(:,3);
x=x-x(1,1);
y=y-y(1,1);
z=z-z(1,1);
intensity=(intensity/2.4);
x=x;
y=y*1.5;
%z=z*2;



surface(...
  'XData',[x(:) x(:)],...
  'YData',[y(:) y(:)],...
  'ZData',[z(:) z(:)],...
  'CData',[c(:) c(:)],...
  'FaceColor','none',...
  'EdgeColor','flat',...
  'Marker','none','Linewidth', 0.1);
set(gca,'fontsize',15)
     xlabel('X coordinates (\mum)','FontSize',15)
        ylabel('Y coordinates (\mum)','FontSize',15)
        zlabel('Z coordinates (\mum)','FontSize',15)
        c = colorbar;
    
          c.TickLabels = {''};
        c.Location = 'southoutside';
       set(gca,'fontsize',15);
c.Label.String = 'Trajectory length = 12.1 s';
 grid on
set(gca,'GridLineStyle','--')        
figure;
        plot(intensity);
        set(gca,'fontsize',15)
        xlabel('Time (ms)','FontSize',15)
        ylabel('Total Photon Counts (kcps)','FontSize',15)

if nargin ==5
    switch varargin{1}
        case {'+' 'o' '*' '.' 'x' 'square' 'diamond' 'v' '^' '>' '<' 'pentagram' 'p' 'hexagram' 'h'}
            set(h,'LineStyle','none','Marker',varargin{1})
        otherwise
            error(['Invalid marker: ' varargin{1}])
    end
elseif nargin > 5
    set(h,varargin{:})
end
save('Intensity')
display(num2str(b/1000))