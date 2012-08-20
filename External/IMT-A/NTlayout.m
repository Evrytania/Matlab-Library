% NTLAYOUT Plot of the network layout in cartesian coordinate system
%   NTLAYOUT(LINKPAR) draws the locations of Ms/Bs/Sector and the array orientations  
%   in the (x,y) coordinate system including also the randomly generated active links
%   See [1, Fig. 6.1 and 6.2].
%
%   NTLAYOUT(LINKPAR) is used only when the layout parameters are defined. 
%   It is used the following way: 
%   NTlayout(layout2link(layoutparset(NoMs,NofBs,SectPerBs,K)))
%
%   Ref. [1]: D5.4, "Final Report on Link Level and System Level Channel Models"
%
%   See also LAYOUTPARSET, LINKPARSET.

% Authors: Daniela Laselva (EBIT)

function NTlayout(linkpar)

BsXY=linkpar.BsXY;           % coordinates of BSs
MsXY=linkpar.MsXY;           % coordinates of BSs
ThetaMs=linkpar.ThetaMs;      
ThetaBs=linkpar.ThetaBs;
NofSect=linkpar.NofSect;
MsOmega=linkpar.MsOmega;
BsOmega=linkpar.BsOmega;
Pairing=linkpar.Pairing;
MsDirection=linkpar.MsDirection;
NofMs=size(MsXY,2);          % number of Ms
NofBs=size(BsXY,2);          % number of Bs

flag=1; % if 1 draw Bs broadside, if 2 active links, if 3 Ms'directions

% Get number of links and Ms & Sector indices            
K = sum(Pairing(:));    % number of links
[r c] = find(Pairing);
indMs = c';             % links coupling to MSs
indSect = r';           % links coupling to Sectors (1:NofSect*NofBs)
% Get links coupling to BSs
tmpsum = cumsum(NofSect);  
for k=1:length(r) tmp = find(diff(r(k)<=tmpsum));
    if isempty(tmp) indBs(k)=1; else indBs(k)=tmp+1;end
end
indSectors=mod(reshape(indSect,1,length(indSect)),NofSect(indBs));  % links coupling to Sect (1:NofSect)
indSectors(find(indSectors==0))=NofSect(indBs(find(indSectors==0)));

% Draw locations of BSs and MSs arrays
clf;hold on; grid on 
for i=1:NofMs, locationMs(MsXY(:,i)); end           % Draw MSs arrays
hold on; grid on
for i=1:NofBs, locationBs(BsXY(:,i)); end           % Draw BSs arrays
plot(BsXY(1,:),BsXY(2,:), 'rO','linewidth',3);
xlabel('Cells area (m)','FontSize',12); ylabel('Cells area (m)','FontSize',12);
title('NETWORK LAYOUT: BSs and BS Sectors, MSs and MS Directions, Active Links',...
        'FontSize',10)

% Draw Bs arrays broadside for each sector 
flag=1;  % draw Bs broadside
for i=1:NofBs n=NofSect(i); for nn=1:n arrow(BsXY,MsXY,indBs,indMs,BsOmega,nn,i,flag);...
        end,end

% Draws Bs-Ms active links related to the correspondent Bs sector
flag=2;  % draw active links
for i=1:K arrow(BsXY,MsXY,indBs,indMs,BsOmega,indSectors(i),indBs(i),flag,i),end

% Draw Ms' directions 
flag=3;  % draw Ms'directions
if K<NofMs MsDirection=[MsDirection, 360*(rand(1,NofMs-K)-0.5)]; end % Directions not available for each Ms
for i=1:NofMs arrow(MsXY,BsXY,indMs,indBs,MsDirection,1,i,flag), end
hold off
% End of NTlayout



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that draws arrows' lines 
function arrow(XY1,XY2,ind1,ind2,Omega,nn,i,flag,nlink);

% const
l1    = 10;   % scale factor
l2    = 1;    % scale factor
l3    = 2;    % scale factor

% Convention for Omega depends on the position Bs vs Ms (D5.4,Fig.6.2)
link=min(find(ind1==i));
if link                        % link=1 if a link exists 
    if XY1(1,i)<XY2(1,ind2(link)) OmegaTmp=-Omega(nn,i); else OmegaTmp=Omega(nn,i); end
else OmegaTmp=Omega(nn,i); end % Ms not linked to any Bs, or Bs not covering any link

% Calculate verteces of the arrows v1 &  v2 
sen=sin(OmegaTmp*pi/180); 
cosen=cos(OmegaTmp*pi/180);
v1 = [XY1(1,i)-l1/2*sen XY1(2,i)+l1/2*cosen];                 % vertex1 of the arrow
if flag==2 v2 = [XY2(1,ind2(nlink)) XY2(2,ind2(nlink))+l1/5]; % vertex2 of the arrow
else v2 = [XY1(1,i)+l1/2*sen XY1(2,i)-l1/2*cosen];            % vertex2 of the arrow
end

% Plot arrows: links, Ms' directions and Bs' arrays
if flag==2                                    % Plot active links 
    line([XY2(1,ind2(nlink)) XY1(1,i)-l1/2*sen],...
        [XY2(2,ind2(nlink))+2 XY1(2,i)+l1/2*cosen],'Color','b','LineWidth',2)   
    [xh, yh] = harrow(v1,v2,l1);                    % arrow head for Bs
    plot(xh,yh,'b','linewidth',2);                  % Plot arrow head for links'
    fill(xh, yh,'b')
else                                          % flag = 1 or 3
    [xh, yh] = harrow(v1,v2,l2);                    % arrow head for Bs or Ms
    if flag==1;                               % Plot Bs broadside
        line([XY1(1,i) v1(1)], [XY1(2,i) v1(2)],'Color', 'r','LineWidth',2); hold on;
        plot(xh(:),yh(:),'r','linewidth',3);        % Plot arrow head for Bs 
    else                                      % Plot Ms'direction
        line([v2(1) v1(1)], [v2(2)-2*l3 v1(2)-2*l3],'Color', 'k','LineStyle','--',...
            'LineWidth',1); hold on;
        plot(xh(:),yh(:)-2*l3,'k','linewidth',2);   % Plot arrow head for Ms'direction
    end 
end
% End function arrow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that calculates the arrow head points x-head & y-head
function [xh,yh]=harrow(v1,v2,l);

hsiz = 0.2;   % Size arrow head vs length of the vector
hwid  = 0.5;  % Width of the arrow head base vs length

v  = (v1-v2); 
xh = [v1(1)-hsiz*(v(1)/l+hwid*v(2)/l); v1(1); v1(1)-hsiz*(v(1)/l-hwid*v(2)/l)]; 
yh = [v1(2)-hsiz*(v(2)/l-hwid*v(1)/l); v1(2); v1(2)-hsiz*(v(2)/l+hwid*v(1)/l)]; 
% End function arrowhead


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that draws the location of Ms
function locationMs(XY);

l = 2;  % scale factor
a = [XY(1)-l; XY(2)-2*l]; % vertex 1
b = [XY(1); XY(2)-2*l];   % vertex 2
c = [XY(1); XY(2)];       % vertex 3
d = [XY(1)-l; XY(2)];     % vertex 4
x1 = [ a(1), b(1), c(1), d(1) ];
y1 = [ a(2), b(2), c(2), d(2) ];
fill ( x1, y1, 'c' );     % Ms as rectangular shape
line([c(1)-0.05 c(1)-0.05], [c(2) c(2)+l],'Color', 'k','LineWidth',3 ) % plot drawing of the MS'antenna
% End function locationMs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that draws the location of Bs
function locationBs(XY);

l = 5;  % scale factor
a = [XY(1)-1/3*l; XY(2)-2*l];  % vertex 1
b = [XY(1)+1/3*l; XY(2)-2*l];  % vertex 2  
x1 = [ XY(1), a(1), b(1) ];
y1 = [ XY(2), a(2), b(2) ];
fill ( x1, y1, 'r' );          % BS as triangular shape
% End function locationBs