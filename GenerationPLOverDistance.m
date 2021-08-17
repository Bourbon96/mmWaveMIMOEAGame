%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Wenbo Wang
%
% Original Article:
% [Wang2021] Wenbo Wang and Amir Leshem, "Decentralized Power Allocation
% and Beamforming Using Non-Convex Nash Game for Energy-Aware mmWave Networks"
%
% Additional Reference:
% [1] A. I. Sulyman, A. Alwarafy, G. R. MacCartney, T. S. Rappaport and A. Alsanie, 
%     "Directional Radio Propagation Path Loss Models for Millimeter-Wave Wireless Networks 
%     in the 28-, 60-, and 73-GHz Bands," in IEEE Transactions on Wireless Communications, 
%     vol. 15, no. 10, pp. 6939-6947, Oct. 2016.
%
% License: 
% This program is licensed under the GPLv3 license. If you in any way use this code for 
% research that results in publications, please cite our original article listed above.
 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU General Public License for more details.
%
function [DistanceMatrix, PLMatrix] = GenerationPLOverDistance(SystemSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
%        SystemSettings    : a struct of system settings, e.g., link
%                            number, Rx/Tx antenna numbers, etc
%
% Output: 
%         DistanceMatrix   : Matrix of Tx-Rx distances
%         PLMatrix         : Pathloss matrix, Nlink*Nlink matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Nlink = SystemSettings.Nlink;

width = 40; % hard coded (a small cell)
widthSeperation = width*3.5;
% arrange the area in a number os sub-rectangles, every 5 links per
% sub-rectangles
NumLinkPerArea = 1;
NumArea = ceil(Nlink/NumLinkPerArea);
LinkList = randperm(Nlink);
NGrid = ceil(sqrt(NumArea));

TxPos = zeros(Nlink, 2);
RxPos = zeros(Nlink, 2);
for ii=1:NGrid
    xStart = (ii-1)*widthSeperation;
    for jj=1:NGrid
        IdArea = NGrid*(ii-1) + jj;
        yStart = (jj-1)*widthSeperation;
        if IdArea <= NumArea
            LinkIDRange = (NumLinkPerArea*(IdArea-1)+1:min(NumLinkPerArea*IdArea, Nlink))';    
            szLinkRange = size(LinkIDRange, 1);
            TxPos(LinkList(LinkIDRange), 1) = width * rand(1, szLinkRange) - width / 2 + xStart;
            TxPos(LinkList(LinkIDRange), 2) = width * rand(1, szLinkRange) - width / 2 + yStart;
            
            tmpAngle = 2*pi*rand(1,szLinkRange);
            tmpR = width*sqrt(rand(1,szLinkRange)) + width/4;
            RxPos(LinkList(LinkIDRange), 1) = TxPos(LinkList(LinkIDRange), 1) + tmpR.*cos(tmpAngle);
            RxPos(LinkList(LinkIDRange), 2) = TxPos(LinkList(LinkIDRange), 2) + tmpR.*sin(tmpAngle);
            
%             RxPos(LinkList(LinkIDRange), 1) = width * rand(1, szLinkRange)  - width / 2 + TxPos(LinkList(LinkIDRange), 1);
%             RxPos(LinkList(LinkIDRange), 2) = width * rand(1, szLinkRange)  - width / 2 + TxPos(LinkList(LinkIDRange), 2);
        end    
    end
end

DistanceMatrix = zeros(Nlink, Nlink);
for ii=1:Nlink
    for jj=1:Nlink
        DistanceMatrix(ii,jj) = norm(TxPos(ii,:) - RxPos(jj,:));
    end
end

%% For simple Log-distance PL model (d0=1), may not be used to generate path loss
gamma_pl = 3.0;
fc = 28e9;
fspl = 32.4 + 20 * log10(fc / 1e9); % Eq.(3) of Ref. [1]
PLMatrix = zeros(Nlink, Nlink);
for ii=1:Nlink
    for jj=1:Nlink
        dist = DistanceMatrix(ii,jj);      
        PLdB = fspl + 10 * gamma_pl * log10(dist);  
        PL = 10^(-PLdB/20); % Debugging
        PLMatrix(ii,jj) = PL;
    end
end

% figure;
% title("Topology of the Tx-Rx pairs");
% hold on;
% plot(TxPos(:,1), TxPos(:,2), 'r.', 'MarkerSize', 10);
% plot(RxPos(:,1), RxPos(:,2), 'b*', 'MarkerSize', 10);

end