%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Wenbo Wang
%
% Original Article:
% [Wang2021] Wenbo Wang and Amir Leshem, "Decentralized Power Allocation
% and Beamforming Using Non-Convex Nash Game for Energy-Aware mmWave Networks"
%
% License: 
% This program is licensed under the GPLv3 license. If you in any way use this code for 
% research that results in publications, please cite our original article listed above.
 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU General Public License for more details.
%
function RxBF = RxBFMVDRAllocation(LinkID, InitTxBFs, InitPower, SystemSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
%        LinkID            : ID of the link of interest
%        InitTxBFs         : Intial Tx-beamformer vectors
%        InitPower         : Intial power values (single stream)
%        SystemSettings    : a struct of system settings, e.g., MIMO
%                            channels,SINR thresolds, maximum power values
%
% Output: 
%         RxBF            : Rx-beamforing vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MIMOChannels = SystemSettings.MIMOChannels;
sigma2 = SystemSettings.sigma2;

Nlink = SystemSettings.Nlink;
Nr = SystemSettings.Nr;

%% One step update for beamformers (Lines 5,6 of Alg.2 in [Wang2021])
H_ii = MIMOChannels{LinkID,LinkID};
TxBF_ii = InitTxBFs{LinkID, 1};
Interference = sigma2*eye(Nr);
for jj=1:Nlink
    if jj == LinkID
        continue;
    end
    H_ji = MIMOChannels{jj,LinkID}; %interference channel
    Power_jj = InitPower(jj, 1);
    TxBF_jj = InitTxBFs{jj, 1};
    Interference = Interference + Power_jj*H_ji*(TxBF_jj*TxBF_jj')*H_ji';
end
% Rx-beamforming with MMSE
RxBF = sqrt(InitPower(LinkID))*inv(Interference)*H_ii* TxBF_ii;

end