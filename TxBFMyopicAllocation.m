%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Wenbo Wang
%
% Original Article:
% [Wang2021] Wenbo Wang and Amir Leshem, "Decentralized Power Allocation
%  and Beamforming Using Non-Convex Nash Game for Energy-Aware mmWave Networks"
%
% Reference:
% [1] R. A. Iltis, Seung-Jun Kim and D. A. Hoang, "Noncooperative iterative MMSE beamforming 
%  algorithms for ad hoc networks," in IEEE Transactions on Communications, vol. 54, no. 4, 
%  pp. 748-759, April 2006,
%
% License: 
% This program is licensed under the GPLv3 license. If you in any way use this code for 
% research that results in publications, please cite our original article listed above.
 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU General Public License for more details.

function [TxBFLocal, PowerLocal] = TxBFMyopicAllocation(LinkID, InitRxBFs, InitTxBFs, InitPower, SystemSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% Return a greedy Tx beamformer for link "LinkID", solely based on the local
% information, to maximize the SINR in a noncooperative manner, such that the
% local power can be further reduced
%
% Input: 
%        InitRxBFs         : Intial Rx-beamformer vectors (not used)
%        InitTxBFs         : Intial Tx-beamformer vectors
%        InitPower         : Intial power values (single stream)
%        SystemSettings    : a struct of system settings, e.g., MIMO
%                            channels,SINR thresolds, maximum power values
%
% Output: 
%         TxBFLocal        : Tx-beamforing vectors, Nt*1
%         PowerLocal       : Local power strategy (unchanged)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nlink = SystemSettings.Nlink;
MIMOChannels = SystemSettings.MIMOChannels;
sigma2 = SystemSettings.sigma2;
Nr = SystemSettings.Nr;
SINRThreshold = SystemSettings.SINRThresholds(LinkID);

Hlocal = MIMOChannels{LinkID, LinkID};
Interference = sigma2*eye(Nr); % Interference plus noise Eq.(4) of [1]
for jj=1:Nlink
    if jj ~= LinkID
        H_ji = MIMOChannels{jj,LinkID};
        Power_jj = InitPower(jj, 1);
        TxBF_jj = InitTxBFs{jj, 1};
        Interference = Interference + Power_jj*H_ji*(TxBF_jj*TxBF_jj')*H_ji';
    end
end

% greedy noncooperative SNR maximization beamforming game
Htilde = Hlocal'*inv(Interference)*Hlocal;
%% To maximize w^H*Htilde^{-1}*w
[TxBFLocal, ~] = eigs(Htilde, 1); % Eq. (11) of [1]

% Eq.(11) of [1]
PowerLocal = SINRThreshold/real(TxBFLocal'*Htilde*TxBFLocal);

end