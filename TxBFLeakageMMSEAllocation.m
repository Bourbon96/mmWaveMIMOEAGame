%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Wenbo Wang
%
% Original Article:
% [Wang2021] Wenbo Wang and Amir Leshem, "Decentralized Power Allocation
% and Beamforming Using Non-Convex Nash Game for Energy-Aware mmWave Networks"
%
% References:
% [1] F. Sun and E. de Carvalho, "A Leakage-Based MMSE Beamforming Design for a MIMO 
%     Interference Channel," in IEEE Signal Processing Letters, vol. 19, no. 6, pp. 
%     368-371, June 2012
%
% License: 
% This program is licensed under the GPLv3 license. If you in any way use this code for 
% research that results in publications, please cite our original article listed above.
 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU General Public License for more details.

function [TxBFLocal, PowerLocal] = TxBFLeakageMMSEAllocation(LinkID, InitRxBFs, InitTxBFs, InitPower, SystemSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% (1) The Tx-BF design is a simplified version of Eq. (31), with the hard 
% constraints on the interference leakage relaxed to soft penalties.
% (2) The objective function (cost func.) is then reduced to Eq. (5) of
% Ref. [1]. By doing so, we avoid the process of handling Lagrange
% multipliers. However, we are not able to update the local power precisely
% with this scheme.
%
% Input: 
%        InitRxBFs         : Intial Rx-beamformer vectors
%        InitTxBFs         : Intial Tx-beamformer vectors (not used)
%        InitPower         : Intial power values (single stream)
%        SystemSettings    : a struct of system settings, e.g., MIMO
%                            channels,SINR thresolds, maximum power values
%
% Output: 
%         TxBFLocal        : Tx-beamforing vectors, Nt*1
%         PowerLocal       : Local power strategy (changed due to global Tx MMSE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nlink = SystemSettings.Nlink;
MIMOChannels = SystemSettings.MIMOChannels;
Nt = SystemSettings.Nt;

% Solve the power-constrained optimization problem (closed-form solution)
H_kk = MIMOChannels{LinkID, LinkID};
P_kk = InitPower(LinkID, 1);
RxBF_k = InitRxBFs{LinkID, 1};

TMatrix = zeros(Nt, Nt);
for jj=1:Nlink
    RxBF_jj = InitRxBFs{jj, 1};
    H_kj = MIMOChannels{LinkID, jj};
    TMatrix = TMatrix + P_kk*H_kj'*(RxBF_jj*RxBF_jj')*H_kj;
end

alpha_k = trace(H_kk*H_kk')/P_kk;

TxBFLocalTmp = inv(TMatrix + alpha_k*eye(Nt))*H_kk'*RxBF_k; % Eq.(32) of Ref. [1]
PowerLocal = trace(TxBFLocalTmp'*TxBFLocalTmp);

TxBFLocal = TxBFLocalTmp/sqrt(trace(TxBFLocalTmp'*TxBFLocalTmp)); %norm(TxBFLocal)

end

