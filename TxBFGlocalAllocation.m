%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Wenbo Wang
%
% Original Article:
% [Wang2021] Wenbo Wang and Amir Leshem, "Decentralized Power Allocation
% and Beamforming Using Non-Convex Nash Game for Energy-Aware mmWave Networks"
%
% References:
% [1] H. Shen, B. Li, M. Tao and X. Wang, "MSE-Based Transceiver Designs for the MIMO 
%     Interference Channel," in IEEE Transactions on Wireless Communications, vol. 9, 
%     no. 11, pp. 3480-3489, November 2010.
%
% License: 
% This program is licensed under the GPLv3 license. If you in any way use this code for 
% research that results in publications, please cite our original article listed above.
 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU General Public License for more details.

function [TxBFLocal, PowerLocal] = TxBFGlocalAllocation(LinkID, InitRxBFs, InitTxBFs, InitPower, SystemSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% (1) Return a Tx beamformer for link "LinkID", based on the global
% information, to maximize the constrained sum of MSE in a cooperative manner
% Constraint: ||w_i||^2 = 1.
% (2) This alogrithm is replaced by , and the GNE solution for power is not
% used for sumb-MSE-based Tx beamforming.
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

% Solve the regularized optimization problem (closed-form solution)
H_kk = MIMOChannels{LinkID, LinkID};
P_kk = InitPower(LinkID, 1);
RxBF_k = InitRxBFs{LinkID, 1};

TMatrix = zeros(Nt, Nt);
for jj=1:Nlink
    RxBF_jj = InitRxBFs{jj, 1};
    H_jk = MIMOChannels{jj, LinkID};
    TMatrix = TMatrix + P_kk * H_jk'*(RxBF_jj*RxBF_jj')*H_jk;
end

% SVD of TMatrix
[V_k, Delta_k] = eig(TMatrix);

% Solve a equation to obtain the proper Lagrange multiplier
G_k = V_k'*H_kk'*(RxBF_k*RxBF_k')*H_kk*V_k;

% Nested function to be solved numerically for Lagrange multiplier
    function DiagonalElementSum = PowerConstraint(lambda)
        DiagonalElementSum = -1;
        for ii=1:Nt
            DiagonalElementSum = DiagonalElementSum + G_k(ii, ii)/(lambda + Delta_k(ii, ii))^2;
        end
    end

% The global optimal-MSE transmitter
LambdaCandidates = fsolve(@PowerConstraint, 1, optimset('Display','off'));
Lambda = max(LambdaCandidates(max(imag(LambdaCandidates))<=1e-10 & real(LambdaCandidates)>0));
if size(Lambda, 1) == 0
    % In case the numerical solver cannot find a proper root
    TxBFLocal = sqrt(P_kk)*inv(TMatrix+ 1*eye(Nt))*H_kk'*RxBF_k;
else
    TxBFLocal = sqrt(P_kk)*inv(TMatrix + Lambda(1)*eye(Nt))*H_kk'*RxBF_k;
end

TxBFLocal = TxBFLocal/sqrt(TxBFLocal'*TxBFLocal);

PowerLocal = P_kk;
end

