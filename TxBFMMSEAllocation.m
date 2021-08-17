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

function [TxBFLocal, PowerLocal] = TxBFMMSEAllocation(LinkID, InitRxBFs, InitTxBFs, InitPower, SystemSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% Return a greedy MMSE beamformer for link "LinkID", soley based on the local
% information, to maximize the SINR in a noncooperative manner, such that the
% local power can be further reduced
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
%         PowerLocal       : Local power strategy (unchanged)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LocalRxBF = InitRxBFs{LinkID, 1};
Hlocal = SystemSettings.MIMOChannels{LinkID, LinkID};
PowerLocal = InitPower(LinkID, 1);
Nt = SystemSettings.Nt;

% Eigenvalue decomposition, see Eq. (23)
T_ii = PowerLocal*Hlocal'*(LocalRxBF*LocalRxBF')*Hlocal;
[V, M] = eig(T_ii);
G_ii = V'*Hlocal'*(LocalRxBF*LocalRxBF')*Hlocal*V;

%% Nested function to be solved numerically
    function DiagonalElementSum = PowerConstraint(lambda)
        DiagonalElementSum = -1;
        for ii=1:Nt
            DiagonalElementSum = DiagonalElementSum + real(G_ii(ii, ii))/(lambda + real(M(ii, ii)))^2;
        end
    end

LambdaCandidates = fsolve(@PowerConstraint, 1, optimset('Display','off'));
Lambda = max(LambdaCandidates(max(imag(LambdaCandidates))<=1e-10 & real(LambdaCandidates)>0));

% Eq.(27), [Wang2021]
if size(Lambda, 1) == 0
    % In case the numerical solver cannot find a proper root
    TxBFLocal = sqrt(PowerLocal)*inv(T_ii+ 1*eye(Nt))*Hlocal'*LocalRxBF;
else
    TxBFLocal = sqrt(PowerLocal)*inv(T_ii + Lambda(1)*eye(Nt))*Hlocal'*LocalRxBF;
end

TxBFLocal = TxBFLocal / sqrt(real(TxBFLocal'*TxBFLocal));

end