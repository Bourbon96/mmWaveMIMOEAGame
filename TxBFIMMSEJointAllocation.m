%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Wenbo Wang
%
% Original Article:
% [Wang2021] Wenbo Wang and Amir Leshem, "Decentralized Power Allocation
% and Beamforming Using Non-Convex Nash Game for Energy-Aware mmWave Networks"
%
% References:
% [1] R. A. Iltis, Seung-Jun Kim and D. A. Hoang, "Noncooperative iterative MMSE 
%     beamforming algorithms for ad hoc networks," in IEEE Transactions on 
%     Communications, vol. 54, no. 4, pp. 748-759, April 2006
%
% License: 
% This program is licensed under the GPLv3 license. If you in any way use this code for 
% research that results in publications, please cite our original article listed above.
 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU General Public License for more details.
%
function [RxBFUpdates, TxBFUpdates, PowerUpdates, flagInfeasible, TotalIteration] = TxBFIMMSEJointAllocation(InitRxBFs, ...
                                            InitTxBFs, InitPower, ...
                                            SystemSettings, IterationSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% The algorithm provided by [1] implicitly requires the Rx/Tx to be
% collocated, such that the Rx-Tx duality holds. It does not hold for our
% consided network (a general ad-hoc one), and does not produce convergent
% strategies.
%
% Input: 
%        InitRxBFs         : Intial Rx-beamformer vectors (not used)
%        InitTxBFs         : Intial Tx-beamformer vectors
%        InitPower         : Intial power values (single stream)
%        SystemSettings    : a struct of system settings, e.g., MIMO
%                            channels,SINR thresolds, maximum power values
%        IterationSettings : a struct of setting for terminating iterations
%
% Output: 
%         RxBFs            : Rx-beamforing vectors, Nlink*1 cell, each element Nr*1
%         TxBFs            : Tx-beamforing vectors, Nlink*1 cell, each element Nt*1
%         Powers           : Power values, Nlink*1
%         flagInfeasible   : flag to indicate whether there is a feasible
%                            solution 
%         TotalIteration   : Total number of iterations to convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Iteratively form the BF-vectors for both the Rxs and Txs
Epsilon = 1;
EpsilonThreshold = IterationSettings.EpsilonThreshold;
MaxIterNumber = IterationSettings.MaxIterNumber;

% Nt = SystemSettings.Nt;
Nr = SystemSettings.Nr;
Nlink = SystemSettings.Nlink;
sigma2 = SystemSettings.sigma2;
MIMOChannels = SystemSettings.MIMOChannels;
SINRThresholds = SystemSettings.SINRThreslods;

IterID = 1;

Powers = InitPower;
% RxBeamformers = InitRxBFs;
TxBeamformers = InitTxBFs;
while Epsilon > EpsilonThreshold && IterID <= MaxIterNumber
    RxBFUpdates = cell(Nlink, 1);
    TxBFUpdates = cell(Nlink, 1);
    PowerUpdates = zeros(Nlink, 1);
    for ii=1:Nlink
       %% Update receivers with MVDR beamformers      
        R_i = sigma2*eye(Nr);
        for jj=1:Nlink
            H_ji = MIMOChannels{jj,ii};
            TxBF_j = TxBeamformers{jj, 1}*Powers(jj, 1); %Tx Beamformer of j
            R_i = R_i + H_ji*(TxBF_j*TxBF_j')*H_ji';
        end
        H_ii = MIMOChannels{ii,ii};
        TxBF_i = TxBeamformers{ii, 1}*Powers(ii, 1);
        RxBF_i = R_i\(H_ii*TxBF_i);
        % Normalize
        RxBF_i = RxBF_i/norm(RxBF_i);
        RxBFUpdates{ii, 1} = RxBF_i;
        
        %%Update Tx beamformer
        TxBF_i = conj(RxBF_i);
        TxBFUpdates{ii, 1} = TxBF_i;
    end
    
    %% Update Tx power
    for ii=1:Nlink
        RxBF_i = RxBFUpdates{ii, 1};
        IpN_i = sigma2*(RxBF_i'*RxBF_i);% Interference plus noise
        for jj=1:Nlink
            if jj ~= ii
                TxBF_j = TxBFUpdates{jj, 1};
                H_ji = MIMOChannels{jj,ii};
                P_j = Powers(jj, 1);
                % set P_jj = PMax
                IpN_i = IpN_i + P_j*RxBF_i'*H_ji*(TxBF_j*TxBF_j')*H_ji'*RxBF_i;
            end
        end
        P_i = Powers(ii, 1);
        H_ii = MIMOChannels{ii,ii};
        TxBF_i = TxBeamformers{ii, 1}*P_i;
    	SINR_i =  real(P_i*RxBF_i'*H_ii*(TxBF_i*TxBF_i')*H_ii'*RxBF_i / IpN_i);
        SINRThreshold_i = SINRThresholds(ii, 1);
        PUpdate_i = P_i*SINRThreshold_i/SINR_i;
        PowerUpdates(ii, 1) = PUpdate_i;
    end
    
    %Update for the next iteration
    Powers = PowerUpdates;
    TxBeamformers = TxBFUpdates;
    IterID = IterID + 1;
end

% Note that the algorithm does not detect whether the problem is feasible
% or not
flagInfeasible = false;
TotalIteration = IterID;

end