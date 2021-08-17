%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Wenbo Wang
%
% Original Article:
% [Wang2021] Wenbo Wang and Amir Leshem, "Decentralized Power Allocation
% and Beamforming Using Non-Convex Nash Game for Energy-Aware mmWave Networks"
%
% References:
% [1] H. Huh, A. M. Tulino and G. Caire, "Network MIMO With Linear Zero-Forcing Beamforming: 
%     Large System Analysis, Impact of Channel Estimation, and Reduced-Complexity Scheduling," 
%     in IEEE Transactions on Information Theory, vol. 58, no. 5, pp.
%     2911-2934, May 2012.
% [2] Xiantao Sun, L. J. Cimini, L. J. Greenstein, D. S. Chany and J. Kruysy, "Coordinated 
%     zero-forcing beamforming in multipoint MIMO networks for backhaul applications," 
%     MILCOM 2009 - 2009 IEEE Military Communications Conference, 2009, pp. 1-7, 
%
% License: 
% This program is licensed under the GPLv3 license. If you in any way use this code for 
% research that results in publications, please cite our original article listed above.
 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU General Public License for more details.
%

function [TxBF, Power] = ZeroForcingTXBF(LinkID, InitRxBFs, InitTxBFs, InitPower, SystemSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%
% Input: 
%        InitRxBFs         : Intial Rx-beamformer vectors 
%        InitTxBFs         : Intial Tx-beamformer vectors (not used)
%        InitPower         : Intial power values (single stream)
%        SystemSettings    : a struct of system settings, e.g., MIMO
%                            channels,SINR thresolds, maximum power values
%
% Output: 
%         TxBF             : Tx-beamforing vectors, Nt*1
%         Power            : Local power strategy (unchanged)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MIMOChannel = SystemSettings.MIMOChannels;
Nlink = SystemSettings.Nlink;

% Nr = size(MIMOChannel{1,1}, 1);
Nt = SystemSettings.Nt;

if Nlink >= Nt
    error('Zero-forcing require link number needs to be smaller than antenna number.');
end

%% For Tx-beamformers
HIntf_i = zeros(Nt, Nlink-1); % Equivalent interference channels
Hii = MIMOChannel{LinkID,LinkID};
% Composite Interference Channel (CIC) for each link, Nt*(Nlink-1) matrix
col_id = 1;
for jj=1:Nlink    
    if (jj ~= LinkID)
        RxBFj = InitRxBFs{jj, 1};
        H_ii_jj = MIMOChannel{LinkID,jj}; %zero-forcing the interference
        % equivalent channels to form CIC
        HIntf_i(:, col_id) = (RxBFj'*H_ii_jj).'; % h_{j,i} = [g^H_j*H_{j,i}]^T, Eq.(6) of Ref.[2]
        col_id = col_id + 1;
    end
end
    
HIntf_i = HIntf_i.';
% SVD for CIC, Nt - max(Nt, Nr), to get row- null space basis vectors of HIntf_i
[~, M, Vi] = svd(HIntf_i);
    
diag_m = diag(M);
rank_M = nnz(diag_m);
Vi0 = Vi(:, rank_M+1:end); % column null-space spanned by the basis vectors
    
RxBF_i = InitRxBFs{LinkID, 1};
HPrime = Vi0' * Hii' * (RxBF_i * RxBF_i') * Hii * Vi0;    
[b_opt_i, ~] = eigs(HPrime, 1); % Largest eigenvector
    
TxBF = Vi0 * b_opt_i; %  see Eq. (9) of Ref. [2]

sigma2 = SystemSettings.sigma2;
Interference = sigma2*(RxBF_i'*RxBF_i); % Interference plus noise Eq.(4) of [1]
for jj=1:Nlink
    if jj ~= LinkID
        H_ji = MIMOChannel{jj,LinkID};
        Power_jj = InitPower(jj, 1);
        TxBF_jj = InitTxBFs{jj, 1};
        Interference = Interference + Power_jj*RxBF_i'*H_ji*(TxBF_jj*TxBF_jj')*H_ji'*RxBF_i;
    end
end

SINRThreshold = SystemSettings.SINRThresholds(LinkID);
Power = SINRThreshold*real(Interference)/real(RxBF_i'*Hii*(TxBF*TxBF')*Hii'*RxBF_i);

end