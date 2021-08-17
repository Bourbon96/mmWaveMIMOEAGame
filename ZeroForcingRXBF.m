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

function RxBF = ZeroForcingRXBF(LinkID, InitTxBFs, InitPower, SystemSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   This part of code is not fully tested.
%
% Input: 
%         MIMOChannel           : Channel coeff. matrices, N*N cell, each element Nr*Nt
%         TxBeamformers         : Rx-beamforing vectors, N*1 cell, each element Nr*1 
%                                 (for iterative LZF purpose)
% Output: 
%         RxBeamformers         : Rx-beamforing vectors, N*1 cell, each element Nr*1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MIMOChannel = SystemSettings.MIMOChannels;

    Nlink = SystemSettings.Nlink;
    Nr = SystemSettings.Nr;

    if Nlink >= Nr
        error('Zero-forcing require link number needs to be smaller than antenna number.');
    end

    %% For Rx-beamformers
    HIntf_i = zeros(Nr, Nlink-1);
    Hii = MIMOChannel{LinkID,LinkID};
    % Composite Interference Channel (CIC) for each link, Nr*(Nlink-1)
    col_id = 1;
    for jj=1:Nlink
        if (jj ~= LinkID)
            TxBF_jj = InitTxBFs{jj,1};
            H_jj_ii = MIMOChannel{jj,LinkID};
            HIntf_i(:, col_id) = H_jj_ii * TxBF_jj; % h_{j,i} = H_{j,i}*q_j
            col_id = col_id + 1;
        end
    end
    
    % SVD for CIC, (Nt - max(Nt, Nr), to get column null-space matrix of HIntf_i
    [Ui, M, ~] = svd(HIntf_i);

    % Note that the discussion about the ZF Rx-BF in Ref.s[2] is not correct
    diag_m = diag(M);
    rank_M = nnz(diag_m);
    Ui0 = Ui(:, rank_M+1:end); % matrix spanned by the row null-space vectors
    
    TxBF_i = InitTxBFs{LinkID, 1};
    HPrime = Ui0' * Hii * (TxBF_i * TxBF_i') * Hii' * Ui0; %same as TXBF
    [a_opt_i, ~] = eigs(HPrime, 1);
    
    gi =  Ui0 * a_opt_i;
    RxBF = gi;
end