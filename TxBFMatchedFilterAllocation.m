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

function [TxBFLocal, PowerLocal] = TxBFMatchedFilterAllocation(LinkID, InitRxBFs, InitTxBFs, InitPower, SystemSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% Return a matched filter beamformer for link "LinkID", based on the input InitTxBfs
% See Section IV.C.1.
%
% Input: 
%        LinkID            : The ID of the link of interest
%        InitRxBFs         : Intial Rx-beamformer vectors
%        InitTxBFs         : Intial Tx-beamformer vectors (not used)
%        InitPower         : Intial power values (not used)
%        SystemSettings    : a struct of system settings, e.g., MIMO
%                            channels,SINR thresolds, maximum power values (not used)
%
% Output: 
%         TxBFLocal        : Tx-beamforing vectors, Nt*1
%         PowerLocal       : Local power strategy (unchanged)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RxBFLocal = InitRxBFs{LinkID, 1};

H_ii = SystemSettings.MIMOChannels{LinkID,LinkID};

[TxBFLocal,~] = eigs(H_ii'*(RxBFLocal*RxBFLocal')*H_ii, 1);
PowerLocal = InitPower(LinkID, 1);
end