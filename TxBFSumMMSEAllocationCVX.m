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

function [TxBFLocal, PowerLocal] = TxBFSumMMSEAllocationCVX(LinkID, InitRxBFs, InitTxBFs, InitPower, SystemSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% (1) Return a Tx beamformer for link "LinkID", based on the global CSI and
% MSE criteria, to maximize the constrained sum of MSE in a cooperative manner
% Constraint: ||TxBFLocal||^2 <= PowerLocal.
% (2) Use CVX for the solution of TxBFLocal (Slow)
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

    % Solve the regularized optimization problem (closed-form solution)
    P_kk = InitPower(LinkID, 1);

    %% Solve the quadratic program for Tx-BF searching with CVX
    cvx_clear 
    cvx_begin
        cvx_quiet(true); % Suppress screen output by the solver
        cvx_precision low
        variable TxBF_k(Nt, 1) complex;
        expression ObjFun;
        
        ObjFun = 0;
        for ii=1:Nlink
            RxBF_i = InitRxBFs{ii, 1};
            H_ki = MIMOChannels{LinkID, ii};
            
            A_ki = H_ki'* RxBF_i;
            ObjFun = ObjFun + quad_form(sqrt(P_kk)*TxBF_k'*A_ki-1, 1);
        end
        
        minimize(ObjFun);
        subject to
            trace(TxBF_k'*TxBF_k) <= P_kk;        
    cvx_end

    % The gloval optimal MMSE receiver
    TxBFLocal = TxBF_k;
    TxBFLocal = TxBFLocal/sqrt(trace(TxBFLocal'*TxBFLocal)); %norm(TxBFLocal)
    
    % Update the power profile
    PowerLocal = trace(TxBF_k'*TxBF_k);    
end

