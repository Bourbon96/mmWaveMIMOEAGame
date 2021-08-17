%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Wenbo Wang
%
% Original Article:
% [Wang2021] Wenbo Wang and Amir Leshem, "Decentralized Power Allocation
% and Beamforming Using Non-Convex Nash Game for Energy-Aware mmWave Networks"
%
% References:
% [1] O. E. Ayach, S. Rajagopal, S. Abu-Surra, Z. Pi and R. W. Heath, "Spatially Sparse 
%     Precoding in Millimeter Wave MIMO Systems," in IEEE Transactions on Wireless 
%     Communications, vol. 13, no. 3, pp. 1499-1513, March 2014, 
% [2] C. Miller, P. J. Smith, P. A. Dmochowski, H. Tataria and A. F. Molisch, "Favorable 
%     Propagation with User Cluster Sharing," 2020 IEEE 31st Annual International 
%     Symposium on Personal, Indoor and Mobile Radio Communications, 2020, pp. 1-7
%
% License: 
% This program is licensed under the GPLv3 license. If you in any way use this code for 
% research that results in publications, please cite our original article listed above.
 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU General Public License for more details.
%

function Ch_Matrix = MIMOChannelGenerationMP(Nt, Nr, Ncl, Nray, AS, Distance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
%         Nt            : Number of Tx antennas 
%         Nr            : Number of Rx antennas
%         Ncl           : Number of clusters
%         Nray          : Number of rays per cluster
%         AS            : Fixed angular spread at both the transmitter and receiver
%         Distance      : Distance between the receiver and the transmitter
% Output: 
%         Ch_Matrix     : Nr x Nt, MIMO channel matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate random mean AAOD, EAOD, AAOA, EAOA of each cluster

% Azimuth   : [-180,180] degree
% Elevation : [0,180] degree

% Transmitter(sectorized)
% Azimuth   : 60 degree wide w.r.t 0 degree
minAAOD = -30;
maxAAOD = 30; 

% Elevation : 20 degree wide w.r.t 90 degree
minEAOD = 80;
maxEAOD = 100;
Cluster_AAOD = rand(Ncl, 1) * (maxAAOD - minAAOD) + minAAOD; % [-30,30] degree
Cluster_EAOD = rand(Ncl, 1) * (maxEAOD - minEAOD) + minEAOD; % [80,100] degree

% Receiver(omni-directional)
Cluster_AAOA = (rand(Ncl, 1) - 0.5) * 360; % [-180,180] degree
Cluster_EAOA = rand(Ncl, 1) * 180; % [0,180] degree

%% Generate random AOD, EOD, AOA, EOA of rays per cluster (Laplacian distribution)

b = AS/sqrt(2); % Scaling parameter, degree

Randomness = rand(Nray*Ncl, 1)-0.5;

% Dimension of AOD, EOD, AOA, EOA : (Nray * Ncl) * 1 
RayAngle = @(cluster_angle, NumRay, RandSeed) (repelem(cluster_angle, NumRay, 1) ...
    - b*sign(RandSeed).* log(1-2.*abs(RandSeed)));

Ray_AAOD = RayAngle(Cluster_AAOD, Nray, Randomness);
Ray_EAOD = RayAngle(Cluster_EAOD, Nray, Randomness);
Ray_AAOA = RayAngle(Cluster_AAOA, Nray, Randomness);
Ray_EAOA = RayAngle(Cluster_EAOA, Nray, Randomness);

%% Obtain antenna element position vectors (normalized by half of the wavelength)

% Transmitter
Nt_H = sqrt(Nt); %horizon 
Nt_V = sqrt(Nt); %vertical

X_Tx = zeros(1,Nt);
[Y_Tx, Z_Tx] = meshgrid(0:Nt_H-1, 0:Nt_V-1);
TxPos = [X_Tx; Y_Tx(:).'; Z_Tx(:).']; % 3 x Nt

% Receiver
Nr_H = sqrt(Nr); %horizon 
Nr_V = sqrt(Nr); %vertical

X_Rx = zeros(1,Nr);
[Y_Rx,Z_Rx] = meshgrid(0:Nr_H-1, 0:Nr_V-1);
RxPos = [X_Rx; Y_Rx(:).'; Z_Rx(:).']; % 3 x Nr

%% Obtain array response vectors at the transmitter and receiver
SphericalUnitVecTx = get_spherical_unit_vector(Ray_EAOD, Ray_AAOD); % 3*Nray*Ncl
SphericalUnitVecRx = get_spherical_unit_vector(Ray_EAOA, Ray_AAOA); % 3*Nray*Ncl

%gg = sin(pi/Nt); % for UCA
gg = 1;

ArrayResponse_TX = (1/sqrt(Nt))*exp(1i*pi/gg*TxPos.'*SphericalUnitVecTx); % Nt*Nray*Ncl
ArrayResponse_RX = (1/sqrt(Nr))*exp(1i*pi/gg*RxPos.'*SphericalUnitVecRx); % Nr*Nray*Ncl

%% Generate complex path gain, using the exponential model
stdv = 10;
Alpha = exp(-1i*2*pi*(randn(1,Nray*Ncl)*stdv+Distance)/0.005);
% Alpha = sqrt(1/2) * (randn(1,Nray*Ncl) + 1i*randn(1,Nray*Ncl));

%% Generate MIMO channel matrix. See Eq. (3) of [Wang2021MIMO]
Ch_Matrix = sqrt((Nt*Nr)/(Nray*Ncl)) * ArrayResponse_RX * diag(Alpha) * ArrayResponse_TX';

end

function SUV = get_spherical_unit_vector(theta,phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%       theta, phi : M * 1
% Output:
%       SUV: 3 * M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For UCA, there is a shift for every antenna 2*pi*m/M for m=1,...,M, see [2]
M = size(phi, 1);
%phi_m = 2*180*(1:M)'/M;
phi_m = 0;
SUV = [(sind(theta).*cosd(phi - phi_m)).';...
                       (sind(theta).*sind(phi - phi_m)).';...
                       (cosd(theta)).']; 
end