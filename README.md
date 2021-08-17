# Decentralized Power Allocation and Beamforming Using Non-Convex Nash Game for Energy-Aware mmWave Networks

## Introduction
This project provides the simulation code in Matlab for the research paper [Wang2021]. Please feel free to contact me for further information/discussion regarding the technical details and implementation.

In this project, we focus on the problem of joint beamforming control and power allocation in the ad-hoc mmWave network. Over the shared spectrum, a number of multi-input-multi-output (MIMO) links attempt to minimize their supply power by simultaneously finding the locally optimal power allocation and beamformers in a self-interested manner. Our algorithm design is featured by separating the adaptation of power-levels and beamformer filters into two iterative sub-stages at each link, using a framework of generalized Nash game. For details about the theoretical proofs regarding the convergence of the two-stage algorithm, please refer to Sections III and IV of [Wang2021]. Several transmit beamforming schemes requiring different levels of information exchange are compared.

**Note**: The full Matlab implementation is available upon request via email before the completion of the review for our submitted paper. Currently, Tx-beamforming modules and part of the channel generation module is available.

## How to Use
The simulations can be started by running the .m files sarting with "main-".

## Content and Code Organization
- Simulation entrance
   - "mainMeshSimulation.m": to generate the mesh figures for Figure 1 of [Wang2021].
   - "mainConvergenceDemonstration": to generate the convergence curves in Figure 1 of [Wang2021].
   - "mainJointPowerBeamformerAllocation.m": to generate simulations results in Figures 2-4.
- Channel matrices generation
   - "MIMOChannelGeneration.m", "MIMOChannelGenerationMP.m", "GenerationPLOverDistance.m".
- Beamformer initialization
   - "InitializeBFs_QoS.m": to find the allowable SINR levels with best effort and initialize Tx/Rx beamformers.
   - "InitializeBFsMatchedFilter.m": to initialize Tx/Rx beamformers with myopic matched filters.
   - "InitializeBFsSimple.m": to initialize Tx/Rx beamformers with an MF-like primitive methods.
- Instantiation of joint allocation algorithms
   - "JointAllocationGame.m": the framework of the proposed algorithm in [Wang2021].
      - "OneStepJointAllocation.m": Single-step iteration of the proposed algorithm. Different Tx/Rx beamforming schemes are allowed to be used for strategy updating.
      - "OneStepBR.m": Single-step best response for power allocation sub-game.
   - "JointAllocationFullCoordination.m": the framework of the coordinated algorithm proposed by [1], [2].
      - "OneStepCoordinatedDualAllocation.m", "OneStepJointAllocation.m": Single-step iteration.
   - "JointAllocationIterativeZF.m": the framework of Coordinated algorithm proposed by [3], [4].
- Tx-beamforming Methods
   - All the .m files starting with the prefix of "TxBF-".
- Rx-beamforming Methods
  - All the .m files starting with the prefix of "RxBF-".
- Utility and miscellaneous functions
  - "PlotResults.m": for figure generation and saving in "mainJointPowerBeamformerAllocation.m".
  - "SpectrumEfficiency.m", "GetSupplyPower.m".

## Acknowledgement
The channel generation code is based on the model provided by [5].

## Original Article
[Wang2021]  Wenbo Wang and Amir Leshem, "Decentralized Power Allocation and Beamforming Using Non-Convex Nash Game for Energy-Aware mmWave Networks".

## References for the Implemented Algorithms
[1] H. Dahrouj and W. Yu, "Coordinated beamforming for the multicell multi-antenna wireless  system," IEEE Transactions on Wireless Communications, vol. 9, no. 5, pp. 1748–1759, 2010.

[2] D. H. N. Nguyen and T. Le-Ngoc, "Multiuser downlink beamforming in multicell wireless systems: A game theoretical approach, ”IEEE Transactions on Signal Processing, vol. 59, no. 7, pp. 3326–3338, 2011.

[3] H. Huh, A. M. Tulino and G. Caire, "Network MIMO With Linear Zero-Forcing Beamforming: Large System Analysis, Impact of Channel Estimation, and Reduced-Complexity Scheduling," in IEEE Transactions on Information Theory, vol. 58, no. 5, pp. 2911-2934, May 2012.

[4] Xiantao Sun, L. J. Cimini, L. J. Greenstein, D. S. Chany and J. Kruysy, "Coordinated zero-forcing beamforming in multipoint MIMO networks for backhaul applications," MILCOM 2009 - 2009 IEEE Military Communications Conference, 2009, pp. 1-7,

[5] O. E. Ayach, S. Rajagopal, S. Abu-Surra, Z. Pi and R. W. Heath, "Spatially Sparse Precoding in Millimeter Wave MIMO Systems," in IEEE Transactions on Wireless Communications, vol. 13, no. 3, pp. 1499-1513, March 2014,
