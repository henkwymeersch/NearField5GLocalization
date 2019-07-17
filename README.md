# Near Field 5G Localization
This code computes position error bounds for near-field localization with a linear array. The code also comprises a spectral method for sub-array joint localization and synchronization.


## Summary
The GeneratePEB.m is used to compute the position error bound (PEB) for a scenario with a uniform linear array (ULA) under a variety of models: 
* the standard far-field, narrowband model
* the near-field, narrowband model
* the far-field wideband model
* the general model, applicable to far-field and near-field, wideband and narrowband. The PEB for the case where the clock bias is known or is unknown are visualized. 
The matlab code JointLocSync.m generates an OFDM signal and determines the observation matrix at the array (of size number of subcarriers x number of antenna elements). Then the matrix is processed to recover the transmitter's location and clock bias. The code evaluates the performance for different locations of the user. 

## Main parameters
Both codes have the same parameters
```
K = 1024;                       % number of subcarriers - 1
N = 128;                        % number of antennas -1
fc = 28;                        % carrier [GHz]
c = 0.3;                        % speed of light [m/ns]
lambda = c/fc;                  % carrier wavelength [m]                
Delta = lambda/2;               % antenna spacing in [m]
W = 0.1;                        % bandwidth [GHz]
Pt = 1;                           % transmit power [mW]
N0 = 290*1e3*1.381e-23*1e9;     % noise PSD in mW/GHz (290 Kelvin * Boltzmann constant in W/Hz)
B = 20;                         % user clock bias [m]
xUE = [1 8];                    % user location [m,m]
```

## Authors
The code was developed by 
* **[Henk Wymeersch](https://sites.google.com/site/hwymeers/)**

For more information, please see **[here](https://arxiv.org/abs/??.??)**
```
@article{Wym2019,
  title={Near-Field Joint Localization and Synchronization},
  author={Wymeersch, Henk},
  journal={arXiv preprint arXiv:??.??},
  year={2019}
}
```
## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
