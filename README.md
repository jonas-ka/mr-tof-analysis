# MR-ToF analysis

Created on 07 September 2019 for the ISOLTRAP experiment

@author: Jonas Karthein<br>
@contact: jonas.karthein@cern.ch<br>
@license: MIT license

### References
[1]: https://doi.org/10.1038/nature12226
[2]: https://doi.org/10.1103/PhysRevLett.120.232501
[3]: https://doi.org/10.1088/1674-1137/41/3/030002
[4]: https://www.etp-ms.com/technology/new_detectors
[5]: https://www.fastcomtec.com/products/software/

[1] F. Wienholtz, et al. Nature 498, 346 (2013).<br>
[2] M. Mougeot, et al. Phys. Rev. Lett. 120, 232501 (2018).<br>
[3] W.J. Huang, et al. Chinese Phys. C 41, 030002 (2017).<br>

### Application
The code was used to analyse data for the following publication:

[6] M. Mougeot, et al. in preparation (2019)<br>

### Introduction
The following code was written to reconstruct raw, unbinned multi-reflection time-of-flight (MR-ToF) data taken by a MagnetTOF electron multiplier particle detector by ETP Ion Detect [4] and using the MPANT software [5] for data aquisition. Besides a regular spectral time-of-flight (ToF) analysis with unbinned, maximum likelihood estimation and diverse cuts, the code includes a more advanced automated statistical cutting tool for ToF-drift detection based on rolling averaging. Finally, a mass value can be calculated from the fitted peaks by utilizing the Atomic Mass Evaluation (AME16) database [3].

### Required software and libraries
The following code was written in Python 3.7. The required libraries are listed below with a rough description for their task in the code. It doesn't claim to be a full description of the library.

- pandas (data storage and calculation)<br>
- numpy (calculation)<br>
- matplotlib (plotting)<br>
- mle (maximum likelihood estimation)<br>
- scipy (chi square fitting)<br>
- configparser (configuration file processing)<br>
- peakutils (peak finder)<br>
- jupyter (Python notebook environment)<br>

All packages can be fetched using pip:

`pip3 install pandas numpy matplotlib mle scipy configparser peakutils jupyter`
