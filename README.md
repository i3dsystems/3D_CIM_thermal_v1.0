# 3D_CIM_thermal_v1.0
 
The 3D CIM thermal framework was developed by [Dr. Muhannad Bakir's group](https://bakirlab.gatech.edu/) (Georgia Institute of Technology). The model is made publicly available on a non-commercial basis. Copyright of the model is maintained by the developers, and the model is distributed under the terms of the [Creative Commons Attribution-NonCommercial 4.0 International Public License](http://creativecommons.org/licenses/by-nc/4.0/legalcode)

This research is supported by ASCENT, one of the SRC/DARPA JUMP centers.

If you use the tool or adapt the tool in your work or publication, you are required to cite the following reference: 

**A. Kaul; Y. Luo; X. Peng; M. Manley; Y.-C. Luo; S. Yu; M. S. Bakir, "3-D Heterogeneous Integration of RRAM-Based Compute-In-Memory: Impact of Integration Parameters on Inference Accuracy," in _IEEE Transactions on Electron Devices (TED)_, 2022, doi: [10.1109/TED.2022.3231570](https://ieeexplore.ieee.org/abstract/document/10003124).**

If you have any questions or comments on the framework, please contact  :man: [Dr. Muhannad Bakir](mailto:mbakir@ece.gatech.edu), and if you have technical questions or comments, please contact :man: [Ankit Kaul](mailto:ankit.kaul@gatech.edu) or :woman: [Madison Manley](mailto:madison.manley@gatech.edu).

## Files
1. Manual: Documents/3D_CIM_thermal_v1.0 Manual.pdf
2. Thermal model and RRAM device model python wrapper: 'thermal_model.py'
3. NeuroSim python wrapper: 'get_inference.py'

## Installation steps (Linux)
1. Get 3D_CIM_thermal from GitHub
```
git clone https://github.com/i3dsystems/3D_CIM_thermal_v1.0.git
```

2. Get DNN+NeuroSim tool from GitHub
```
git clone https://github.com/neurosim/DNN_NeuroSim_V1.3.git
```

3. Download MATLAB engine library

For usage, please refer to the manual at /Documents/3D_CIM_thermal_v1.0 Manual.pdf.

## References
[1] A. Kaul; Y. Luo; X. Peng; M. Manley; Y.-C. Luo; S. Yu; M. S. Bakir, "3-D Heterogeneous Integration of RRAM-Based Compute-In-Memory: Impact of Integration Parameters on Inference Accuracy," in _IEEE Transactions on Electron Devices (TED)_, 2022, doi: [10.1109/TED.2022.3231570](https://ieeexplore.ieee.org/abstract/document/10003124).

[2] Y. Zhang, Y. Zhang, and M. S. Bakir, “Thermal design and constraints for heterogeneous integrated chip stacks and isolation technology using air gap and thermal bridge,” IEEE Transactions on Components, Packaging and Manufacturing Technology, vol. 4, no. 12, pp. 1914–1924, 2014. [Online]. Available: 10.1109/TCPMT.2014.2364742

[3] P.-Y. Chen and S. Yu, “Compact modeling of rram devices and its applications in 1t1r and 1s1r array design,” IEEE Trans. Electron Devices, vol. 62, no. 12, pp. 4022–4028, 2015. [Online]. Available: 10.1109/TED.2015.2492421

[4] Y.-H. Chen, J. Emer, and V. Sze, “Eyeriss: A spatial architecture for energy-efficient dataflow for convolutional neural networks,” in Proceedings of the 43rd International Symposium on Computer Architecture, ser. ISCA ’16. IEEE Press, 2016, p. 367–379. [Online]. Available: https://doi.org/10.1109/ISCA.2016.40

[5] X. Peng, S. Huang, Y. Luo, X. Sun, and S. Yu, “Dnn+neurosim: An end-to-end benchmarking framework for compute-in-memory accelerators with versatile device technologies,” in Proc. IEEE Int. Elec. Devices Meeting, 2019, pp. 32.5.1–32.5.4. [Online]. Available: 10.1109/IEDM19573.2019.8993491
