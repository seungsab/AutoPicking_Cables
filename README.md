# **[Release original source code: Tutorial using Python]**

## **[Python code (Jupyter notebook)](./00_Auto%20Peak%20Picking%20based%20on%20Domain%20Knowledge%20for%20Cables.ipynb.ipynb) (Uploaded on 2024.01.04)**

----

# Automated Peak-Picking for Cable Structures

Fully-automated peak-picking method using modified AMPD (Automatic multiscale-based peak detection) [1] and MAD (Median Absolute Deviation) after Baseline-correction [2].

The source codes in this repository are uploaded for the replication of results section.
- S.S. Jin, S.H. Jeong*, S.H. Sim, D.W. Seo and Y.S. Park (2021), Fully automated peak-picking method for an autonomous stay-cable monitoring system in cable-stayed bridges, Automation in Construction, 126, 103628, https://doi.org/10.1016/j.autcon.2021.103628.

All source codes were developed by Dr. Jin.

The results of the Faster R-CNN were obtained using Tensorflow in PYthon by Dr. Jeong.

The source codes are only working with the uploaded datasets, since the algorithm is patent pending.

If you are interesting in our method, please contact to "seungsab@gmail.com".

**[Source code]**

Three automatic peak-pikcing methods can be performed with the results of the Faster R-CNN by running the following m-file
- Main_Run_automatic_peak_picking.m: Run four automatic peak-pikcing methods for each class

**[Data]**

Three levels of difficulty for automatic peak-picking are defined.
- f_PSD_data_class_easy_data.mat: 20 PSDs for easy class
- f_PSD_data_class_intermediate_data.mat: 20 PSDs for intermediate class
- f_PSD_data_class_hard_data.mat: 20 PSDs for hard class

The raw signals and their PSD were plotted by jpg-files in the subfolder of "Raw_Data"

The peaks selected by Faster R-CNN (tensorflow in python)
- Faster_RCNN_Box_final.mat

The peaks selected by five experts (manually peak-picking)
- Peak_expert_Easy.mat: 
- Peak_expert_Intermediate.mat: 
- Peak_expert_Hard.mat: 

The execution time (second) for peak detection by the pre-trained Faster R-CNN.
- time_consumption_easy.txt: execution time for easy class
- time_consumption_inter.txt: execution time for intermediate class
- time_consumption_hard.txt: execution time for hard class

**[Result]**

All results are in the folders as "Result" with subfolders as follows:
- Easy: 20 results in jpg-file
- Intermediate: 20 results in jpg-file
- Hard: 20 results in jpg-file

If you have any questions, please feel free to contact me: seungsab@gmail.com (Seung-Seop Jin).

Note) These source codes are developed in MATLAB 2020a. They may not be working below MATLAB R2020a.

[1] Automatic Multiscale-based Peak Detection (AMPD)
- F. Scholkmann, J. Boss, M. Wolf, An Efficient Algorithm for Automatic Peak Detection in Noisy Periodic and Quasi-Periodic Signals, Algorithms 5 (4) (2012) 588-603.

[2] Baseline estimation for Raman Spectra
- H.G. Schulze, R.B. Foist, K. Okuda, A. Ivanov, R.F.B. Turner, A Small-Window Moving Average-Based Fully Automated Baseline Estimation Method for Raman Spectra, Applied Spectroscopy 66 (7) (2012) 757-764.
