# Intensify3D - Tutorial

## Examples for usage modes

In this section we will demonstrate  the use of Intensify3D under 2 examples:

1) Images with clear signal+backgroud and no non-relevant regions in the image

2) Images with no clear signal with non-relevant regions in the image (Tissue area)

3) Images with clear signal with non-relevant regions in the image (Tissue area)



# Example 1 - 3D Speres Images with clear signal+backgroud

<img src="SyntheticDataExample.png?raw=true." data-canonical-src="SyntheticDataExample.png?raw=true" />


 
This dataset was generated artificially, and a gradient of imaging intensity was generated in all xyz dimensions.
Different samples were generated by "different acquisition parameters" hence the baseline intensity is different. 
This Image does not have non-relevant regions so the option of **Tissue Detection is set to: No Detection**.
The Maximum Background Intensity was automatically selected and seems fitting since it detects a few pixels from the background. 
Apart from that we Assume that the background composition should ideally be similar at different depths so the in-stack normalization option is set to **Z-Normalization type is: Semi-Quantile**   
Between stacks, we chose a weaker assumption where we "only" want for the general intensity levels to be matched by Contrast-Stretching.

The following image shows the result of this run showing representative images from the image stacks. 

Find this dataset in the following link: [Artificial 3D spheres](https://drive.google.com/file/d/1tfSCFSalF4edfj0q2lJxSFqGs1gkj2XH/view?usp=sharing)


<img align="left" src="Montage2-01.jpg?raw=true." data-canonical-src="Montage.jpg?raw=true" />


# Example 2 - Auto-Fluorescent  Scanning of an iDISCO brain 

<img src="GUI_iDISCO_AutoFlu3.jpg?raw=true." data-canonical-src="GUI_iDISCO_AutoFlu3.jpg?raw=true" />

Down-sampled light sheet scan of an iDISCO cleared mouse hemisphere. Scanning iDISCO samples in GFP excitation/emission settings produces a distinct auto-fluorescent signal. while informative, this signal is easily photobleached and is prune for distortions. 
Here, since there is no clear definition of signal and background I set the MBI to max value. To accurately estimate the tissue area I set the background detection to threshold mode. in terms of z normalization the composition of imaged components varies along the z axis so semi-quantile is a too strong assumption so I selected contrast-stretch. Since this is a single sample I chose "No Normalization" for between sample normalization.

Find this dataset in the following link: [iDISCO Auto-Flu](https://drive.google.com/file/d/11k61eBUM8aNUg5Gf73U8hbwQbGFErtJ_/view?usp=sharing)

The result of this normalization is illustrated in the video below:


 ![Alt text](iDISCOHemi.gif?raw=true "Optional Title")


# Example 3 - Alexa647 Scan of an iDISCO brain
 
<img src="GUI_iDISCO.jpg?raw=true." data-canonical-src="GUI_iDISCO.jpg?raw=true" />
 
The following example is a cortical region of an iDISCO cleared mouse with staining against endogenously expressed tdTomato with alexa647. This sample has a clearly defined signal and a portion of the image that does represent any imaged tissue. For this reason I selected the tissue detection option, a threshold cutoff and a spatial filter that is large enough to accurately estimate the background with minimal interference of the signal. Here I selected Semi-quantile normalization for Z and contrast-stretch for between stack normalization. 


The result of this normalization is illustrated below:

<img src="iDISCO samples.jpg?raw=true." data-canonical-src="iDISCO samples.jpg?raw=true" />

 
 All rights reserved. No part of this software may be reproduced, 
 distributed, or transmitted in any form or by any means, including photocopying,
 recording, or other electronic or mechanical methods,
 without the prior written permission of the publisher,
 except in the case of a citation and certain other
 noncommercial uses permitted by copyright law.

