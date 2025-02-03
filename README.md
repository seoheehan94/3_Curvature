# Computational Measurement of Curvature: Aligning Image Analysis with Human Behaviour and Neural Data
This study investigates how curvature influences visual perception and cognition by developing a novel method for computing curvature from contours in drawings. I compare this approach with curved Wavelet filters and human curvature ratings to analyze how curvature is represented in natural scenes. Using neural data, I examine how sensitivity to different degrees of curvature is distributed across the visual cortex. These findings provide new insights into computational curvature analysis and its role in visual processing.

[Ppt link](https://drive.google.com/file/d/1O5kGABUJNYFPZ00kdIj87_l478FKHiJw/view?usp=sharing)

---------
**Compute Curvilinearity**
 - contour curvature: [MLV toolbox](https://github.com/bwlabToronto/MLV_toolbox)
 - filter curvature: Wavelet filtering was based on an algorithm described by Kru¨ger and colleagues. (1996)

   (Kru¨ger NP, Peters G, Malsburg CVD (1996) Object recognition with an autonomously learned representation based on banana wavelets. Technical Report. Bochum, Germany: Institut für Neuroinformatik.) 

    Benchmark model:https://github.com/cechava/Rectilinearity_Toolbox

    (Nasr, S., Echavarria, C.E. and Tootell, R.B., 2014. Thinking outside the box: rectilinear shapes selectively activate scene-selective cortex. Journal of Neuroscience, 34(20), pp.6721-6735.)

---------
**curvature_contour**
- Image-computable model that predicts different degrees of curvature preference in the visual cortex using contour curvature

---------
**curvature_filter**
- Image-computable model that predicts different degrees of curvature preference in the visual cortex using filter curvature
