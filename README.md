# Mocap-Data-Completion
(1) we utilize motion database to span two expressive spaces, namely, the frame space and trajectory space, to successfully process the new-come motion;

(2) improve the linear interpolation techniques applied to the gaps filling problem via a coordinate transformation, which depends on the pairwise training samples;

(3) utilize the high covariance both between frame coordinates and between marker coordinates to recover the missing data from the information provided by the available frames and markers respectively. 

## T method:
T0: gap interpolation on 3D 

T0_norepeat: gap interpolation with preprocessing on 3D data

T0_2D: gap interpolation on 2D data

## F method:
F0: gap interpolation on 3D data

F0_norepeat: gap interpolation with preprocessing on 3D data

F0_2D: gap interpolation on 2D data

## Extreme cases
T: joint interpolation

F1: frame interpolation on 3D data

F1_norepeat: frame interpolation with preprocessing on 3D data

F1_2D: frame interpolation on 2D data

## For rendering
testBVH: read bvh file and extract the xyz axis of all the joints over the whole video, and further render the motion

## Experiments
test4data: comparative exp against weighted PCA

test_marker_num: evaluate performance with increasing missing markers
