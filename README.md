# Automated_GUI_based_heart_rate_analyzer

Progress Report : UI + Simple EVM

I have managed to circumvent the problem of performing the whole EVM on video for analysis, and now relying only on Band Pass Filter with FFT peak analysis to extract the data. Also, have managed to work only on the selected region(both on time axis and spatially, cropping the video) of the video rather than the whole video, as a result getting more reliable results each iteration and saving some computation time.

How to test robustness or accuracy of results ?
Results are reproducible multiple times, no random arte-effects arising
Ground truth can be a good reference point to test against - that can only be found when applied to a known database

Here are the new features provided for the UI development module of the signal from video for Heart Rate Analysis :

Use Video Analyzer application in MATLAB for analyzing/previewing the video and deciding the start and stop frames to be analyzed
Added simple User Interface feature to crop video (time cropping), and select region of interest of video, then process only that region of video (spatial cropping)
Comment : Even 1 second can be a long time for small resolution of activity , when person falls,etc so keeping the option of choosing frame numbers rather than seconds - video analyzer toolbox gives exact frame number to user 
Added simple UI to choose the color plane from a dialog box to choose R/G/B and see effects [which is better, more robust ?]
Added simple User Interface feature to add or set features for the band pass filter to be used like low, high frequencies, frame overlap (sliding window length), etc. 
Default values also set for user 
Condensed all files into 1 file and rest of the video cropping, etc things as functions 
Now 3 function files exist within 1 main file, hence very easy for the user to navigate to the matlab folder and run the video, and then just enter values on the prompt
Results shown very lucidly and interactively
Added drawnow feature to update figures as analysis of video frame by frame is done so as to provide a more immersive & interactive experience for the user
Performed 25 Iterations with various combinations :
Small versus Large ROI size
Short versus Long window sampling size [0.5 to 2.5 seconds]
R vs G vs B color channels
Small versus large range of BP Frequencies
Light conditions : well lit versus dim light
Human heart rate condition : after exercise versus normal 
Choose an IIR Butterworth filter versus FIR filter : lesser computations as order of complexity of IIR is less , and Flat pass-stop bands of the IIR filter means lesser ripples of edge frequencies 

To be analyzed :
EVM video of the same ? Does it work similarly or differently ? How much of accuracy difference between these two methods ?


Advantages :

DFT versus FFT : dft has o(n2) but fft is o(nlog(n))
Correlation as a measure of similarity: Goal was to find the tone most similar to the signal in a frequency range around the FFT peak, by measuring the similarity (correlation) between the signal and a series of reference complex tones/other frequncies
 So this method is equivalent to computing the FFT on a zero-padded version of the signal, which is a commonly used technique to smooth the FFT data, and picking only a frequency range. If the original 6-second signal in the window is zero-padded up to 60 seconds, the tone frequencies used in this method are found among the orthonormal frequencies of the 60-second FFT. Then, the method is equivalent to applying the DFT definition on the zero-padded window for certain frequencies only

To do further :

Plot graph of frequency versus amplitude : can be interesting to see what are the various amplitudes of various frequency components in the signal & hence draw a conclusion for the mean or most dominant frequencies present
Plot time series versus amplitude plot for whole of the cropped video 
Plot the time series versus amplitude plot for the whole of cropped video + EVM of video for comparison of the two
Give a video player GUI to the user for selecting the time (in seconds and not frame number)
Give a gui to browse the video file from the folder directory
Change small things in the dialog boxes
See the app interface given by prof martin and take some ideas from it

