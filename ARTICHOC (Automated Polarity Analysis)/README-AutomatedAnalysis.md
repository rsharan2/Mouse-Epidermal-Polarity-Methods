# ARTICHOC
## AUTOMATED RECOGNITION OF TWO IMMUNOSTAINS FOR COMPUTING HAIR FOLLICLE ORIENTATION AND COORDINATION

### 1. Preparing images for input
ARTICHOC takes as input average projection images in the .h5 format. To prepare
images for input, the user has two options:

  ##### A. In Fiji, 
  1) Duplicate the P-cadherin and Sox9 channels (if needed — if the image only
contains two channels, this may not be necessary). The order is important; ensure that
Ch0=Pcad and Ch1=Sox9. 2) Generate an Average z-projection using all z-planes with
Pcad or Sox9 signal. Image>Stack>Z Project>Average Intensity 3) Save file in .h5
format. File>Save as>”h5 (new or replace)”. The popup window will show “Dataset
names template;” make sure these names match the format /t{t}/channel{c}

 #####  B. Use the ARTICHOC_convert2h5.mlx script; 
 this can be run with the code hidden if desired. To set up for use, follow instructions at the top of the script, including
downloading the MIJI dependency to interface with ImageJ, adjusting the path names
to match your computer’s folders, and selecting your input file type. After set-up,
simply click run. This will pull up a menu, where you can select one or multiple files
from a given folder. You will need to click ‘save’ on each image, but otherwise the
script will do all the work to save the .h5 images in the same folder as your input
images.

### 2. Using ARTICHOC_analyzer.mlx to get angle and downgrowth data
The script can be run as-is and should give results, which are output as a data
table saved to the same folder as the input image and four figures in the MATLAB window:
the main output, showing arrows overlayed on the original image; the segmentation of both channels; the
“rings” generated during downgrowth analysis; and a dot plot of angles along the A-P axis to
quickly check overall trends. These can be each opened in a separate window by hovering over
each figure and clicking the arrow that appears, allowing the user to zoom in on a smaller
region. Some optimization is needed to get the best results given that size, resolution, and background may
vary between datasets. All parameters that the user will need to change are shown as sliders;
the easiest way to see all of these is to hide the code (under View>Hide Code). These fall
under three main sections — P-cadherin segmentation, Sox9 segmentation, and Analysis —
and two main classes — Intensity threshold and Size. In all cases, adjust the threshold parameters first. 
Size parameters can be adjusted if adjusting the thresholding does not improve results. The size parameters 
should maintain a similar ratio to the default values. 

  Note: ARTICHOC_analyzer_AdaptiveThresholding.mlx is identical except that it uses adaptive thresholding with a 
sensitivity value, rather than a fixed threshold value. It may perform better on some unevenly illuminated datasets.

### 3. Use ARTICHOC_hand_correction to correct and complete data.
ARTICHOC_analyzer may miss some data. The hand correction script gives the user the option to 
correct inaccurate data points, add angles for follicles that were overlooked in the initial output, and 
add new downgrown follicles. This is implemented by a MATLAB-ImageJ interface. Before the first use, the user
must download dependencies and add file locations to the script.

  ##### A. Download dependencies:
   i. MIJI library at https://imagej.net/Miji 
  
   ii. ARROWS function at https://www.mathworks.com/matlabcentral/fileexchange/37371-arrowsgeneralized-
and-vectorized-2-d-arrows-plot

  ##### B. Enter the path to Fiji script files and MIJI (may match the examples below)
	
   Windows:
		
		C:\Users\[username]\Programs\Fiji.app\scripts
		
		C:\Users\[username]\[your folder path]\mij\mij.jar
	
   Mac:
		
		/Users/[username]/Applications/Fiji.app/scripts
		
		/Users/[username]/[your folder path]/mij/mij.jar

		
   Choose the starting directory to select data files; this can be the same or different from MATLAB’s default directory.
	 
  ##### C. Ensure that your data and image files are in the same folder, and choose whether to:
		
   i. Manually select both the data and image files
		
   ii. Select the data file and auto-import the image file by matching the name
		
   iii. Adjust the data file prefix and suffix to match your dataset (if processed
by the ARTICHOC pipeline they should match the example given)
	
  ##### D. It may be necessary to allocate more memory to MATLAB.
Preferences>General>Java Heap Memory. 4096MB appears to be plenty; 1152MB is
not enough for working with larger images (2500 x 8000 pixels).

#### After this initial set-up, hand correction is implemented with three steps:
	
  ##### A. Run the code, and wait for a few seconds. 
  MATLAB will open an instantiation of
ImageJ and prepare the image and tools for analysis. All standard ImageJ menu options
are available. The channels tool will open automatically, so the user can view the
composite image if desired. A pop-up in the bottom right corner contains instructions
and a “continue” button (Figure 13).
	
  ##### B. Add new data points for polarized follicles.
		
   i. The “arrow” tool is selected to allow the user to add angle
measurements. Zoom in on a region and draw an arrow over the desired
follicle towards the P-cadherin (green) region. This can be done to over-write existing data 
(by drawing over an arrow or circle that is already present) or to add new data.

   ii. Press Command-M (Mac) or Ctrl-M (Windows) to measure each point. It may
be necessary to enter full screen so this doesn’t minimize the image window.
These arrows will be automatically added to the overlay.

   iii. When finished, press continue.
	
  ##### C. Add new data points for downgrown follicles.
		
   i. The “multipoint” tool will be selected; simply click
anywhere on a downgrown follicle to add it. Click on all downgrown follicles before continuing.

   ii. To check the results and/or add more points, press 'Continue Adding Points.' This will update 
the overlay and allow you to repeat the process in steps 3-4. To save a data table of what is 
currently on your screen, press Done. A new table called “[original file name] corrected.mat” 
will be saved to the same folder as the original files. Click 'Don't Save' if a pop-up asks to save the measurements or
changes to the image.
		
### 4. Figure creation
To create any of the figures shown here, use the ARTICHOC_figures.mlx script. This
script has several options laid out as check boxes and drop-down menus. To use,
select the desired plots from the six options: anterior and posterior rose plots, angle vs Anterior-Posterior
position dot plot (n=1 or n>1), Voronoi diagram overlay, angle map, and local order map. Select from several options 
(color, whether to include a color bar, number of neighbors to include in order calculation). Finally, choose whether or not to save
the generated plots to files, and set the destination folder. By default, these plots will be saved
to the same folder as the data file.
