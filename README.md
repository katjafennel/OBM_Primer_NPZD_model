# OBM_Primer_NPZD_model

A simple 1-dimensional NPZD model in `Matlab`, based on Kuhn et al. (2015, doi: _[10.1016/j.pocean.2015.07.004](http://memg.ocean.dal.ca/memg/pubs/Kuhn_et_al_PiO_2015.pdf)_), with a Graphical User Interface (GUI).

The model is representative of a location at 50 degree N in the North Atlantic Ocean and the mixed layer evolution from this location is imposed. The model is run for 2 years, but only the 2nd year is shown in the auto-generated plots. The satellite-observed surface phytoplankton evolution at this location is shown for comparison in the surface property plot.

The GUi is called as follows from the Matlab command line:
`>> GUI_NPZD`

A window with the GUI will pop up that allows the user to modify 4 paramters, the latitude of solar forcing, the initial nutrient concentration, the maximum phytoplankton growth rate, and the maximum zooplankton grazing rate. These parameters can be adjusted by either using the slider to adjust values or by typing them directly into the corresponding editbox. If the number entered is outside the range of the slider, the number will be reset.

The buttons below the four sliders allow the user to run the model and quit the GUI. The evolution of state variables at the surface layer will be written into a Matlab file along with essential meta-information.
