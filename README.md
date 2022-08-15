# OBM_Primer_NPZD_model

A simple, 1-dimensional NPZD model in `Matlab`, based on Kuhn et al. (2015, doi: _[10.1016/j.pocean.2015.07.004](http://memg.ocean.dal.ca/memg/pubs/Kuhn_et_al_PiO_2015.pdf)_), with a Graphical User Interface (GUI).

The model is representative of a location at 50 degree N in the North Atlantic Ocean and the mixed layer evolution from this location is imposed. The model is run for 2 years, but only the 2nd year is shown in the auto-generated plots. The satellite-observed surface phytoplankton evolution at this location is shown for comparison in the surface property plot.

The GUI is called from the Matlab command line as follows:
`>> GUI_NPZD`

A window will pop up that allows the user to modify 4 parameters: the latitude of solar forcing, the initial nutrient concentration, the maximum phytoplankton growth rate, and the maximum zooplankton grazing rate. These parameters can be adjusted by using the sliders or typing new values directly into the corresponding editbox. If the value entered is outside the allowable range indicated by the slider, it will be reset to the corresponding allowable maximum or minimum.

The buttons below the four sliders allow the user to run the model or quit the GUI. The evolution of state variables at the surface layer will be written into a Matlab file along with essential meta-information.

## Feedback

Feel free use GitHub's issues feature (link above) to provide feedback.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
