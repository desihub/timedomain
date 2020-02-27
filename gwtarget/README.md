# Gravitational Wave Scripts

## GW target notebook

The notebook in this folder will estimate the 50%/90% credible localization
interval of a GW alert and calculate all DR8 targets within the interval
given an r-band magnitude threshold.

## GW plots

Grab LIGO reconstruction from GraceDB (currently defaults to bayestar
reconstruction, not LALInference) and plot arbitrary credible intervals on a
HEALPix Mollweide map with pixels outside the DESI field of view masked out.

Example usage:

    python gw_skymap.py S191204r -d -t "BBH (>99%)"

To see options, run

    python gw_skymap.-h

and you will see usage instructions:

    usage: gw_skymap.py [-h] [-d] [-t EVENT_TYPE] event_id
    
    GW Event Plotter
  
    positional arguments:
      event_id              LVC event ID
    
    optional arguments:
      -h, --help            show this help message and exit
      -d, --display         Display plot
      -t EVENT_TYPE, --type EVENT_TYPE
                            Event type [BBH, BNS, ..]
