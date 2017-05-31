[![DOI](https://zenodo.org/badge/92933745.svg)](https://zenodo.org/badge/latestdoi/92933745)

A collection of pRF scripts. The code was used for the pRF-mapping in:

* Ekman, Kok & de Lange (2017). Time-compressed preplay of anticipated events in human primary visual cortex. *Nature Communications*. 8, 15276 doi: 10.1038/ncomms15276 https://www.nature.com/articles/ncomms15276

* Ekman, Roelfsema & de Lange (in prep). Automatic spread of the attentional field from a spatial cue to an underlying object in primary visual cortex.



## Stimulus code

Matlab scripts originally developed by Dumoulin & Wandell 2008 *NeuroImage*.
Please note that these scripts were written in 2013 and a newer version is available [here](https://github.com/vistalab/vistadisp/tree/master/Applications2/Retinotopy/standard). If you're starting from scratch, that ~~might be~~ *is* a better place to look.

* the stimulus scripts require Matlab (tested with matlab2011b) and [Psychtoolbox3](http://psychtoolbox.org/).
* start the menu with ``run('ret_initialise.m')``.
* choose "Translating Bars 8 pass (Dumoulin)" and load your scanner calibration (e.g. "Siemens Skyra"). Note that you can create your own calibration; have a look at ``/MRstim/Displays/Skyra/displayParams.m`` for an example.
* stimulus presentation will start with the first MRI trigger.

We run 4 blocks (stimulation in each block is identical), but more blocks are always better.

## Analysis code

The analysis is based on [MrVista](https://github.com/vistalab) and requires [SPM8](http://www.fil.ion.ucl.ac.uk/spm/software/spm8/); other SPM versions might work as well.

* Preprocessing SPM/FSL
* ... #TODO

## Citing

If you use these scripts please consider citing the following articles:

```
Dumoulin & Wandell (2008). Population receptive field estimates in human visual
cortex. Neuroimage 39, 647–660. http://www.sciencedirect.com/science/article/pii/S1053811907008269

Ekman, Kok & de Lange (2017). Time-compressed preplay of anticipated events in
human primary visual cortex. *Nature Communications*. 8, 15276
doi: 10.1038/ncomms15276 https://www.nature.com/articles/ncomms15276
```

## Related

analyzepRF, popeye, samprf
