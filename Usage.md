# Introduction #

The manual follows here


# Details #


# Usage #

## Inputs ##


<dem\_file>

 ESRI ascii Grid with the elevation height. Must have the same extent, resolution and coordinate system as the patch\_map


<patch\_map>

 ESRI ascii Grid with the patches identified by unique ID's

## Outputs ##


<output\_file>

 Writes a table to the output file where the landscape metrics for each patch ID are listed

## Comand ##
```
./patch3d <dem_file> <patch_map> <output_file>
```