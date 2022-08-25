# Flaugergues

This repository contains digitized records and reconstructed sunspot positions from observations by amateur astronomer [H. Flaugergues](https://en.wikipedia.org/wiki/Honor%C3%A9_Flaugergues) in the period 1795–1830.

Destription of the method is presented in the paper [Sunspot positions 
from observations by Flaugergues in the Dalton minimum]()

## Contents
* [records]() - files with digitized records, separated by observation schemes
* coordinates.csv - reconstructed coordinates for each record
* sunspots.csv - sunspot positions averaged over multiple records for the day

## Data description

coordinates.csv

| Column name   |      Description      |
|----------|----------|
| obs_id |  Unique id of the record | 
| recorded_time |  Original start time of the record | 
| utc_time |  Estimated UTC time corresponding to the recorded_time |
| name | Sunspot name |
| lat |  Latitude | 
| long |  Carrington longitude | 
| cmd | Central meridian distance |
| type | Type of measurements: "hv" for horizontal and vertical wires, "ob" for oblique wires, "rh" for rhomboid scheme |
| x | x-coordinate in the primary coordinate system (see the paper for definition) |
| y | y-coordinate in the primary coordinate system (see the paper for definition) | 
| xs | x-coordinate in the coordinate system aligned with the line of solar motion |
| ys | y-coordinate in the coordinate system aligned with the line of solar motion | 
| l0 | Solar L0 angle for the time of observation | 
| b0 | Solar B0 angle for the time of observation |
| p | Solar P angle for the time of observation |
| incl | Angle between the line of solar motion and the parallel line of the telescope for the time and place of observation |
| incl_obs | Angle between the line of solar motion and the parallel line of the telescope derived from the record (if applicable) |
| err | Positioning error (if applicable) |
| r_sun | Solar "radius" measured in time units |
| modified | True, if original record was modified to correct possible miswritings | 
| complemented | True, if we complemented the original record with additional information that allows reconstruction of coordinates | 
| outlier | True, if we assume the record is invalid | 

sunspots.csv

| Column name   |      Description      |
|----------|----------|
| date |  Date of observation | 
| name |  Sunspot name | 
| lat_mean | Mean latitude |
| lat_std | Standard deviation for the latitude | 
| long_mean | Mean Carrington longitude |
| long_std | Standard deviation for the longitude | 
| cmd_mean | Mean central meridian distance |
| cmd_std | Standard deviation for the central meridian distance |
| n | Number of records used in averaging |
| confirmed | True, if there are multiple records on the same day or there are sunspots with similar coordinates before or after this day |

For description of the digitized data see the folder [records]().

## Citing this work
```
Illarionov, E., Arlt, R. Sunspot positions from observations by Flaugergues in the Dalton minimum. 2022
