# Flaugergues

This repository contains digitized records and reconstructed sunspot positions from observations by French astronomer [H. Flaugergues](https://en.wikipedia.org/wiki/Honor%C3%A9_Flaugergues) in the period 1788–1830.

Destription of the method is presented in the paper [Sunspot positions 
from observations by Flaugergues in the Dalton Minimum](https://doi.org/10.1093/mnras/stad1489)

## Contents
* [records](/records) - files with digitized records, separated by observation schemes
* [utils](/utils) - source functions for calculation of heliographic coordinates
* coordinates.csv - reconstructed coordinates for each record
* sunspots.csv - sunspot positions averaged over multiple records for the day

## Data description

coordinates.csv

| Column name   |      Description      |
|----------|----------|
| obs_id |  Unique id of the record (note the the record can contain several spots) | 
| recorded_time |  Original start time of the record |
| sidereal | True, if recorded time is sidereal | 
| utc_time |  Estimated UTC time corresponding to the recorded_time |
| name | Sunspot name |
| lat |  Latitude (degrees)| 
| long |  Carrington longitude (degrees)| 
| cmd | Central meridian distance (degrees)|
| outlier | True, if we assume the record is invalid | 
| type | Type of measurements: "hv" for horizontal and vertical wires, "ob" for oblique wires, "rh" for rhomboid scheme, "ec" for eclipse |
| x | x-coordinate in the coordinate system aligned with solar motion (see the paper for definition) |
| y | y-coordinate in the coordinate system aligned with solar motion (see the paper for definition) |  
| rising | True, if we assume that the record is before noon |
| l0 | Solar L0 angle (degrees) for the time of observation | 
| b0 | Solar B0 angle (degrees) for the time of observation |
| p | Solar P angle (degrees) for the time of observation |
| incl | Angle (degrees) between the line of solar motion and the parallel line of the telescope for the time and place of observation |
| incl_obs | Angle (degrees) between the line of solar motion and the parallel line of the telescope derived from the record (if applicable) |
| incl_diff | Calculated absolute difference (in degrees) between inclination of the line of solar motion at the beginning of the observation and at the end of the observation |
| err | Positioning error (degrees, if applicable) |
| r_sun | Solar "radius" measured in time units (seconds) |
| modified | True, if original record was modified to correct possible miswritings | 
| complemented | True, if we complemented the original record with additional information that allows reconstruction of coordinates | 

sunspots.csv

| Column name   |      Description      |
|----------|----------|
| date |  Date yyyy-mm-dd of observation | 
| name |  Sunspot name | 
| lat_mean | Mean latitude |
| lat_std | Standard deviation for the latitude | 
| long_mean | Mean Carrington longitude |
| long_std | Standard deviation for the longitude | 
| cmd_mean | Mean central meridian distance |
| cmd_std | Standard deviation for the central meridian distance |
| n | Number of records used in averaging |
| reliable | True, if there are sunspots with similar coordinates before or after this day |

For description of the digitized data see the folder [records](/records).

## Cite
```
Egor Illarionov, Rainer Arlt, Sunspot positions from observations by Flaugergues in the Dalton Minimum, Monthly Notices of the Royal Astronomical Society, 2023; https://doi.org/10.1093/mnras/stad1489
