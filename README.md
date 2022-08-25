# Flaugergues

This repository contains digitized records and reconstructed sunspot positions from observations by amateur astronomer [H. Flaugergues](https://en.wikipedia.org/wiki/Honor%C3%A9_Flaugergues) in the period 1795–1830.

Destription of the method is presented in the paper [Sunspot positions 
from observations by Flaugergues in the Dalton minimum]()

## Contents
* [records]() - files with digitized records, separated by observation schemes
* coordinates.csv - reconstructed coordinates for each record
* sunspots.csv - sunspot positions averaged over multiple records

## Data description

coordinates.csv

| Column name   |      Description      |  Data type |
|----------|----------|------|
| Obs_Id |  Unique id of the record | int |
| Recorded_Time |  Original start time of the record | str |
| UTC_Time |  Estimated UTC time corresponding to the Recorded Time | str |
| Column_1 |  Coordinate β | float |
| Column_2 |  Coordinate λ' | float |
