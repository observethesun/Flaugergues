# Digitized records

This folder contains digitized records separated in different files according to the method of observation:
* records_hv - horizontal and vertical wires
* records_ob - oblique wires
* records_rh - rhomboid scheme
* records_st - single transits
* records_w3 - 3-wires scheme 

## Data description
| Column name   |      Description      |
|----------|----------|
| obs_id |  Unique id of the record | 
| recorded_time |  Original time | 
| modified_time |  Modified time to correct possible miswritings |
| object| Object type: "S" for the solar disk, "t" for sunspot (tache in French) |
| event|  Name of the contacting wire. "?" means that the wire was not specified in the records. | 
| rising | True, if we assume that the record is before noon |
| pole|  Indication of the solar pole (following the original notation). "N" is for the north pole, "S" is for the south pole, "SP" is for the superieur bord (we assume it is the same as "N"), "IN" is for the inférieur bord (we assume it is the same as "S"). Black cells mean that the pole was not specified in the records and we could not unambiguously identify it from the records before or after this day.  | 

