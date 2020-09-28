# STRAF to Tercen operator

##### Description

`straf2tercen` operator allows one to import STRAF input data into Tercen.

##### Usage

Input projection|.
---|---
`y-axis`        | numeric, input data, per cell 

Output relations|.
---|---
`median`        | numeric, median of the input data

##### Details
The operator takes all the values of a cell and returns the value which is the median.The computation is done per cell. There is one value returned for each of the input cell.

#### References


##### See Also

[mean_operator](https://github.com/tercen/mean_operator)
#### Examples
