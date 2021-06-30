# Infrequently asked questions

## How is the percentile reported for genes with zero expression when more than one percentile of the samples in the compendium have no expression?

Answered by Ellen

Short answer: For gene A with 0 expression in a focus sample, the reported percentile for that sample+gene will be half (rounded up to the nearest full percent) of the highest percentile in the reference compendium for gene A which has a value of 0. For example, in polyA v9, gene AC006050.2 has 0 through 79 percentiles reported as 0.0, so any sample with AC006050.2 = 0 (eg, TH34_1447_S02) will report a percentile of 39.

Long answer:  0 routinely spans multiple percentiles. So we run into the dilemma of what to do when we have all these equal values (ie 0)  and we need to find our focus's rank within that big pile of identical things -- do we put it at the beginning, at the end, or in the middle? CARE chooses to put it at the middle.

This calculation is done in step 4 of CARE.
See https://github.com/UCSC-Treehouse/CARE/blob/f2adf596723e15e0dc4abd768a9c3c0915784253/care/4.0_outlier-analysis.ipynb 

Here's where it happens in the code precisely. This function is called once for every gene when a focus sample is analyzed
```
# Function for quantile expression. Used with pd.apply. takes a gene column
# from the binned cohort, the entire sample data frame, and the sample id
def get_sample_quantile(cohort_column, sample_df, sample_id):
    sample_value = sample_df.loc[cohort_column.name][sample_id] # a single number
    column_values = cohort_column.values # an array of the percentile values, 0 to 100
    
    leftmost = bisect.bisect_left(column_values, sample_value)
    rightmost = bisect.bisect_right(column_values, sample_value)
    
    return int((float(leftmost) + float(rightmost) - 1)/2)

cohort_column : The gene's column in the percentiles.hd5 file
sample_df, sample_id : the sample expression (all genes) and id
Going through line by line:
sample_value =...  get the focus sample's expression value for the desired gene (eg, 0.0)
column_values =... format the gene's column as an array so we can use bisect on it
```

Then, "bisect" is an algorithm - if we have a sorted array column_values  (eg the compendium percentiles) and a single value sample_value (ie the focus value), then bisect locates the insertion point for sample_value in column_values to maintain sorted order. 
We bisect twice (bisect_left and bisect_right)  to find the leftmost place we could possibly put sample_value and the rightmost. For almost every values this will be the same because the two adjacent percentiles to the value will be lower and higher than it, so there is only one reasonable place to put it.
However when 0 spans many percentiles, the leftmost will be 0 and the rightmost will be one past the "highest" 0.
For AC006050.2 (v9), here's a snippet of column_values:

```
0.00    0.0

0.01    0.0

0.02    0.0

...

0.77    0.000000

0.78    0.000000

0.79    0.000000

0.80    0.014498

0.81    0.028569
```


So bisect_right will give 80. 
Finally we do ((0  + 80) - 1 ) /2 and round up to get this "median" percentile that represents 0 expression, in this case 39.
(The minus one is because bisect is returning the index of insertion, not the index of the value already in the list so it's one higher than you would expect. For example if the sample values high and past the 100th %ile, bisect will return 101, not 100, so you need to go back one to get the percentile it's in.)
