# Pathway analysis with scBonita


## To perform pathway analysis, scBONITA uses the rules generated in Step 2. 

In addition, scBONITA requires:

* a **metadata file** specifiying the treatments/experimental variables for each cell and 
* a **contrast file** specifying the pairs of treatments to be compared.


## The pathway analysis script has the following arguments:

* **dataFile** 
    
    Specify the name of the file containing processed scRNA-seq data


* **conditions**
    
    Specify the name of the file containing cell-wise condition labels, ie, metadata for each cell in the training dataset. The columns are condition variables and the rows are cells. The first column must be cell names that correspond to the columns in the training data file. The column names must contain the variables specified in the contrast file (see contrast --help for more information).


* **contrast**
    
    A text file where each line contains the two conditions (corresponding to column labels in the conditions file) are to be compared during pathway analysis.


* **conditions_separator**
    
    Separator for the conditions file. Must be one of 'comma', 'space' or 'tab' (spell out words, escape characters will not work).",


* **contrast_separator**
    
    Separator for the contrast file. Must be one of 'comma', 'space' or 'tab' (spell out words, escape characters will not work).",

* **pathwayList**

    Paths to GRAPHML files that should be used for scBONITA analysis. Usually networks from non-KEGG sources, saved in GRAPHML format. The default value is to use sll the networks initially used for rule inference.
   

## Example usage with the provided example files in the `data` folder:

> `python3.6 pathwayAnalysis.py --dataFile "data/trainingData.csv" --conditions "data/conditions.txt" --contrast "data/contrast.txt" --conditions_separator "comma" --contrast_separator "comma" --pathwayList "hsa00010"`

## Output files from scBONITA Pathway Analysis

1. A comma-separated (CSV) file named as

> **pvalues + contrast[0] + _vs_ + contrast[1] + .csv**

For example, if the conditions to be compared (and specified in the contrasts file) are 'control' and 'treatment', the output file of scBONITA pathway analysis will be:

> **pvalues_control_vs_treatment.csv**

2. For each network, a file ending in **_importanceScores.csv**

This file contains a table with the following columns:

* **Node**: Gene name
* **Importance Score**: Importance score for the node in the network, calculated using the provided training dataset
* **ObsERS**: Observed size of the equivalent rule set or ERS. This is the number of possible equally valid Boolean rules for this node.
* **MaxERS**: Maximum possible size of the ERS for that node. This is (2^n) - 1, where n is the number of incoming edges for the node in the network.
* Two columns named as the two conditions to be compared/specified in the contrast file: the expression of the cells in those two conditions
* **Upregulated_in_condition** A boolean column stating whethere the node is up regulated in 'condition'

3. For each network, a GRAPHML file ending in "_IS" with node attributes similar to the ones in the importance scores CSV

## View and analyze the output of scBONITA

#### Load required packages


```python
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from singleCell import *
```

#### Example pvalues file


```python
pvalues = pd.read_csv("data/pvalues_Control_vs_Treatment.csv")
pvalues.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Pathway ID</th>
      <th>Pathway Name</th>
      <th>P value</th>
      <th>Contrast</th>
      <th>Upregulated_in_Control</th>
      <th>Adjusted P value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>data/hsa00010</td>
      <td>data/hsa00010</td>
      <td>0.158147</td>
      <td>Control_vs_Treatment</td>
      <td>False</td>
      <td>0.158147</td>
    </tr>
  </tbody>
</table>
</div>



#### Make bubbleplots using generated pvalues


```python
from singleCell import *

singleCell.makeBubblePlots(
    pd.read_csv("data/pvalues_Control_vs_Treatment.csv"),
    adjPValueThreshold=1,
    wrap=25,
    height=8,
    width=10,
    palette="colorblind",
    saveAsPDF=True,
    outputFile="example_PA_bubbleplot.pdf",
)
```


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <ipython-input-23-0f1bbf6b27ea> in <module>
          1 from singleCell import *
    ----> 2 singleCell.makeBubblePlots(pd.read_csv("data/pvalues_Control_vs_Treatment.csv"),
          3                     adjPValueThreshold=1,
          4                     wrap=25,
          5                     height=8,


    /gpfs/fs2/scratch/mpalshik/scBONITA-main/scBonita_package/src/scBONITA/singleCell.py in makeBubblePlots(pvalues, adjPValueThreshold, wrap, height, width, palette, saveAsPDF, outputFile)
       1519             the adjusted p-value threshold below which dysregulated pathways are shown on the bubbleplot
       1520         wrap: int
    -> 1521             wrap pathway names at this value on the y axis of the bubbleplot
       1522         height: float
       1523             height of image in inches


    ~/.local/lib/python3.8/site-packages/seaborn/_decorators.py in inner_f(*args, **kwargs)
         44             )
         45         kwargs.update({k: arg for k, arg in zip(sig.parameters, args)})
    ---> 46         return f(**kwargs)
         47     return inner_f
         48 


    ~/.local/lib/python3.8/site-packages/seaborn/relational.py in scatterplot(x, y, hue, style, size, data, palette, hue_order, hue_norm, sizes, size_order, size_norm, markers, style_order, x_bins, y_bins, units, estimator, ci, n_boot, alpha, x_jitter, y_jitter, legend, ax, **kwargs)
        806 
        807     variables = _ScatterPlotter.get_semantics(locals())
    --> 808     p = _ScatterPlotter(
        809         data=data, variables=variables,
        810         x_bins=x_bins, y_bins=y_bins,


    ~/.local/lib/python3.8/site-packages/seaborn/relational.py in __init__(self, data, variables, x_bins, y_bins, estimator, ci, n_boot, alpha, x_jitter, y_jitter, legend)
        585         )
        586 
    --> 587         super().__init__(data=data, variables=variables)
        588 
        589         self.alpha = alpha


    ~/.local/lib/python3.8/site-packages/seaborn/_core.py in __init__(self, data, variables)
        603     def __init__(self, data=None, variables={}):
        604 
    --> 605         self.assign_variables(data, variables)
        606 
        607         for var, cls in self._semantic_mappings.items():


    ~/.local/lib/python3.8/site-packages/seaborn/_core.py in assign_variables(self, data, variables)
        666         else:
        667             self.input_format = "long"
    --> 668             plot_data, variables = self._assign_variables_longform(
        669                 data, **variables,
        670             )


    ~/.local/lib/python3.8/site-packages/seaborn/_core.py in _assign_variables_longform(self, data, **kwargs)
        901 
        902                 err = f"Could not interpret value `{val}` for parameter `{key}`"
    --> 903                 raise ValueError(err)
        904 
        905             else:


    ValueError: Could not interpret value `log10pvals` for parameter `x`


#### Example importance scores file


```python
importanceScores = pd.read_csv(
    "data/hsa00010.graphml_processed.graphml_importanceScores.csv"
)
importanceScores.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Node</th>
      <th>Importance Score</th>
      <th>ObsERS</th>
      <th>MaxERS</th>
      <th>Control</th>
      <th>Treatment</th>
      <th>Upregulated_in_Control</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>DLAT</td>
      <td>0.003994</td>
      <td>127</td>
      <td>127</td>
      <td>0.012987</td>
      <td>0.043478</td>
      <td>False</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GALM</td>
      <td>0.881890</td>
      <td>16</td>
      <td>127</td>
      <td>0.038961</td>
      <td>0.043478</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>DLD</td>
      <td>0.003994</td>
      <td>127</td>
      <td>127</td>
      <td>0.116883</td>
      <td>0.260870</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>BPGM</td>
      <td>0.463312</td>
      <td>12</td>
      <td>127</td>
      <td>0.025974</td>
      <td>0.043478</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>MINPP1</td>
      <td>0.362319</td>
      <td>3</td>
      <td>7</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
</div>



#### Plot the equivalent rule set sizes for this network


```python
importanceScores.loc[importanceScores["MaxERS"] == 127].hist(column="ObsERS")
plt.xlabel("Observed ERS")
plt.ylabel("Frequency")
plt.title("ERS of nodes with in-degree >= 3)")
plt.show()
plt.clf()
```


    
![png](output_15_0.png)
    



    <Figure size 432x288 with 0 Axes>



```python
importanceScores.loc[importanceScores["MaxERS"] == 7].hist(column="ObsERS")
plt.xlabel("Observed ERS")
plt.ylabel("Frequency")
plt.title("ERS of nodes with in-degree = 2)")
plt.show()
plt.clf()
```


    
![png](output_16_0.png)
    



    <Figure size 432x288 with 0 Axes>


#### Visualize the output network in external software such as CytoScape or Gephi


```python
graph = nx.read_graphml("data/hsa00010_IS_Control_vs_Treatment.graphml")
```


```python
pd.DataFrame.from_dict(dict(graph.nodes(data=True)), orient="index")
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>DLAT</td>
      <td>DLD</td>
      <td>{'color': 'purple', 'subtype': 'compound', 'ty...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>DLAT</td>
      <td>POR</td>
      <td>{'signal': 'a'}</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GALM</td>
      <td>ADPGK</td>
      <td>{'color': 'purple', 'subtype': 'compound', 'ty...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GALM</td>
      <td>HK1</td>
      <td>{'signal': 'a'}</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GALM</td>
      <td>HK2</td>
      <td>{'signal': 'a'}</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>117</th>
      <td>ACSS1</td>
      <td>ALDH2</td>
      <td>{'signal': 'a'}</td>
    </tr>
    <tr>
      <th>118</th>
      <td>ACSS1</td>
      <td>ALDH1B1</td>
      <td>{'signal': 'a'}</td>
    </tr>
    <tr>
      <th>119</th>
      <td>ACSS1</td>
      <td>ALDH9A1</td>
      <td>{'signal': 'a'}</td>
    </tr>
    <tr>
      <th>120</th>
      <td>ACSS1</td>
      <td>ALDH3A2</td>
      <td>{'signal': 'a'}</td>
    </tr>
    <tr>
      <th>121</th>
      <td>ACSS1</td>
      <td>ALDH3B1</td>
      <td>{'signal': 'a'}</td>
    </tr>
  </tbody>
</table>
<p>122 rows × 3 columns</p>
</div>




```python
pd.DataFrame.from_dict(graph.edges(data=True))
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>name</th>
      <th>type</th>
      <th>importanceScore</th>
      <th>Observed ERS</th>
      <th>Max ERS</th>
      <th>relativeAbundance</th>
      <th>Control</th>
      <th>Treatment</th>
      <th>Upregulated_in_Control</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>DLAT</th>
      <td>DLAT</td>
      <td>gene</td>
      <td>0.003994</td>
      <td>127</td>
      <td>127</td>
      <td>-0.030491</td>
      <td>0.012987</td>
      <td>0.043478</td>
      <td>False</td>
    </tr>
    <tr>
      <th>GALM</th>
      <td>GALM</td>
      <td>gene</td>
      <td>0.881890</td>
      <td>16</td>
      <td>16</td>
      <td>-0.004517</td>
      <td>0.038961</td>
      <td>0.043478</td>
      <td>False</td>
    </tr>
    <tr>
      <th>DLD</th>
      <td>DLD</td>
      <td>gene</td>
      <td>0.003994</td>
      <td>127</td>
      <td>127</td>
      <td>-0.143986</td>
      <td>0.116883</td>
      <td>0.260870</td>
      <td>False</td>
    </tr>
    <tr>
      <th>BPGM</th>
      <td>BPGM</td>
      <td>gene</td>
      <td>0.463312</td>
      <td>12</td>
      <td>12</td>
      <td>-0.017504</td>
      <td>0.025974</td>
      <td>0.043478</td>
      <td>False</td>
    </tr>
    <tr>
      <th>MINPP1</th>
      <td>MINPP1</td>
      <td>gene</td>
      <td>0.362319</td>
      <td>3</td>
      <td>3</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>False</td>
    </tr>
    <tr>
      <th>AKR1A1</th>
      <td>AKR1A1</td>
      <td>gene</td>
      <td>0.003994</td>
      <td>127</td>
      <td>127</td>
      <td>0.124788</td>
      <td>0.298701</td>
      <td>0.173913</td>
      <td>True</td>
    </tr>
    <tr>
      <th>TPI1</th>
      <td>TPI1</td>
      <td>gene</td>
      <td>0.000000</td>
      <td>7</td>
      <td>7</td>
      <td>0.009599</td>
      <td>0.792208</td>
      <td>0.782609</td>
      <td>True</td>
    </tr>
    <tr>
      <th>GPI</th>
      <td>GPI</td>
      <td>gene</td>
      <td>0.335502</td>
      <td>16</td>
      <td>16</td>
      <td>0.107849</td>
      <td>0.194805</td>
      <td>0.086957</td>
      <td>True</td>
    </tr>
    <tr>
      <th>ADPGK</th>
      <td>ADPGK</td>
      <td>gene</td>
      <td>1.000000</td>
      <td>1</td>
      <td>1</td>
      <td>0.042349</td>
      <td>0.259740</td>
      <td>0.217391</td>
      <td>True</td>
    </tr>
    <tr>
      <th>ALDH9A1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.072464</td>
      <td>7</td>
      <td>7</td>
      <td>-0.039526</td>
      <td>0.090909</td>
      <td>0.130435</td>
      <td>False</td>
    </tr>
    <tr>
      <th>AHR</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.507246</td>
      <td>1</td>
      <td>1</td>
      <td>0.021457</td>
      <td>0.064935</td>
      <td>0.043478</td>
      <td>True</td>
    </tr>
    <tr>
      <th>PDHB</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.507246</td>
      <td>1</td>
      <td>1</td>
      <td>0.025409</td>
      <td>0.155844</td>
      <td>0.130435</td>
      <td>True</td>
    </tr>
    <tr>
      <th>GAPDH</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.239644</td>
      <td>56</td>
      <td>56</td>
      <td>-0.051948</td>
      <td>0.948052</td>
      <td>1.000000</td>
      <td>False</td>
    </tr>
    <tr>
      <th>PGM2</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.428571</td>
      <td>5</td>
      <td>5</td>
      <td>-0.013552</td>
      <td>0.116883</td>
      <td>0.130435</td>
      <td>False</td>
    </tr>
    <tr>
      <th>MDH1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.507246</td>
      <td>1</td>
      <td>1</td>
      <td>0.064370</td>
      <td>0.194805</td>
      <td>0.130435</td>
      <td>True</td>
    </tr>
    <tr>
      <th>MDH2</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.507246</td>
      <td>1</td>
      <td>1</td>
      <td>-0.010164</td>
      <td>0.337662</td>
      <td>0.347826</td>
      <td>False</td>
    </tr>
    <tr>
      <th>POR</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.072464</td>
      <td>7</td>
      <td>7</td>
      <td>-0.013552</td>
      <td>0.116883</td>
      <td>0.130435</td>
      <td>False</td>
    </tr>
    <tr>
      <th>ALDOA</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.000000</td>
      <td>1</td>
      <td>1</td>
      <td>0.018069</td>
      <td>0.844156</td>
      <td>0.826087</td>
      <td>True</td>
    </tr>
    <tr>
      <th>ALDOC</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.464976</td>
      <td>1</td>
      <td>1</td>
      <td>0.012987</td>
      <td>0.012987</td>
      <td>0.000000</td>
      <td>True</td>
    </tr>
    <tr>
      <th>ALDH2</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.072464</td>
      <td>7</td>
      <td>7</td>
      <td>0.003388</td>
      <td>0.220779</td>
      <td>0.217391</td>
      <td>True</td>
    </tr>
    <tr>
      <th>ALDH1B1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.072464</td>
      <td>7</td>
      <td>7</td>
      <td>0.012987</td>
      <td>0.012987</td>
      <td>0.000000</td>
      <td>True</td>
    </tr>
    <tr>
      <th>ALDH3A2</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.072464</td>
      <td>7</td>
      <td>7</td>
      <td>0.090909</td>
      <td>0.090909</td>
      <td>0.000000</td>
      <td>True</td>
    </tr>
    <tr>
      <th>ALDH3B1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.072464</td>
      <td>7</td>
      <td>7</td>
      <td>-0.048560</td>
      <td>0.168831</td>
      <td>0.217391</td>
      <td>False</td>
    </tr>
    <tr>
      <th>ADH5</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.003994</td>
      <td>127</td>
      <td>127</td>
      <td>0.034444</td>
      <td>0.077922</td>
      <td>0.043478</td>
      <td>True</td>
    </tr>
    <tr>
      <th>PDHA1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.507246</td>
      <td>1</td>
      <td>1</td>
      <td>-0.083004</td>
      <td>0.090909</td>
      <td>0.173913</td>
      <td>False</td>
    </tr>
    <tr>
      <th>LDHA</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.072464</td>
      <td>7</td>
      <td>7</td>
      <td>-0.145680</td>
      <td>0.506494</td>
      <td>0.652174</td>
      <td>False</td>
    </tr>
    <tr>
      <th>LDHB</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.072464</td>
      <td>7</td>
      <td>7</td>
      <td>0.120271</td>
      <td>0.337662</td>
      <td>0.217391</td>
      <td>True</td>
    </tr>
    <tr>
      <th>LDHC</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.072464</td>
      <td>7</td>
      <td>7</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>False</td>
    </tr>
    <tr>
      <th>PKM</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.003994</td>
      <td>127</td>
      <td>127</td>
      <td>-0.076228</td>
      <td>0.532468</td>
      <td>0.608696</td>
      <td>False</td>
    </tr>
    <tr>
      <th>ENO1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.255620</td>
      <td>64</td>
      <td>64</td>
      <td>0.079616</td>
      <td>0.688312</td>
      <td>0.608696</td>
      <td>True</td>
    </tr>
    <tr>
      <th>ENO2</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.259614</td>
      <td>63</td>
      <td>63</td>
      <td>0.025974</td>
      <td>0.025974</td>
      <td>0.000000</td>
      <td>True</td>
    </tr>
    <tr>
      <th>ENO3</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.259614</td>
      <td>63</td>
      <td>63</td>
      <td>0.012987</td>
      <td>0.012987</td>
      <td>0.000000</td>
      <td>True</td>
    </tr>
    <tr>
      <th>PGAM1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.255620</td>
      <td>64</td>
      <td>64</td>
      <td>-0.028797</td>
      <td>0.623377</td>
      <td>0.652174</td>
      <td>False</td>
    </tr>
    <tr>
      <th>PFKL</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.048309</td>
      <td>7</td>
      <td>7</td>
      <td>-0.027103</td>
      <td>0.233766</td>
      <td>0.260870</td>
      <td>False</td>
    </tr>
    <tr>
      <th>PFKM</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.066425</td>
      <td>7</td>
      <td>7</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>False</td>
    </tr>
    <tr>
      <th>PFKP</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.066425</td>
      <td>7</td>
      <td>7</td>
      <td>-0.104461</td>
      <td>0.025974</td>
      <td>0.130435</td>
      <td>False</td>
    </tr>
    <tr>
      <th>FBP1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.066425</td>
      <td>7</td>
      <td>7</td>
      <td>-0.118012</td>
      <td>0.142857</td>
      <td>0.260870</td>
      <td>False</td>
    </tr>
    <tr>
      <th>PGM1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.332126</td>
      <td>3</td>
      <td>3</td>
      <td>0.008470</td>
      <td>0.051948</td>
      <td>0.043478</td>
      <td>True</td>
    </tr>
    <tr>
      <th>HK1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.003661</td>
      <td>127</td>
      <td>127</td>
      <td>0.029362</td>
      <td>0.246753</td>
      <td>0.217391</td>
      <td>True</td>
    </tr>
    <tr>
      <th>HK2</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.003994</td>
      <td>127</td>
      <td>127</td>
      <td>-0.009034</td>
      <td>0.077922</td>
      <td>0.086957</td>
      <td>False</td>
    </tr>
    <tr>
      <th>HK3</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.944882</td>
      <td>8</td>
      <td>8</td>
      <td>-0.178995</td>
      <td>0.168831</td>
      <td>0.347826</td>
      <td>False</td>
    </tr>
    <tr>
      <th>HKDC1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.003994</td>
      <td>127</td>
      <td>127</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>False</td>
    </tr>
    <tr>
      <th>G6PC3</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.003994</td>
      <td>127</td>
      <td>127</td>
      <td>-0.026539</td>
      <td>0.103896</td>
      <td>0.130435</td>
      <td>False</td>
    </tr>
    <tr>
      <th>PGK1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.439347</td>
      <td>8</td>
      <td>8</td>
      <td>-0.033315</td>
      <td>0.662338</td>
      <td>0.695652</td>
      <td>False</td>
    </tr>
    <tr>
      <th>PCK2</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.507246</td>
      <td>1</td>
      <td>1</td>
      <td>-0.030491</td>
      <td>0.012987</td>
      <td>0.043478</td>
      <td>False</td>
    </tr>
    <tr>
      <th>ACSS2</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.507246</td>
      <td>1</td>
      <td>1</td>
      <td>0.064935</td>
      <td>0.064935</td>
      <td>0.000000</td>
      <td>True</td>
    </tr>
    <tr>
      <th>ACSS1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.507246</td>
      <td>1</td>
      <td>1</td>
      <td>-0.017504</td>
      <td>0.025974</td>
      <td>0.043478</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
</div>


