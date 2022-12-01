# Community detection in empirical kinase networks

This repository contains the code to build empirical kinase networks and identify kinase signalling communities from quantitative phosphoproteomics data, as published in [placeholder].

- [Contributors](#contributors)
- [Prerequisites](#prerequisites)
- [Comparative wuantitative phosphoproteomics input](#comparative-quantitative-phosphoproteomics-input)
- [Calculate kinase-kinase interactions z-scores from phosphoproteomics data with KSEAR+](#calculate-kinase-kinase-interactions-z-scores-from-phosphoproteomics-data-with-ksear)
- [Network construction and community detection in Python](#network-construction-and-community-detection-in-python)
- [Example Analyses](#example-analyses)
- [Graph visualizations](#graph-visualizations)
- [References](#references)

## Contributors

All code in this repository was developed by [Celia Colomina Basanta](https://github.com/celiaccb) in collaboration with Marya Bazzi.

If you have any question regarding the computational analyses, you can raise an issue in this repository or contact [Celia Colomina Basanta](mailto:celiacolominab@gmail.com?subject=[KinaseCommunitiesDetection]) 

The [KSEAR+](https://github.com/CutillasLab/KSEA_plus) computational tool used to calculate the kinase-kinase interactions z-scores that serve as input for the pipeline was developed by [Pedro R. Cutillas group](https://github.com/CutillasLab).

## Prerequisites

### Access to github through GitHub Desktop or command line

We recommend that you run the pipeline by cloning this github repository, either using [GitHub Desktop](https://desktop.github.com/) or through your command line, in which case you might need to [install git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) depending on your operating system.

### KSEAR+

⚠ KSEAR+ is not compatible with mac os operating systems at the moment

Follow the instructions in the [KSEA_plus repository](https://github.com/CutillasLab/KSEA_plus) to install the KSEAR+ application.

### Python

Download Python 3 ([macOS](https://www.python.org/downloads/macos/), [Windows](https://www.python.org/downloads/windows/), [Linux/UNIX](https://www.python.org/downloads/source/), [others](https://www.python.org/download/other/))

If you want to run the jupyter notebook version of this pipeline (KinaseNetworks.ipynb), install [jupyter notebook.](https://jupyter.org/install)

The python packages listed below need to be installed prior to running the code in this repository. You can install all at once by running the following line on your command line (note that you need to have installed Python 3 first):
```
pip3 install pandas networkx scipy plotly pyvis statistics python-louvain matplotlib scikit-network
```

List of packages:
- pandas
- networkx
- NumPy (unless you have installed the Anaconda distribution of Python)
- SciPy
- plotly
- pyvis
- statistics
- louvain
- matplotlib
- scikit-network

<!-- ### R (optional)

Only necessary if you want to reproduce the network visualizations from the article (igraph_network_viz.R), in which case you will need to have [R](https://cran.r-project.org/) and [RStudio Desktop](https://www.rstudio.com/products/rstudio/#rstudio-desktop) installed. -->

## Comparative quantitative phosphoproteomics input

[KSEAR+](#calculate-kinase-kinase-interactions-z-scores-from-phosphoproteomics-data-with-ksear) calculates the difference in activity of pairwise kinase interactions based on the overall phosphorylation changes in treatment compared to control of the phosphosites that are identified of said pairwise kinase interactions.

See **144972_1_supp_335636_prw1dz_csv.csv** in the input folder of this repository for an example of the correct formatting for the comparative quantitative phosphoproteomics dataset to be considered a correct input by KSEAR+

## Calculate kinase-kinase interactions z-scores from phosphoproteomics data with KSEAR+

After downloading the KSEAR+ executable from the [KSEA_plus github repository](https://github.com/CutillasLab/KSEA_plus), open the application in your computer, upload your [comparative quantitative phosphoproteomics input](#comparative-quantitative-phosphoproteomics-input) and select the edges database.

The application returns six data frames: 

(1)z-scores 

The difference in the overall activity of kinase pairwise interactions in treatment compared to the control, as calculated based on the phosphorylation changes of its mapped substrates, using the [KSEA](#ksear-publications) method 

(2)distance

The difference in the overall activity of kinase pairwise interactions in treatment compared to the control, as calculated alternatively by simply getting the median of the phosphorylation changes of phosphosites in the kinase group minus the median of all sites

(3)pvalues 

The z-score p-value

(4)m

The number of phosphosites in the [comparative quantitative phosphoproteomics input](#comparative-quantitative-phosphoproteomics-input) that were identified as substrates of the pairwise kinase interaction, and based on which the z-score was calculated.

(5)q

The number of substrates (see m above) that are high in expression.

(6) sites

The list of matched substrates for the pairwise kinase interaction (see (4) m).

For more in-depth explanations, see the [KSEA publications](#ksear-publications).

For the following analyses you only need to save the (1)z-scores dataset as a csv (only csv will work with the code), see the **PROJECT_DATASET_2.csv** from [example 1](#example-analyses) in the input folder of this repository, or **example2_edges_zscores.csv** from [example2](#example-analyses) in the the input folder of this repository for examples of how your saved csv should look like.

## Network construction and community detection in Python

All code necessary to construct networks from edges z-scores and identify kinase signalling networks can be found in the KinaseNetworks script, which is available as a jupyter notebook (KinaseNetworks.ipynb) or a python script (KinaseNetworks.py).

The input for this script is the kinase-kinase interactions z-scores dataset generated from a quantitative phosphoproteomics dataset in the [previous step](##calculate-kinase-kinase-interactions-z-scores-from-phosphoproteomics-data). 
See, for example, PROJECT_DATASET_2.csv in the input folder, which is used for the first [example analysis](#example-analyses), detailed under the '2 - Analysis and Results' section of KinaseNetworks.

### Overview

The first section of this script ('1 - Functions') contains the four python functions we developed, which do the following:

1.  Given a kinase-kinase interactions z-scores dataset, separate these into each of the cell treatment, in each treatment extract the kinase-kinase interactions with a negative z-score (i.e, their activity is estimated to be decreased in that treatment compared to a control), and the individual kinases that are involved in these kinase-kinase interactions (function *dataset_info*).

2. For each treatment, construct a network from the downregulated kinase-kinase interactions extracted in step 1, and extract kinase communities from them with the louvain algorithm (function *converge_asso*).

3. From the kinase communities extracted for a given treatment in the previous step, the user can specify a kinase of interest, the function will return a list of all the kinases that belong to the same community as the kinase of interest (function *community_contents*).

4. Run step 2 and step 3 for a number of iterations specified by the user, the function returns the kinases (if any) that are assigned to the same community as the kinase of interest in some iterations, but are assigned to a different community in other iterations (function *rep_check*).

To see these functions applied to real-life datasets, see the [example analyses.](#example-analyses)

### Running the script

To run this script, clone this repository to GitHub Desktop or your command line, then open the repository in your command line.

- You can access the jupyter notebook in your browser by running ```jupyter notebook``` in your command line, then navigating to KinaseNetworks.ipynb
- If you are only interested in reproducing the published analyses, you can run all functions at once by placing your cursor in the cell above '2 - Analysis and Results' and clicking Cell>Run All Above, then skip to [Example Analyses](#example-analyses).

<!-- - If you prefer to run the python script, simply run ```python KinaseNetworks.py``` on your command line. -->

## Example Analyses

### 1 - Leukemic cell line treated with kinase inhibitors that target signaling pathways PI3K/AKT/mTOR and MEK/ERK 

Phosphoproteomics data of P31/FUJ cells treated separately with the kinase inhibitors GDC0941, AZD5363, trametinib and GDC0994, whose main targets are PIK3CA, AKT1/2, MAP2K1 ans MAPK1/3, respectively (this data is not publicily available at the moment).

The fold difference in the phosphorylation levels of an array of phosphosites between P31/FUJ cells treated separately with the kinase inhibitors GDC0941, AZD5363, trametinib and GDC0994 (whose main targets are PIK3CA, AKT1/2, MAP2K1 ans MAPK1/3, respectively) compared untreated P31/FUJ cells was calculated (this data is not publicily available at the moment). 

The difference in activity (i.e. z-score) of kinase-kinase interactions in treatment compared to control was calculated using KSEAR+ ( as [described above](#calculate-kinase-kinase-interactions-z-scores-from-phosphoproteomics-data-with-ksear) ) on the data described above.

The resulting z-scores file **PROJECT_DATASET_2.csv** can be found in the input folder in this repository.

We constructed negative networks (constructed based on the kinase pairwise interactions with a positive z-score (i.e. upregulated in treatment compared to control)) and positive networks (constructed based on the kinase pairwise interactions with a positive z-score (i.e. upregulated in treatment compared to control)) for each of the treatments in **PROJECT_DATASET_2.csv**, and applied community detection methods to each networks to identify signalling pathways, as can be found in the Examples section of the **KinaseNetworks.ipynb** jupyter notebook in this repository.

<!-- In KinaseNetworks.ipynb, you can run this analysis by placing your cursor on the last cell of section 2.5, then click Cell>Run All Above, which will return the following outputs: -->

[placeholder]

### 2 - Oestrogen receptor positive breast cancer cells treated with mTOR inhibitor rapamycin and oestrogen

[Cuesta et al. (2019)](#dataset-for-example-analysis-2) performed a phosphoproteome analysis of ER-positive MCF7 breast cancer cells treated with estrogen or oestrogen and the mTORC1 inhibitor rapamycin. They calculated the fold difference in the phosphorylation levels of ~9500 phosphosites in MCF7 cells treated with oestrogen and the mTORC1 inhibitor rapamycin compared to untreated MCF7 cells, which can be found in their Supplemental Data as **144972_1_supp_335636_prw1dz.xlsx**, and in the input folder of this repository as **144972_1_supp_335636_prw1dz_csv.csv**.

We calculated the z-scores of kinase pairwise interactions representative of the difference in activity in MCF7 cells treated with oestrogen and the mTORC1 inhibitor rapamycin compared to untreated MCF7 cells using [KSEAR+]((#calculate-kinase-kinase-interactions-z-scores-from-phosphoproteomics-data-with-ksear) on the **144972_1_supp_335636_prw1dz_csv.csv** dataset, resulting in the **example2_edges_zscores.csv** that can be found in the input folder of this repository.

We constructed a network based on the kinase pairwise interactions with a negative z-score, which in this case means their activity is calculated to be decreased in treatment compared to control, and identified tight-knit communities of kinases from it, as shown in the Examples section of the **KinaseNetworks.ipynb** jupyter notebook in this repository.

## Graph visualizations

<!-- For the [example analyses](#example-analyses) shown above, the communities of interest were visually represented as a graph in the [published manuscript](#manuscript).

You can reproduce these graphs by running the igraph_network_viz.R script (from RStudio):
- If you have cloned the repository into GitHub Desktop, click on 'Repository' on your navigation bar, then click 'Open in RStudio'.
- Otherwise, copy the repository URL, then open RStudio and navigate to File>New Project>Version Control>Git. Paste the repository URL and click Create Project. -->

## References

### Manuscript

[placeholder]


### Dataset for Example Analysis 2

Cuesta R, Gritsenko MA, Petyuk VA, Shukla AK, Tsai CF, Liu T, McDermott JE, Holz MK. Phosphoproteome Analysis Reveals Estrogen-ER Pathway as a Modulator of mTOR Activity Via DEPTOR. Mol Cell Proteomics. 2019 Aug;18(8):1607-1618. doi: 10.1074/mcp.RA119.001506. Epub 2019 Jun 12. PMID: 31189691; PMCID: PMC6683011.

### KSEAR+ publications

Casado P, Rodriguez-Prados JC, Cosulich SC, Guichard S, Vanhaesebroeck B, Joel S, Cutillas PR. Kinase-substrate enrichment analysis provides insights into the heterogeneity of signaling pathway activation in leukemia cells. Sci Signal. 2013 Mar 26;6(268):rs6. doi: 10.1126/scisignal.2003573. PMID: 23532336.

Casado P., Hijazi M., Gerdes H., Cutillas P.R. (2022) Implementation of Clinical Phosphoproteomics and Proteomics for Personalized Medicine. In: Corrales F.J., Paradela A., Marcilla M. (eds) Clinical Proteomics. Methods in Molecular Biology, vol 2420. Humana, New York, NY. [https://doi.org/10.1007/978-1-0716-1936-0_8](https://doi.org/10.1007/978-1-0716-1936-0_8)
