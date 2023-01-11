# Peptimetric - code by Hartman E. and Mahdavi S.
This git repository contains the code and data for a research project conducted by Erik Hartman, Simon Mahdavi, Sven Kjellström and Artur Schmidtchen. The project is hosted by heroku and live on https://peptimetric.herokuapp.com/.

### To run peptimetric locally
If servers are down or you would like to deposit large files to peptimetric - it can be cloned and used locally. This requires python >= 3.7 and poetry (https://python-poetry.org/docs/). 

#### Steps:
1. Clone the repo
2. Open the repo and run 	`poetry install` (this will install the required packages locally).
3. Run `webapp_main.py`. This will open a local server (http://127.0.0.1:8050/) running peptimetric.


### Background
Peptimetric is a web app for visualizing and exploring differences in peptidomic datasets retrived from MS and MSMS. Peptimetric is souly written in Python and uses the Dash library from plotly for visualization. It is dependent on a local proteome database to fetch FASTA sequences.

### Example data 
The data used to illustrate the usage of Peptimetric was retrieved from a study by Van et al (Peptidomic Analysis of Urine from Youths with Early Type 1 Diabetes Reveals Novel Bioactivity of Uromodulin Peptides In Vitro, 2020). The data was searched using PEAKS Xpro 
and the resulting files are stored in the example-files directory.

## Contact information 
For any questions regarding Peptimetric, please contact:
Erik Hartman - erik.hartman@hotmail.com
Simon Mahdavi - simonmahdavi@msn.com
