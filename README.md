# VarEffectViz
VarEffectViz is a visualizer that presents variant effect interpretation for the BRCA1 gene using multiple sources of data, including allele count in UKBiobank, the CADD computational score, AlphaFold protein structure prediction and the results of saturation genome editing experiments.

The SGE function score in the vizualisation is the SGE function score is the normalized negative function score obtained using min-max normalization on function score from Findlay *et al.*, 2018. The cumulative score is computed as follow:

$cumulative\_{score} = \frac{1}{3} \left(\frac{1}{AC} + SGE_{score} + \frac{CADD_{score}}{50}\right)$

where $AC$ is allele count, $SGE_{score}$ is the normalized negative function score obtained using min-max normalization on function score from Findlay *et al.*, 2018, and cadd_score is the CADD score divided by the maximum cadd score in this region which is 50. Note that 1/3 is a constant factor used to normalize the sum of the three scores.


## Getting Started

### Online Version
An online version of the app is also available at https://vareffectviz.onrender.com/, but please note that it may be slower due to the limitations of the free deployment. It is recommended to run the app locally for optimal performance.

### Prerequisites

Before running the app, make sure you have the following installed:

- Python 3.10 or later
- Dash 2.7.0 or later
- Pandas 1.5.0 or later
- Plotly 5.14.0 or later

### Installation

1. Clone this repository to your local machine.
2. Navigate to the project directory in your terminal.
3. Install the required dependencies using pip:

```pip install -r requirements.txt```

### Running the App

1. Navigate to the project directory in your terminal.
2. Run the app by executing `app.py`:

```python app.py```

3. Open a web browser and navigate to `http://localhost:8050` to view the app.

## Usage

The app allows users to interact with and visualize BRCA1 variants. Users can select different options from dropdown menus and update the visualization by selecting subset of variants or differents annotations.

## Contributing

Contributions are welcome! If you find a bug or have a feature request, please open an issue on the GitHub repository.
