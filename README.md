# Early prediction of Drug-Induced Liver Injury

The DILI-Predictor predicts 10 features related to DILI toxicity including in-vivo and in-vitro and physicochemical parameters. It has been developed by the Broad Institute using the DILIst dataset (1020 compounds) from the FDA and achieved an accuracy balance of 70% on a test set of 255 compounds held out from the same dataset. The authors show how the model can correctly predict compounds that are not toxic in human despite being toxic in mice.

This model was incorporated on 2024-02-19.

## Information
### Identifiers
- **Ersilia Identifier:** `eos5gge`
- **Slug:** `dili-predictor`

### Domain
- **Task:** `Annotation`
- **Subtask:** `Activity prediction`
- **Biomedical Area:** `ADMET`
- **Target Organism:** `Homo sapiens`
- **Tags:** `Toxicity`, `Metabolism`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `10`
- **Output Consistency:** `Fixed`
- **Interpretation:** Prediction of 10 DILI-related endpoints. The most important is the first, DILI. Threshold for DILI active is set at 0.16 by the authors.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| dili | float | high | Classification score for drug-induced liver injury data |
| diverse_dili_c | float | high | Classification score for transient liver function abnormalities from Diverse DILI dataset C |
| besp | float | high | Classification score for bile salt export pump inhibition related to liver toxicity |
| mitotox | float | high | Classification score for mitochondrial toxicity leading to cellular energy failure |
| reactive_metabolite | float | high | Classification score for formation of reactive metabolites causing liver protein damage |
| human_hepatotoxicity | float | high | Classification score for hepatotoxic effects observed in human studies |
| animal_hepatotoxicity_a | float | high | Classification score for hepatic histopathologic effects from chronic animal studies |
| animal_hepatotoxicity_b | float | high | Classification score for hepatocellular hypertrophy effects in animals |
| preclinical_hepatotoxicity | float | high | Classification score for comprehensive preclinical hepatotoxicity assessment |
| diverse_dili_a | float | high | Classification score for diverse DILI dataset A encompassing heterogeneous liver toxicity data |


### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos5gge](https://hub.docker.com/r/ersiliaos/eos5gge)
- **Docker Architecture:** `AMD64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos5gge.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos5gge.zip)

### Resource Consumption
- **Model Size (Mb):** `67`
- **Environment Size (Mb):** `1707`
- **Image Size (Mb):** `923.16`

**Computational Performance (seconds):**
- 10 inputs: `36.85`
- 100 inputs: `962.01`
- 10000 inputs: `-1`

### References
- **Source Code**: [https://github.com/Manas02/dili-pip](https://github.com/Manas02/dili-pip)
- **Publication**: [https://pubs.acs.org/doi/10.1021/acs.chemrestox.4c00015](https://pubs.acs.org/doi/10.1021/acs.chemrestox.4c00015)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2024`
- **Ersilia Contributor:** [Zainab-ik](https://github.com/Zainab-ik)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [None](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos5gge
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos5gge
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
