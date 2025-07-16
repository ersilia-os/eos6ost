# REINVENT 4 LibInvent

REINVENT 4 LibInvent creates new molecules by appending R groups to a given input. If the input SMILES string contains specified attachment points, it is directly processed by LibInvent to generate new molecules. If no attachment points given, the model try to find potential attachment points, and iterates through different combinations of these points. It passes each combination to LibInvent to generate new molecules.

This model was incorporated on 2024-04-18.Last packaged on 2025-07-16.

## Information
### Identifiers
- **Ersilia Identifier:** `eos6ost`
- **Slug:** `reinvent4-libinvent`

### Domain
- **Task:** `Sampling`
- **Subtask:** `Generation`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Similarity`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `1000`
- **Output Consistency:** `Variable`
- **Interpretation:** Model generates up to 1000 similar molecules per input molecule.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| smi_000 | string |  | Generated compound index 0 using pre-trained LibInvent model |
| smi_001 | string |  | Generated compound index 1 using pre-trained LibInvent model |
| smi_002 | string |  | Generated compound index 2 using pre-trained LibInvent model |
| smi_003 | string |  | Generated compound index 3 using pre-trained LibInvent model |
| smi_004 | string |  | Generated compound index 4 using pre-trained LibInvent model |
| smi_005 | string |  | Generated compound index 5 using pre-trained LibInvent model |
| smi_006 | string |  | Generated compound index 6 using pre-trained LibInvent model |
| smi_007 | string |  | Generated compound index 7 using pre-trained LibInvent model |
| smi_008 | string |  | Generated compound index 8 using pre-trained LibInvent model |
| smi_009 | string |  | Generated compound index 9 using pre-trained LibInvent model |

_10 of 1000 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos6ost](https://hub.docker.com/r/ersiliaos/eos6ost)
- **Docker Architecture:** `AMD64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos6ost.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos6ost.zip)

### Resource Consumption
- **Model Size (Mb):** `176`
- **Environment Size (Mb):** `8349`
- **Image Size (Mb):** `8463.22`

**Computational Performance (seconds):**
- 10 inputs: `112.46`
- 100 inputs: `-1`
- 10000 inputs: `-1`

### References
- **Source Code**: [https://github.com/MolecularAI/REINVENT4](https://github.com/MolecularAI/REINVENT4)
- **Publication**: [https://chemrxiv.org/engage/chemrxiv/article-details/65463cafc573f893f1cae33a](https://chemrxiv.org/engage/chemrxiv/article-details/65463cafc573f893f1cae33a)
- **Publication Type:** `Preprint`
- **Publication Year:** `2023`
- **Ersilia Contributor:** [ankitskvmdam](https://github.com/ankitskvmdam)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [Apache-2.0](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos6ost
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos6ost
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
