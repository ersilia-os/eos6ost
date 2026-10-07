# REINVENT 4 LibInvent

Grows new R-groups onto a scaffold instead of editing a whole molecule. An input that already carries attachment points is decorated directly; otherwise it is reduced to its Murcko scaffold, discarding substituents, stereochemistry and functional groups, and attachment points are enumerated on free scaffold carbons in combinations of one to three. The RNN generator, shipped with REINVENT 4 and pre-trained on ChEMBL 27, is fed a randomised SMILES each run, so the 100 structures returned vary between runs and need not resemble the query.

This model was incorporated on 2024-04-18.Last packaged on 2026-09-29.

## Information
### Identifiers
- **Ersilia Identifier:** `eos6ost`
- **Slug:** `reinvent4-libinvent`

### Domain
- **Task:** `Sampling`
- **Subtask:** `Generation`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Compound generation`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `100`
- **Output Consistency:** `Variable`
- **Interpretation:** 100 molecules built by decorating the input's Murcko scaffold with R-groups, so its own substituents are lost.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| smi_00 | string |  | Generated compound index 0 using pre-trained LibInvent model |
| smi_01 | string |  | Generated compound index 1 using pre-trained LibInvent model |
| smi_02 | string |  | Generated compound index 2 using pre-trained LibInvent model |
| smi_03 | string |  | Generated compound index 3 using pre-trained LibInvent model |
| smi_04 | string |  | Generated compound index 4 using pre-trained LibInvent model |
| smi_05 | string |  | Generated compound index 5 using pre-trained LibInvent model |
| smi_06 | string |  | Generated compound index 6 using pre-trained LibInvent model |
| smi_07 | string |  | Generated compound index 7 using pre-trained LibInvent model |
| smi_08 | string |  | Generated compound index 8 using pre-trained LibInvent model |
| smi_09 | string |  | Generated compound index 9 using pre-trained LibInvent model |

_10 of 100 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos6ost](https://hub.docker.com/r/ersiliaos/eos6ost)
- **Docker Architecture:** `AMD64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos6ost.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos6ost.zip)

### Resource Consumption
- **Model Size (Mb):** `175`
- **Environment Size (Mb):** `2367`
- **Image Size (Mb):** `2521.33`

**Computational Performance (seconds):**
- 10 inputs: `37.98`
- 100 inputs: `1024.41`
- 10000 inputs: `-1`

### References
- **Source Code**: [https://github.com/MolecularAI/REINVENT4](https://github.com/MolecularAI/REINVENT4)
- **Publication**: [https://doi.org/10.1186/s13321-024-00812-5](https://doi.org/10.1186/s13321-024-00812-5)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2024`
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
