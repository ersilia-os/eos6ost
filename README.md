# REINVENT 4 LibInvent

REINVENT 4 LibInvent creates new molecules by appending R groups to a given input. If the input SMILES string contains specified attachment points, it is directly processed by LibInvent to generate new molecules. If no attachment points given, the model try to find potential attachment points, and iterates through different combinations of these points. It passes each combination to LibInvent to generate new molecules.

## Identifiers

* EOS model ID: `eos6ost`
* Slug: `reinvent4-libinvent`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Generative`
* Output: `Compound`
* Output Type: `String`
* Output Shape: `List`
* Interpretation: Model generates ~1000 similar molecules per input molecule.

## References

* [Publication](https://chemrxiv.org/engage/chemrxiv/article-details/65463cafc573f893f1cae33a)
* [Source Code](https://github.com/MolecularAI/REINVENT4)
* Ersilia contributor: [ankitskvmdam](https://github.com/ankitskvmdam)

## Ersilia model URLs
* [GitHub](https://github.com/ersilia-os/eos6ost)
* [AWS S3](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos6ost.zip)
* [DockerHub](https://hub.docker.com/r/ersiliaos/eos6ost) (AMD64)

## Citation

If you use this model, please cite the [original authors](https://chemrxiv.org/engage/chemrxiv/article-details/65463cafc573f893f1cae33a) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a Apache-2.0 license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!