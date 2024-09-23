# Early prediction of Drug-Induced Liver Injury

The DILI-Predictor predicts 10 features related to DILI toxicity including _in-vivo_ and _in-vitro_ and physicochemical parameters. It has been developed by the Broad Institute using the DILIst dataset (1020 compounds) from the FDA and achieved an accuracy balance of 70% on a test set of 255 compounds held out from the same dataset. The authors show how the model can correctly predict compounds that are not toxic in human despite being toxic in mice.

## Identifiers

* EOS model ID: `eos5gge`
* Slug: `dilipred`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Classification`
* Output: `Probability`
* Output Type: `Float`
* Output Shape: `List`
* Interpretation: Prediction of 10 DILI-related endpoints. The most important is the first, DILI. Threshold for DILI active is set at 0.16 by the authors.

## References

* [Publication](https://pubs.acs.org/doi/10.1021/acs.chemrestox.4c00015)
* [Source Code](https://github.com/Manas02/dili-pip)
* Ersilia contributor: [Zainab-ik](https://github.com/Zainab-ik)

## Ersilia model URLs
* [GitHub](https://github.com/ersilia-os/eos5gge)
* [AWS S3](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos5gge.zip)
* [DockerHub](https://hub.docker.com/r/ersiliaos/eos5gge) (AMD64, ARM64)

## Citation

If you use this model, please cite the [original authors](https://pubs.acs.org/doi/10.1021/acs.chemrestox.4c00015) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a None license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!