## *P. putida* genome-scale metabolic models

- The main model: *i*JN1463 ( [Nogales *et al.* (2020)](https://doi.org/10.1111/1462-2920.14843)
    - Curated *i*JN1463 (i.e. *i*JN1445, [Bergen et al. (2025)](https://www.biorxiv.org/content/10.1101/2025.02.17.638771v1.abstract)) used in [Puiggené et al. (2025)](https://www.cell.com/trends/biotechnology/fulltext/S0167-7799(25)00216-1).
    - Find more information on the curation of the model in [here](https://github.com/danbergen42/Co-Substrate-Free_Lignin_Valorisation).

## Errors/inaccuracy in the model *i*JN1463:

- Serine hydroxymethyltransferases (GlyA-I/II, GHMT2) are reversible in *P. putida* [Turlin *et al.* (2022)](https://doi.org/10.1016/j.ymben.2022.10.008)
- Homoserine dehydrogenase (HSDy) produces homoserine from aspartate-semialdehyde irreversibly in *P. putida* as we have proven experimentally [Puiggené et al. (2025)](https://www.cell.com/trends/biotechnology/fulltext/S0167-7799(25)00216-1).
- Formaldehyde dismutase (FALDM) reaction was deactivated given that it was never proven *in silico* nor *in vivo* that *P. putida* harbors a formaldehyde dismutase [Turlin and Puiggené *et al.* (2023)](https://doi.org/10.1128/msystems.00004-23)
- Given that formatotrophy via native PurU does not occur in *P. putida*, PurU (GARFT) was set irreversible [Puiggene *et al.* (2025)], [Turlin *et al.* (2022)](https://doi.org/10.1016/j.ymben.2022.10.008)
- Autotrophy based on CO<sub>2</sub> assimilation via the lipoamide-dependent reaction complexes (AKGDa, PDHa) was prevented by setting the reactions irreversible. 
- Finally, genes harbored in the TOL plasmid (reactions with pWW0 gene IDs) were deleted since *P. putida* KT2440 strain does not contain the plasmid.
- Nickel was removed from the biomass function, given that *P. putida* does not have that requirement for growth [Hartmans *et al.* (1989)](https://doi.org/10.1128/aem.55.11.2850-2855.1989)

## Other curations (upon necessity)

- Methanol oxidation, namely via PQQ-dependent methanol dehydrogenases, at high concentrations ([MeOH] > 125 mM) occurs in *P. putida* [Turlin and Puiggene *et al.* (2023)](https://doi.org/10.1128/msystems.00004-23), [Puiggené et al. (2025)](https://www.cell.com/trends/biotechnology/fulltext/S0167-7799(25)00216-1).
