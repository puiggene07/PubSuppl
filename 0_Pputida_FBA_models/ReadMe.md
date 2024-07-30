## *P. putida* genome-scale metabolic models

- The main model: *i*JN1463 ( [Nogales *et al.* (2020)](https://doi.org/10.1111/1462-2920.14843)
- Curated iJN1463 used in [Puiggene et al. (2024)]

## Errors/inaccuracy in the model *i*JN1463:

- Serine hydroxymethyltransferases (GlyA-I/II, GHMT2) are reversible in *P. putida* [Turlin *et al.* (2022)](https://doi.org/10.1016/j.ymben.2022.10.008)
- Homoserine dehydrogenase (HSDy) produces homoserine from aspartate-semialdehyde irreversibly in *P. putida* as we have proven experimentally [Puiggene *et al.* (2024)]
- Formaldehyde dismutase (FALDM) reaction was deactivated given that it was never proven *in silico* nor *in vivo* that *P. putida* harbors a formaldehyde dismutase [Turlin and Puiggene *et al.* (2023)](https://doi.org/10.1128/msystems.00004-23)
- Given that formatotrophy via native PurU does not occur in *P. putida*, PurU (GARFT) was set irreversible [Puiggene *et al.* (2024)], [Turlin *et al.* (2022)](https://doi.org/10.1016/j.ymben.2022.10.008)
- Autotrophy based on CO<sub>2</sub> assimilation via the lipoamide-dependent reaction complexes (AKGDa, PDHa) was prevented by setting the reactions irreversible. 
- Finally, genes harbored in the TOL plasmid (reactions with pWW0 gene IDs) were deleted since *P. putida* KT2440 strain does not contain the plasmid.
- Nickel was removed from the biomass function, given that *P. putida* does not have that requirement for growth [Hartmans *et al.* (1989)](https://doi.org/10.1128/aem.55.11.2850-2855.1989)

## Other curations (upon necessity)

- Methanol oxidation, namely via PQQ-dependent methanol dehydrogenases, at high concentrations ([MeOH] > 125 mM) occurs in *P. putida* [Turlin and Puiggene *et al.* (2023)](https://doi.org/10.1128/msystems.00004-23), [Puiggene *et al.* (2024)]
