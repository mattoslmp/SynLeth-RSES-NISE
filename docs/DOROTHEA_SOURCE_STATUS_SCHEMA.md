# DoRothEA source-status schema

`data/raw/human_functional_evidence/omnipath_dorothea_status.json` records:

- `status`: `downloaded`, `cache`, `local_file` or `unavailable`;
- `available`: whether a validated regulatory table was available;
- `rows`: number of validated TF-target records;
- `source`: endpoint or local file used;
- `attempts`: each endpoint, organism parameter, result and error;
- `scientific_interpretation`: explicit missingness statement when unavailable.

Pair evidence contains `regulatory_source_available` and
`regulatory_source_status`. Missing DoRothEA coverage leaves
`component_regulatory_network`, `regulator_jaccard` and `shared_regulators`
unobserved; it is not encoded as negative evidence.
