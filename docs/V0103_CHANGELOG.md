# RSES-Onco v0.10.3

- Reuse complete STRING aggregate/status output without repeating 513 partner calls.
- Add resilient DoRothEA acquisition with validated cache, local TSV and multiple
  official OmniPath service addresses.
- Record every DoRothEA attempt and preserve persistent service outages as explicit
  missing regulatory coverage by default.
- Add strict DoRothEA mode and local-file support through environment variables.
- Add regression tests for cache reuse, local input, outage missingness and strict
  failure.
- Update the end-to-end recovery protocol and software version metadata.
