# v0.10.3 recovery note

RSES-Onco v0.10.3 adds resilient OmniPath/DoRothEA acquisition after the validated
STRING stage completed but the OmniPath interaction service returned persistent
HTTP 502 responses.

See:

- `docs/DOROTHEA_RECOVERY_WORKFLOW.md`
- `docs/END_TO_END_ARTICLE_PROTOCOL_V0103_ADDENDUM.md`

The default workflow preserves a complete OmniPath outage as missing regulatory
coverage, never as zero evidence. Set `DOROTHEA_STRICT=1` when regulatory evidence
must be complete before downstream analysis.
