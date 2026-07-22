# OmniPath outage interpretation

HTTP 5xx responses from the OmniPath web service are acquisition failures, not
negative biological evidence. RSES-Onco v0.10.3 therefore records the source as
unavailable and leaves regulatory-network score components missing unless strict
mode is enabled.
