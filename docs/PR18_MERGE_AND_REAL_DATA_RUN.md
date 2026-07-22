# PR #18 integration and real-data execution

PR #18 is merged into `main` at commit `1ed801dff5a57c6e25b3f13760aed0fc4f794ce0`.

Update a local WSL checkout:

```bash
NEW="/mnt/c/Users/Microsoft/Desktop/SynLeth-RSES-Onco-v010"
cd "$NEW" || exit 1

git fetch origin
git checkout main
git pull --ff-only origin main

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rses-onco
python -m pip install -e .

printf 'HEAD: '
git rev-parse HEAD
grep '^version' pyproject.toml
```

The expected version is `0.10.7`. The real-data rerun must recalculate the ranking with `scoring_semantics_version=eligibility-aware-v1` before publication assets are accepted.
