# Reuse existing results

Edit server_config.sh or pass `--config FILE`. Defaults use testdata/ and Python from the active environment.

```bash
bash usage_scenarios/run_case.sh reuse_all
```

| Case | Additional inputs |
|---|---|
| auto | ARG database |
| orfs | ARG database, ORF GFF and proteins |
| arg_hits | ORF GFF, proteins and raw ARG hits |
| ani | ARG database and FastANI results |
| blastn | ARG database and BLASTn results |
| reuse_all | BLASTn, FastANI, ORF GFF/proteins and ARG hits |
| reuse_tsv | Same as reuse_all, using ORF TSV |
| filtered_hits | Same as reuse_all, using prefiltered ARG hits |
| hgt_only | MAGs and taxonomy |

testdata supports reuse_all, reuse_tsv and filtered_hits. Set real input paths for the other cases. Output: scenario_results/CASE/. Use `--dry-run` to print the command.

Export inputs from a completed run:

```bash
python3 usage_scenarios/prepare_existing_inputs.py --run-dir /path/to/run -o reusable_inputs
```

Set the config to the exported paths. [Input formats](../docs/INPUT_FORMATS.md).
