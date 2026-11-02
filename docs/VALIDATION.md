# Tests

```bash
python3 -m unittest discover -s tests -v
```

37 tests cover fragment selection, coordinates, ARG annotation, optional-input reuse and file handling. External-tool calls use mocks. The synthetic example checks the main table and both FASTA exports against expected files.

Real-data runs using existing Conda environments produced matching counts for automatic annotation, imported ORF/ARG hits and fully reused inputs: 18 fragment pairs, 9 pairs with a complete ARG on at least one endpoint, and 15 overlapping ARG ORFs. A fresh environment.yml installation is a separate check.
