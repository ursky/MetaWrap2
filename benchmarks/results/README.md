# benchmarks/results/

The committed results history: one TSV per dataset, appended to by `metawrap2 benchmark run` and
`metawrap2 benchmark score`.

Committed rather than gitignored so that a changed default arrives in review next to the numbers
that justify it, and so a regression is a diff rather than a recollection.

Columns are `dataset  comparison  variant  sample  metric  value  unit  source`; see
[../README.md](../README.md) for why the table is long rather than wide, and why cost lives in it
alongside accuracy.

Nothing is here yet. The CAMI2 GI run needs its 95.4 GB fetched first:

```
metawrap2 benchmark fetch cami2_gi --dest benchmarks/data/cami2_gi --yes
```
