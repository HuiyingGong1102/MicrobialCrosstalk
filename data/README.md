# Input data contract

The local delivery contains copies of the three supplied workbooks. `.xlsx` files
are ignored by Git; either distribute them separately or deliberately change the
ignore rule when preparing the public repository.

`k40.RData` is a byte-identical copy of the supplied best clustering result from
`D:/Doctoral/CAF2/f_b2/AR-cluster/k40.RData`. It is also ignored by Git. It contains
`return_object$clustered_data`, a 563-by-31 matrix with OTU row names and a named
`cluster` column, assigning every retained OTU to one of 40 modules.

| Workbook | Rows in supplied file | Required fields |
|---|---:|---|
| fungi_otu.xlsx | 210 | sample, Treatment, Month, OTU count columns |
| bacteria_out.xlsx | 210 | sample, Treatment, Month, OTU count columns |
| soil_properties.xlsx | 210 | sample, Treatment, Month, CMF, NMF, PMF |

Sample IDs must be unique and identical across files. Rows are matched by `sample`,
not position. Treatment and month metadata must agree. Count columns begin with
`OTU`, are numeric, finite, and nonnegative. Additional metadata is permitted.

The supplied six treatments are CKN, CKR, T30N, T30R, T45N, and T45R. The default
analysis uses April through September (180 rows; 30 per treatment) and excludes
December by name, reproducing the original study subset without positional row
deletion. SMF is not present in the supplied workbook.
