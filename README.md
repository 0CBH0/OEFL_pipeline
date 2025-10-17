# OEFL_pipeline
## Analysis scripts for single-cell long-read sequencing of olfactory sensory neuron isoform profiles.


First run the shell script in pre-processing to handle the raw data.

Then oe_fl_fig_xxx.R could be uesd for figure generating.

- oe_fl_fig_01.R is used for figure 1 generating.

- oe_fl_fig_02.R is used for figure 2 generating.

- oe_fl_fig_03.R is used for figure 3 generating.

- oe_fl_fig_04.R is used for figure 4 generating.

- oe_fl_fig_S01.R is used for figure S1 generating.

- oe_fl_fig_S02.R is used for figure S2 generating.

- oe_fl_fig_S03.R is used for figure S3 generating.

- oe_fl_fig_S04.R is used for figure S4 generating.

```bash
├── pre-processing
│   ├── run_oefl.sh
│   └── sub-script
│       ├── isoquant_oe.py
│       ├── isoquant_osn.py
│       ├── calcReads_oe.R
│       └── calcReads_osn.R
├── oe_fl_fig_01.R
├── oe_fl_fig_02.R
├── oe_fl_fig_03.R
├── oe_fl_fig_04.R
├── oe_fl_fig_S01.R
├── oe_fl_fig_S02.R
├── oe_fl_fig_S03.R
└── oe_fl_fig_S04.R
```
