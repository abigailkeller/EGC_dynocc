# Occupancy data objects

Run `Rscript code/data_prep/02_BayesianTrapData.R` from the project root.
This also builds the separate WSG objects. Both products retain June–September
and the existing study extent, with a continuous year axis from 2014 through
the latest retained individual-trap year. Unsampled cells are `NA`, not zero.

- `PresenceArrayBinary.rds`: individual-trap binary observations, arranged as
  `[site, year, replicate]`. Existing source and trap-type filters are retained;
  earlier available records, including the separate DNWR 2017 files, are included.
- `TrapTypeArraysNew.rds`: matching individual-trap one-hot indicators.
- `cpue_fukui.rds`: individual-trap Fukui CPUE by site and year.
- `WSGPresenceArrayBinary.rds`: a **separate** binary array in the same layout.
  A replicate is one raw sampling-event row; total catch greater than zero is 1.
  Events are ordered by date and original row within each site/year.
- `WSGTrapEffortArray.rds`: recorded total trap count in matching cells.
- `WSGTrapCountArrays.rds`: named Fukui/Shrimp/Minnow arrays in the same layout,
  currently `NA` because the raw file does not identify the trap mix. These are
  counts, not the one-hot indicators used for individual traps. Confirm the
  protocol before assigning three Minnow and three Fukui traps to an event;
  some events have total effort other than six.
- `WSGSamplingEvents.rds` / `.csv`: event records, dates, coordinates, source row,
  catch, effort, and replicate indices for auditing and future spatial mapping.

WSG site labels use `WSG::<SiteID>` to preserve the source's identifiers. The
year axes match, but the site axes and maximum replicate counts differ between
objects. Match sites explicitly before sharing latent occupancy states; array
row positions are not shared site identifiers.

This change prepares data only. Existing fitting scripts still use positional
year slices and connectivity inputs beginning in 2018. Update those selections,
site/zone mappings, historical covariates, and add the pooled-event likelihood
before fitting the extended joint model. WSG must use an event detection
probability, not an individual-trap probability. The older summary script 01
produces a different, aggregated product and is not an input to this workflow.
