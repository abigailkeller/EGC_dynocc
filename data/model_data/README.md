# Occupancy data objects

Run `Rscript code/data_prep/02_BayesianTrapData.R` from the project root.
This also builds the separate WSG objects. Both products retain June–September
and use a continuous year axis from 2014 through
the latest retained individual-trap year. Unsampled cells are `NA`, not zero.
Files with the original names retain the original study-extent filter. Files
ending in `_expanded` use the same processing but include every site with usable
coordinates, including sites outside the original study extent.
`DFO::possibleGPSerror` is explicitly excluded because its source label flags
the coordinates as unreliable.

- `PresenceArrayBinary.rds`: individual-trap binary observations, arranged as
  `[site, year, replicate]`. Existing source and trap-type filters are retained;
  earlier available records, including the separate DNWR 2017 files, are included.
- `TrapTypeArraysNew.rds`: matching individual-trap one-hot indicators.
- `cpue_fukui.rds`: individual-trap Fukui CPUE by site and year.
- `PresenceArrayBinary_expanded.rds`, `TrapTypeArraysNew_expanded.rds`, and
  `cpue_fukui_expanded.rds`: expanded counterparts with identical definitions.
- `PresenceArrayBinary_WSG.rds`: a **separate** binary array in the same layout.
  A replicate is a calendar month, with axis labels `6`, `7`, `8`, `9`
  (June–September). Each site/year/month is 1 if any retained record has positive
  catch, 0 if sampled with no detections, and `NA` if unsampled. Multiple records
  in the same month are pooled into one observation.
  The usual protocol is three Minnow plus three Fukui traps; separate effort
  and trap-count arrays are no longer produced. Recorded effort remains in the
  audit table because some events report fewer than six traps.
- `PresenceArrayBinary_WSG_expanded.rds`: expanded WSG monthly presence data.
- `SamplingEvents_WSG.rds` / `.csv`: event records, dates, coordinates, source row,
  catch, effort, month, and replicate indices (1–4 for June–September) for
  auditing and future spatial mapping. Rows remain raw events, so multiple rows
  can map to the same monthly array cell.
- `site_zone_map.csv`: zone assignments for the individual-trap site axis.
- `site_zone_map_WSG.csv`: zone assignments for the separate WSG site axis,
  including both the `WSG::<SiteID>` key and readable site name. Run
  `Rscript code/data_prep/03_map_model_sites_to_zones.R` to rebuild both tables.
- `site_zone_map_expanded.csv` and `site_zone_map_WSG_expanded.csv`: matching
  zone tables for the expanded array site axes. Sites beyond the zone polygons
  retain `NA` zone assignments rather than being assigned arbitrarily.
- `site_coordinates_expanded.csv`: all retained individual-trap and WSG sites,
  their longitude (`x`) and latitude (`y`), and an `in_study_extent` audit flag.
- `site_coordinates_outside_extent.csv`: the outside-only subset of that list.
- `expanded_data_validation.csv` and `expanded_data_flags.csv`: reproducible QA
  results from `code/data_prep/validate_expanded_model_data.R`.

WSG site labels use `WSG::<SiteID>` to preserve the source's identifiers. The
year axes match, but the site axes and maximum replicate counts differ between
objects. Match sites explicitly before sharing latent occupancy states; array
row positions are not shared site identifiers.

This change prepares data only. Existing fitting scripts still use positional
year slices and connectivity inputs beginning in 2018. Update those selections,
site/zone mappings, historical covariates, and add the pooled-month likelihood
before fitting the extended joint model. WSG must use a monthly detection
probability, not an individual-trap probability. The older summary script 01
produces a different, aggregated product and is not an input to this workflow.
